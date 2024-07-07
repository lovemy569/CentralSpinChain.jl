"""
Creates a function to define the initial state of the spin chain

Parameters:
- N: Length of the spin chain (excluding the central spin)
- center_spin_initial_state: Initial state of the central spin
- chain_initial_state: Initial state of the spin chain

Returns:
- A function that returns the spin state at position `n`.
"""
function create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    return n -> n == N + 1 ? center_spin_initial_state : (
        chain_initial_state == "Néel_state" ? (isodd(n) ? "Up" : "Dn") :
        chain_initial_state == "all_up" ? "Up" :
        chain_initial_state == "all_dn" ? "Dn" :
        chain_initial_state == "all_x+" ? "X+" :
        chain_initial_state == "all_x-" ? "X-" : error("Invalid chain initial state")
    )
end

"""
Creates a list of gates for the spin chain

Parameters:
- N: Length of the spin chain
- s: Spin sites in the chain
- Jx: Exchange coefficient in the X direction
- Jy: Exchange coefficient in the Y direction
- Jz: Exchange coefficient in the Z direction
- tau: Time step
- F: Linear field coefficient
- W: Range for the random field coefficients
- specified_sites: Pairs of spin sites to be excluded from gate operations

Returns:
- An array of ITensor gates for the time evolution of the spin chain.
"""
function create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, specified_sites)
    num_gates = N - length(specified_sites) + N  # Calculate total number of gates
    gates = Vector{ITensor}(undef, num_gates)  # Pre-allocate the gates array
    gate_idx = 1  # Initialize gate index
    
    # Precompute all periodic boundary conditions
    next_site = [mod(j, N) + 1 for j in 1:N]

    # Create gate operations for spin chain coupling
    for j in 1:N
        s1 = s[j]
        s2 = s[next_site[j]]  # Periodic boundary condition: couple with the next spin
        if (j, next_site[j]) in specified_sites || (next_site[j], j) in specified_sites
            continue  # Skip specified pairs of spin sites
        end
        # Define Heisenberg Hamiltonian term
        hj = Jx * op("Sx", s1) * op("Sx", s2) + 
             Jy * op("Sy", s1) * op("Sy", s2) +
             Jz * op("Sz", s1) * op("Sz", s2)
        # Create time evolution operator
        gates[gate_idx] = exp(-im * tau * hj)
        gate_idx += 1  # Update gate index
    end

    # Generate uniformly distributed random field
    h = rand([-W, W], N)
    for j in 1:N
        # Linear field and random field term in the Z direction
        fz = F * j * op("Sz", s[j]) + h[j] * op("Sz", s[j])
        # Create time evolution operator
        gates[gate_idx] = exp(-im * tau * fz)
        gate_idx += 1  # Update gate index
    end

    return gates  # Return the array of gate operations
end

"""
Add central spin coupling terms to the list of gate operations

Parameters:
- gates: Array of gate operations for the spin chain
- N: Length of the spin chain
- s: Spin sites in the chain
- center_spin: Central spin site
- gamma: Central spin coupling coefficient
- Cx: Coupling coefficient for X direction between central spin and chain spins
- Cy: Coupling coefficient for Y direction between central spin and chain spins
- Cz: Coupling coefficient for Z direction between central spin and chain spins
- tau: Time step

This function iterates over each spin site in the chain, calculates the coupling term with the central spin, 
and adds the corresponding time evolution operator to the gate operations array.
"""
function add_center_spin_coupling!(gates, N, s, center_spin, gamma, Cx, Cy, Cz, tau)
    # Pre-compute constants
    gamma_over_N = gamma / N

    # Iterate over each spin site in the chain
    for j in 1:N
        # Calculate the coupling Hamiltonian term between the central spin and the chain spin
        hjc = gamma_over_N * (
            Cx * (op("Sx", s[j]) * op("Sx", center_spin)) +
            Cy * (op("Sy", s[j]) * op("Sy", center_spin)) +
            Cz * (op("Sz", s[j]) * op("Sz", center_spin))
        )
        # Create the time evolution operator
        Gjc = exp(-im * tau * hjc)
        # Add the time evolution operator to the gate operations array
        push!(gates, Gjc)
    end
end

"""
Calculate dynamic coupling coefficients

Parameters:
- t: Current time
- Jx: Exchange coefficient in the X direction
- Jy: Exchange coefficient in the Y direction
- Jz: Exchange coefficient in the Z direction
- period: Period of the dynamic coupling

Returns:
- Returns the dynamically adjusted exchange coefficients for the X, Y, and Z directions.
"""
function coupling_strength(t, Jx, Jy, Jz, period)
    half_period = period / 2  # Precompute half of the period
    in_first_half = mod(t, period) < half_period  # Check if in the first half of the period
    
    Jx_dynamic = in_first_half ? Jx : 0.0
    Jy_dynamic = in_first_half ? Jy : 0.0
    Jz_dynamic = in_first_half ? Jz : 0.0
    
    return Jx_dynamic, Jy_dynamic, Jz_dynamic
end

"""
Calculate dynamic coupling coefficients for the central spin

Parameters:
- t: Current time
- Cx: Coupling coefficient for the X direction between the central spin and chain spins
- Cy: Coupling coefficient for the Y direction between the central spin and chain spins
- Cz: Coupling coefficient for the Z direction between the central spin and chain spins
- center_spin_period: Period of the dynamic coupling for the central spin

Returns:
- Returns the dynamically adjusted coupling coefficients for the X, Y, and Z directions.
"""
function center_spin_coupling_strength(t, Cx, Cy, Cz, center_spin_period)
    half_period = center_spin_period / 2  # Precompute half of the period
    in_first_half = mod(t, center_spin_period) < half_period  # Check if in the first half of the period
    
    Cx_dynamic = in_first_half ? Cx : 0.0
    Cy_dynamic = in_first_half ? Cy : 0.0
    Cz_dynamic = in_first_half ? Cz : 0.0
    
    return Cx_dynamic, Cy_dynamic, Cz_dynamic
end

"""
Apply dynamic coupling to the quantum state

Parameters:
- psi: Current quantum state
- s: Array of spin sites
- specified_sites: Specified pairs of spin sites
- t: Current time
- gamma: Central spin coupling coefficient
- Jx: Exchange coefficient in the X direction
- Jy: Exchange coefficient in the Y direction
- Jz: Exchange coefficient in the Z direction
- Cx: Coupling coefficient for the X direction between central spin and chain spins
- Cy: Coupling coefficient for the Y direction between central spin and chain spins
- Cz: Coupling coefficient for the Z direction between central spin and chain spins
- tau: Time step
- period: Period of the dynamic coupling
- center_spin_period: Period of the dynamic coupling for the central spin

Returns:
- Returns the quantum state after applying dynamic coupling.
"""
function apply_dynamic_coupling(psi, s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, tau, period, center_spin_period)
    # Precompute dynamic coupling coefficients
    Jx_dynamic, Jy_dynamic, Jz_dynamic = coupling_strength(t, Jx, Jy, Jz, period)
    Cx_dynamic, Cy_dynamic, Cz_dynamic = center_spin_coupling_strength(t, Cx, Cy, Cz, center_spin_period)
    center_spin = s[end]

    # Handle dynamic coupling on the chain
    for (j1, j2) in specified_sites
        s1, s2 = s[j1], s[j2]
        hj = Jx_dynamic * op("Sx", s1) * op("Sx", s2) + Jy_dynamic * op("Sy", s1) * op("Sy", s2) + Jz_dynamic * op("Sz", s1) * op("Sz", s2)
        G_dynamic = exp(-im * tau * hj)
        psi = apply(G_dynamic, psi; cutoff=1e-9)
    end

    # Handle dynamic coupling with the central spin
    gamma_over_len_s = gamma / length(s)
    for j in 1:length(s) - 1
        sj = s[j]
        hjc = gamma_over_len_s * (Cx_dynamic * op("Sx", sj) * op("Sx", center_spin) +
                                  Cy_dynamic * op("Sy", sj) * op("Sy", center_spin) +
                                  Cz_dynamic * op("Sz", sj) * op("Sz", center_spin))
        Gjc_dynamic = exp(-im * tau * hjc)
        psi = apply(Gjc_dynamic, psi; cutoff=1e-9)
    end

    return psi
end

"""
Calculate the entanglement entropy of the given quantum state psi at site n

Parameters:
- psi: Current quantum state
- n: Index of the site to calculate the entanglement entropy

Returns:
- Returns the entanglement entropy value.
"""
function measure_EE(psi, n)
    orthogonalize!(psi, n)  # Orthogonalize psi at site n
    if n == 1
        _, S = svd(psi[n], (siteinds(psi, n),))  # Perform SVD and get singular values
    else
        _, S = svd(psi[n], (linkind(psi, n - 1), siteinds(psi, n)); cutoff=1e-18, alg="qr_iteration")  # Perform SVD and get singular values
    end
    SvN = 0.00
    for n in 1:dim(S, 1)
        p = S[n, n]^2
        SvN -= p * log(p)
    end
    return SvN  # Return entanglement entropy
end

"""
Perform time evolution

Parameters:
- N: Length of the spin chain
- ttotal: Total time
- tau: Time step
- psi: Initial quantum state
- gates: Array of gate operations
- s: Array of spin sites
- specified_sites: Specified pairs of spin sites
- Jx: Exchange coefficient in the X direction
- Jy: Exchange coefficient in the Y direction
- Jz: Exchange coefficient in the Z direction
- Cx: Coupling coefficient for the X direction between central spin and chain spins
- Cy: Coupling coefficient for the Y direction between central spin and chain spins
- Cz: Coupling coefficient for the Z direction between central spin and chain spins
- gamma: Central spin coupling coefficient
- dynamic_period: Period of the dynamic coupling
- center_spin_period: Period of the dynamic coupling for the central spin

Returns:
- Returns arrays of time, Sz, Sy, Sx, entropy, CSz, and CSz values.
"""
function time_evolution(N, ttotal, tau, psi, gates, s, specified_sites, Jx, Jy, Jz, Cx, Cy, Cz, gamma, dynamic_period, center_spin_period)
    Sz_array = Float64[]
    Sy_array = Float64[]
    Sx_array = Float64[]
    Entropy_array = Float64[]
    t_array = Float64[]

    num_steps = Int(ttotal / tau) + 1
    CSz_array = zeros(Float64, num_steps, N)
    CSz = zeros(Float64, num_steps, N)

    step_index = 1
    @showprogress "Time Evolution" for t in 0:tau:ttotal
        Sz = 2 * expect(complex(psi), "Sz"; sites = N + 1)
        Sy = 2 * expect(complex(psi), "Sy"; sites = N + 1)
        Sx = 2 * expect(complex(psi), "Sx"; sites = N + 1)
        push!(t_array, t)
        push!(Sz_array, Sz)
        push!(Sy_array, Sy)
        push!(Sx_array, Sx)

        for j in 1:N
            sz_j = expect(psi, "Sz"; sites = j)
            CSz[step_index, j] = 2 * sz_j
            CSz_array[step_index, j] = (-1) ^ (j + 1) * (sz_j + 0.5)
        end
        step_index += 1

        SvN = measure_EE(psi, N)
        push!(Entropy_array, SvN)

        psi = apply_dynamic_coupling(psi, s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, tau, dynamic_period, center_spin_period)
        psi = apply(gates, psi; cutoff=1e-9)
        normalize!(psi)
    end

    return t_array, Sz_array, Sy_array, Sx_array, Entropy_array, CSz_array, CSz
end

"""
Create and display plots

Parameters:
- t_array: Array of time values
- Sz_array: Array of Sz expectation values
- Sy_array: Array of Sy expectation values
- Sx_array: Array of Sx expectation values
- CSz_array: Array of central spin Sz values
- CSz: Central spin Sz values
- Entropy_array: Array of entanglement entropy values

This function creates three plots: spin plot, imbalance plot, and entropy plot, and combines them for display.
"""
function create_plots(t_array, Sz_array, Sy_array, Sx_array, CSz_array, CSz, Entropy_array)
    # Create spin plot
    spin_plot = plot(
        t_array, Sz_array, 
        label=L"$\langle S_z \rangle$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Central Spin"
    )
    plot!(spin_plot, t_array, Sy_array, label=L"$\langle S_y \rangle$", legend=:best, color=:green, line=(:dash, 1.5))
    plot!(spin_plot, t_array, Sx_array, label=L"$\langle S_x \rangle$", legend=:best, color=:red, line=(:dot, 1.5))

    # Calculate imbalance
    imbalance_numerator = sum(CSz_array, dims=2)
    imbalance_denominator = sum(CSz ./ 2 .+ 0.5, dims=2)
    imbalance = imbalance_numerator ./ imbalance_denominator
    imbalance_plot = plot(
        t_array, imbalance, 
        label=L"$I=(S_{z,\uparrow}^{e}-S_{z,\uparrow}^{o})/(S_{z,\uparrow}^{o}+S_{z,\uparrow}^{e})$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Imbalance"
    )

    # Create entropy plot
    entropy_plot = plot(
        t_array, Entropy_array, 
        label=L"$S_{central}=-\sum_np_n\log p_n$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Entropy"
    )

    # Combine and display all plots
    combined_plot = plot(spin_plot, imbalance_plot, entropy_plot, layout=(1, 3), size=(1080, 320))
    display(combined_plot)
end

"""
Run simulation

Parameters:
- N: Length of the spin chain
- Jz: Exchange coefficient in the Z direction
- Jy: Exchange coefficient in the Y direction
- Jx: Exchange coefficient in the X direction
- Cz: Coupling coefficient for the Z direction between central spin and chain spins
- Cy: Coupling coefficient for the Y direction between central spin and chain spins
- Cx: Coupling coefficient for the X direction between central spin and chain spins
- F: Linear field coefficient
- W: Random field coefficient
- gamma: Central spin coupling coefficient
- tau: Time step
- ttotal: Total time
- cutoff: Cutoff precision
- center_spin_initial_state: Initial state of the central spin
- chain_initial_state: Initial state of the spin chain
- specified_sites: Specified pairs of spin sites
- dynamic_period: Period of the dynamic coupling
- center_spin_period: Period of the dynamic coupling for the central spin

This function initializes the spin chain, creates gate operations, performs time evolution, and plots the results.
"""
function run_simulation(N, Jz, Jy, Jx, Cz, Cy, Cx, F, W, gamma, tau, ttotal, cutoff, center_spin_initial_state, chain_initial_state, specified_sites, dynamic_period, center_spin_period)
    # Initialize spin sites
    s = siteinds("S=1/2", N + 1)
    initial_state_function = create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    psi = MPS(s, initial_state_function)

    # Create gate operations
    gates = create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, specified_sites)
    add_center_spin_coupling!(gates, N, s, s[N+1], gamma, Cx, Cy, Cz, tau)

    # Perform time evolution
    t_array, Sz_array, Sy_array, Sx_array, Entropy_array, CSz_array, CSz = time_evolution(N, ttotal, tau, psi, gates, s, specified_sites, Jx, Jy, Jz, Cx, Cy, Cz, gamma, dynamic_period, center_spin_period)

    # Plot results
    create_plots(t_array, Sz_array, Sy_array, Sx_array, CSz_array, CSz, Entropy_array)
end



