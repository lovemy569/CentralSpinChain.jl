using ITensors
using Plots
using HDF5
using Dates

"""
Creates a function to define the initial state of the spin chain

Parameters:
- N: Length of the spin chain (excluding the central spin)
- center_spin_initial_state: Initial state of the central spin
- chain_initial_state: Initial state of the spin chain

Returns:
- A function that returns the spin state at position n.
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
Create a list of gates for the spin chain

Parameters:
- N: Length of the spin chain
- s: Spin sites in the chain
- Jx: Exchange coefficient in the X direction
- Jy: Exchange coefficient in the Y direction
- Jz: Exchange coefficient in the Z direction
- tau: Time step
- F: Linear field coefficient
- W: Range for the random field coefficients
- α: Coefficient for the non-uniform field
- specified_sites: Pairs of spin sites to be excluded from gate operations
- periodic: Whether to use periodic boundary conditions, true for periodic, false for open boundary

Returns:
- An array of ITensor gates for the time evolution of the spin chain.
"""
function create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, α, specified_sites, periodic::Bool)
    center_spin_site = N + 1
    count = 0
    for (i, j) in specified_sites
        if i == center_spin_site || j == center_spin_site
            count += 1
        end
    end

    num_gates = (periodic ? N : N - 1) - length(specified_sites) + N + count # Calculate total number of gates
    gates = Vector{ITensor}(undef, num_gates)  # Pre-allocate the gates array
    gate_idx = 1  # Initialize gate index

    # Precompute the next site indices
    next_site = periodic ? [mod(j, N) + 1 for j in 1:N] : [j + 1 for j in 1:N-1]

    # Create gate operations for spin chain coupling
    for j in 1:(periodic ? N : N-1)  # For open boundary, exclude the last site
        s1 = s[j]
        s2 = s[next_site[j]]  # Next spin site
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
    L = N  # Chain length
    for j in 1:N
        # Linear field and random field term in the Z direction, including non-uniform term
        Wj = -F * j + α * j^2 / (L)^2
        # Wj = (j % 2 == 1) ? 0.0 : -5.0
        fz = Wj * op("Sz", s[j]) + h[j] * op("Sz", s[j])
        # Create time evolution operator
        gates[gate_idx] = exp(-im * tau * fz)
        gate_idx += 1  # Update gate index
    end

    return gates  # Return the array of gate operations
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
    if period == 0
        return Jx, Jy, Jz
    end 

    half_period = period / 2  # Precompute half of the period
    in_first_half = mod(t, period) < half_period  # Check if in the first half of the period
    
    Jx_dynamic = in_first_half ? 0.0 : Jx
    Jy_dynamic = in_first_half ? 0.0 : Jy
    Jz_dynamic = in_first_half ? 0.0 : Jz
    
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
    if center_spin_period == 0
        return Cx, Cy, Cz
        # return 0,0,0
    end

    half_period = center_spin_period / 2  # Precompute half of the period
    in_first_half = mod(t, center_spin_period) < half_period  # Check if in the first half of the period
    
    Cx_dynamic = in_first_half ? 0.0 : Cx
    Cy_dynamic = in_first_half ? 0.0 : Cy
    Cz_dynamic = in_first_half ? 0.0 : Cz
    
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
- Ω: Coefficient for the additional term in the Z direction for the central spin

Returns:
- Returns the quantum state after applying dynamic coupling.
"""
function apply_dynamic_coupling(psi, s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, Ω, tau, period, center_spin_period)
    # 计算动态耦合系数
    Jx_dynamic, Jy_dynamic, Jz_dynamic = coupling_strength(t, Jx, Jy, Jz, period)
    Cx_dynamic, Cy_dynamic, Cz_dynamic = center_spin_coupling_strength(t, Cx, Cy, Cz, center_spin_period)
    center_spin = s[end]  # 中心自旋的站点索引

    # 处理自旋链上的动态耦合
    for (j1, j2) in specified_sites
        if j1 != length(s) && j2 != length(s)
            s1, s2 = s[j1], s[j2]
            hj = Jx_dynamic * op("Sx", s1) * op("Sx", s2) + Jy_dynamic * op("Sy", s1) * op("Sy", s2) + Jz_dynamic * op("Sz", s1) * op("Sz", s2)
            G_dynamic = exp(-im * tau * hj)
            psi = apply(G_dynamic, psi; cutoff=1e-9)
        end
    end

    # 处理中心自旋的动态耦合
    # gamma_over_len_s = gamma / (length(s) - 1)
    gamma_over_len_s = 1
    specified_center_spin_sites = Set((j1 == length(s) ? j2 : j1) for (j1, j2) in specified_sites if j1 == length(s) || j2 == length(s))
    
    for j in specified_center_spin_sites
        sj = s[j]
        hjc_dynamic = gamma_over_len_s * (Cx_dynamic * op("Sx", sj) * op("Sx", center_spin) +
                                          Cy_dynamic * op("Sy", sj) * op("Sy", center_spin) +
                                          Cz_dynamic * op("Sz", sj) * op("Sz", center_spin))
        Gjc_dynamic = exp(-im * tau * hjc_dynamic)
        psi = apply(Gjc_dynamic, psi; cutoff=1e-9)
    end

    # 处理中心自旋的静态耦合
    for j in 1:length(s) - 1
        if j in specified_center_spin_sites
            continue
        end
        sj = s[j]
        hjc_default = gamma_over_len_s * (Cx * op("Sx", sj) * op("Sx", center_spin) +
                                          Cy * op("Sy", sj) * op("Sy", center_spin) +
                                          Cz * op("Sz", sj) * op("Sz", center_spin))
        Gjc_default = exp(-im * tau * hjc_default)
        psi = apply(Gjc_default, psi; cutoff=1e-9)
    end

    # 中心自旋的额外项 \Omega * σz
    omega_z_term = Ω * op("Sz", center_spin)
    G_omega_z = exp(-im * tau * omega_z_term)
    psi = apply(G_omega_z, psi; cutoff=1e-9)

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
function time_evolution(N, ttotal, tau, psi, gates, s, specified_sites, Jx, Jy, Jz, Cx, Cy, Cz, Ω, gamma, dynamic_period, center_spin_period)
    Sz_array = Float64[]
    Sy_array = Float64[]
    Sx_array = Float64[]
    Entropy_array = Float64[]
    Center_Entropy_array = Float64[]
    t_array = Float64[]
    tcorr_array = ComplexF64[]
    zz_corr_array = Vector{Vector{Float64}}()
    sp_corr_array = Vector{Vector{Float64}}()
    sm_corr_array = Vector{Vector{Float64}}()

    psi_z = deepcopy(psi)
    opp = op("Sz",s[N+1])
    A_z = opp * psi_z[N+1]
    noprime!(A_z)
    psi_z[N+1] = A_z
    os = OpSum()
    os += "Sz", N+1
    H = MPO(os, s)

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
            sz_j = expect(complex(psi), "Sz"; sites = j)
            CSz[step_index, j] = 2 * sz_j
            CSz_array[step_index, j] = (-1) ^ (j + 1) * (sz_j + 0.5)
        end
        step_index += 1

        SvN = measure_EE(psi, div(N, 2))
        SvN_center = measure_EE(psi, N)
        push!(Entropy_array, SvN)
        push!(Center_Entropy_array, SvN_center)

        zz_corr = real(correlation_matrix(psi, "Sz", "Sz")[N+1,:])
        sp_corr = real(correlation_matrix(psi, "S+", "S+")[N+1,:])
        sm_corr = real(correlation_matrix(psi, "S-", "S-")[N+1,:])

        push!(zz_corr_array, vec(zz_corr))
        push!(sp_corr_array, vec(sp_corr))
        push!(sm_corr_array, vec(sm_corr))

        psi = apply_dynamic_coupling(psi, s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, Ω, tau, dynamic_period, center_spin_period)
        psi = apply(gates, psi; cutoff=1e-9)

        # psi = apply_dynamic_coupling(psi, s, specified_sites, t,
        #                             gamma, Jx, Jy, Jz, Cx, Cy, Cz,
        #                             Ω, tau/2, dynamic_period, center_spin_period)
        # # 2) 整步静态门 H2（τ）
        # psi = apply(gates, psi; cutoff=1e-9)
        # 3) 半步动态耦合 H1（τ/2），时刻 t+τ/2
        # psi = apply_dynamic_coupling(psi, s, specified_sites, t + tau/2,
        #                              gamma, Jx, Jy, Jz, Cx, Cy, Cz,
        #                              Ω, tau/2, dynamic_period, center_spin_period)

        psi_z = apply_dynamic_coupling(psi, s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, Ω, tau, dynamic_period, center_spin_period)
        psi_z = apply(gates, psi; cutoff=1e-9)

        normalize!(psi)
        normalize!(psi_z)

        tcorr = inner(psi', H, psi_z)
        push!(tcorr_array, tcorr)
    end

    corrz_mean = hcat(CSz, Sz_array) .* Sz_array
    zz_matrix = Array(hcat(zz_corr_array...)') - corrz_mean / 4
    sp_matrix = Array(hcat(sp_corr_array...)')
    sm_matrix = Array(hcat(sm_corr_array...)')
    return t_array, Sz_array, Sy_array, Sx_array, Entropy_array, Center_Entropy_array, CSz_array, CSz, zz_matrix, sp_matrix, sm_matrix, tcorr_array
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
function create_plots(t_array, Sz_array, Sy_array, Sx_array, CSz_array, CSz, Entropy_array, Center_Entropy_array)
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
    t_array = t_array .+ 0.01
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
        title="Imbalance",
        xscale=:log10
    )

    # Create entropy plot
    entropy_plot = plot(
        t_array, Entropy_array, 
        label=L"$S_{bath}=-\sum_np_n\log p_n$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Entropy"
    )
    center_entropy_plot = plot(
        t_array, Center_Entropy_array,
        label=L"$S_{central}=-\sum_np_n\log p_n$",
        legend=:best,
        color=:purple,
        line=(:solid, 1.5),
        xlabel=L"$t$",
        ylabel="",
        title="Center Entropy"
    )

    # # Create zz Spatial correlations plot
    # zz_plot = plot(
    #     t_array, zz_matrix[:, 1], 
    #     label="1", 
    #     legend=:best, 
    #     color=:blue, 
    #     line=(:solid, 1.5), 
    #     xlabel=L"$t$", 
    #     ylabel=L"$\langle S^z_i S^z_j \rangle$", 
    #     title="zz Spatial Correlations"
    # )
    # for j in 2:size(zz_matrix, 2) - 1
    #     plot!(zz_plot, t_array, zz_matrix[:, j], label="$j", color=:auto, line=(:solid, 1.5))
    # end

    # zzt_plot = plot(
    #     t_array, real(tcorr_array), 
    #     # label="1", 
    #     legend=:best, 
    #     color=:blue, 
    #     line=(:solid, 1.5), 
    #     xlabel=L"$t$", 
    #     ylabel=L"$\langle S^z(t) S^z(0) \rangle$", 
    #     title="zz Time Correlations"
    # )

    # Combine and display all plots
    combined_plot = plot(spin_plot, imbalance_plot, center_entropy_plot, entropy_plot, layout=(2, 2), size=(1200, 1200))
    display(combined_plot)
    # display(imbalance_plot)
end

"""
Save simulation results

Parameters:
- folder: Path to the folder where data will be saved
- t_array: Array of time values
- Sz_array: Array of Sz expectation values
- Sy_array: Array of Sy expectation values
- Sx_array: Array of Sx expectation values
- Entropy_array: Array of entanglement entropy values
- CSz_array: Array of central spin Sz values
- CSz: Central spin Sz values
- parameters: Dictionary containing the simulation parameters
"""
function save_simulation_results(folder, t_array, Sz_array, Sy_array, Sx_array, Entropy_array, Center_Entropy_array, CSz_array, CSz, zz_matrix, sp_matrix, sm_matrix, tcorr_array, parameters)
    # Create folder if it doesn't exist
    mkpath(folder)

    # Save data to HDF5
    h5file = joinpath(folder, "simulation_results.h5")
    h5open(h5file, "w") do file
        write(file, "t_array", t_array)
        write(file, "Sz_array", Sz_array)
        write(file, "Sy_array", Sy_array)
        write(file, "Sx_array", Sx_array)
        write(file, "Entropy_array", Entropy_array)
        write(file, "Center_Entropy_array", Center_Entropy_array)
        write(file, "CSz_array", CSz_array)
        write(file, "CSz", CSz)
        write(file, "zz_matrix", zz_matrix)
        write(file, "sp_matrix", sp_matrix)
        write(file, "sm_matrix", sm_matrix)
        write(file, "tcorr_array", tcorr_array)
    end

    # Save parameters
    open(joinpath(folder, "parameters.txt"), "w") do io
        for (key, value) in parameters
            println(io, "$key = $value")
        end
    end

    # Create and save plots
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
    savefig(joinpath(folder, "spin_plot.png"))

    imbalance = sum(CSz_array, dims=2) ./ (sum(CSz ./ 2 .+ 0.5, dims=2))
    imbalance_plot = plot(
        t_array, imbalance, 
        label=L"$I=(S_{z,\uparrow}^{e}-S_{z,\uparrow}^{o})/(S_{z,\uparrow}^{o}+S_{z,\uparrow}^{e})$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Imbalance",
        xscale=:log10
    )
    savefig(joinpath(folder, "imbalance_plot.png"))

    entropy_plot = plot(
        t_array, Entropy_array, 
        label=L"$S_{bath}=-\sum_np_n\log p_n$", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel="", 
        title="Entropy"
    )
    savefig(joinpath(folder, "entropy_plot.png"))

    center_entropy_plot = plot(
        t_array, Center_Entropy_array,
        label=L"$S_{central}=-\sum_np_n\log p_n$",
        legend=:best,
        color=:purple,
        line=(:solid, 1.5),
        xlabel=L"$t$",
        ylabel="",
        title="Center Entropy"
    )
    savefig(joinpath(folder, "center_entropy_plot.png"))

    # Create zz Spatial correlations plot
    zz_plot = plot(
        t_array, zz_matrix[:, 1], 
        label="1", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel=L"$\langle S^z_i S^z_j \rangle$", 
        title="zz Spatial Correlations"
    )
    for j in 2:size(zz_matrix, 2) - 1
        plot!(zz_plot, t_array, zz_matrix[:, j], label="$j", color=:auto, line=(:solid, 1.5))
    end
    savefig(joinpath(folder, "zz_spatial_correlations.png"))

    zzt_plot = plot(
        t_array, real(tcorr_array), 
        # label="1", 
        legend=:best, 
        color=:blue, 
        line=(:solid, 1.5), 
        xlabel=L"$t$", 
        ylabel=L"$\langle S^z(t) S^z(0) \rangle$", 
        title="zz Time Correlations"
    )
    savefig(joinpath(folder, "zz_time_correlations.png"))

    combined_plot = plot(spin_plot, imbalance_plot, center_entropy_plot, entropy_plot, zz_plot, zzt_plot, layout=(3, 2), size=(1200, 1200))
    savefig(joinpath(folder, "combined_plot.png"))
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
- periodic: Whether to use periodic boundary conditions, true for periodic, false for open boundary
- α: Coefficient for the non-uniform field

This function initializes the spin chain, creates gate operations, performs time evolution, and plots the results.
"""
function run_simulation(N, Jz, Jy, Jx, Cz, Cy, Cx, Ω, F, W, gamma, tau, ttotal, center_spin_initial_state, chain_initial_state, specified_sites, dynamic_period, center_spin_period, periodic, α)
    # Initialize spin sites
    s = siteinds("S=1/2", N + 1)
    initial_state_function = create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    psi = MPS(s, initial_state_function)

    # Create gate operations
    gates = create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, α, specified_sites, periodic)

    # Perform time evolution
    t_array, Sz_array, Sy_array, Sx_array, Entropy_array, Center_Entropy_array, CSz_array, CSz, zz_matrix, sp_matrix, sm_matrix, tcorr_array = time_evolution(N, ttotal, tau, psi, gates, s, specified_sites, Jx, Jy, Jz, Cx, Cy, Cz, Ω, gamma, dynamic_period, center_spin_period)

    # Plot results
    create_plots(t_array, Sz_array, Sy_array, Sx_array, CSz_array, CSz, Entropy_array, Center_Entropy_array)

    # Save results
    parameters = OrderedDict(
        "N" => N,
        "Jz" => Jz,
        "Jy" => Jy,
        "Jx" => Jx,
        "Cz" => Cz,
        "Cy" => Cy,
        "Cx" => Cx,
        "Ω" => Ω,
        "F" => F,
        "W" => W,
        "gamma" => gamma,
        "α" => α,
        "tau" => tau,
        "ttotal" => ttotal,
        "center_spin_initial_state" => center_spin_initial_state,
        "chain_initial_state" => chain_initial_state,
        "specified_sites" => specified_sites,
        "dynamic_period" => dynamic_period,
        "center_spin_period" => center_spin_period,
        "periodic" => periodic
    )

    # Create timestamp folder
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    folder = joinpath("data", timestamp)
    mkpath(folder)

    # Save simulation results
    save_simulation_results(folder, t_array, Sz_array, Sy_array, Sx_array, Entropy_array, Center_Entropy_array, CSz_array, CSz, zz_matrix, sp_matrix, sm_matrix, tcorr_array, parameters)
end
