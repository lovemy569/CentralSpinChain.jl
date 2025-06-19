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