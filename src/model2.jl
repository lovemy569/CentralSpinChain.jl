using ITensors
using Plots
using HDF5
using Dates

#————————————————— 初态函数 —————————————————#
function create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    return n -> n == N + 1 ? center_spin_initial_state : (
        chain_initial_state == "Néel_state" ? (isodd(n) ? "Up" : "Dn") :
        chain_initial_state == "all_up"    ? "Up" :
        chain_initial_state == "all_dn"    ? "Dn" :
        chain_initial_state == "all_x+"    ? "X+" :
        chain_initial_state == "all_x-"    ? "X-" :
        error("Invalid chain initial state")
    )
end

#————————————————— 动态耦合函数 —————————————————#
function coupling_strength(t, Jx, Jy, Jz, period)
    if period == 0 return Jx, Jy, Jz end
    half = period/2
    on = mod(t,period) < half
    return on ? (0,0,0) : (Jx,Jy,Jz)
end

function center_spin_coupling_strength(t, Cx, Cy, Cz, period)
    if period == 0 return Cx, Cy, Cz end
    half = period/2
    on = mod(t,period) < half
    return on ? (0,0,0) : (Cx,Cy,Cz)
end

function apply_dynamic_coupling(psi, s, specified_sites, t,
                                gamma, Jx, Jy, Jz, Cx, Cy, Cz,
                                Ω, τ, period, center_period)
    # 1) 计算系数
    Jx_d, Jy_d, Jz_d = coupling_strength(t, Jx,Jy,Jz, period)
    Cx_d, Cy_d, Cz_d = center_spin_coupling_strength(t, Cx,Cy,Cz, center_period)
    center_site = s[end]

    # 2) 链-链动态耦合
    for (i,j) in specified_sites
        if i<length(s) && j<length(s)
            h = Jx_d*op("Sx",s[i])*op("Sx",s[j]) +
                Jy_d*op("Sy",s[i])*op("Sy",s[j]) +
                Jz_d*op("Sz",s[i])*op("Sz",s[j])
            G = exp(-1im*τ*h)
            psi = apply(G, psi; cutoff=1e-9)
        end
    end

    # 3) 链-中心动态耦合
    γ = gamma/(length(s)-1)
    center_pairs = Set((i==length(s) ? j : i) for (i,j) in specified_sites if i==length(s)||j==length(s))
    for j in center_pairs
        h = γ*( Cx_d*op("Sx",s[j])*op("Sx",center_site) +
               Cy_d*op("Sy",s[j])*op("Sy",center_site) +
               Cz_d*op("Sz",s[j])*op("Sz",center_site) )
        G = exp(-1im*τ*h)
        psi = apply(G, psi; cutoff=1e-9)
    end

    # 4) 中心自旋 Ω σᶻ
    Gz = exp(-1im*τ*(Ω*op("Sz",center_site)))
    psi = apply(Gz, psi; cutoff=1e-9)

    return psi
end

#————————————————— 纠缠熵测量 —————————————————#
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

#————————————————— 二阶 TEBD + 动态耦合 TDVP 演化 —————————————————#
function time_evolution(N, ttotal, τ, psi, s, specified_sites,
                        Jx, Jy, Jz, Cx, Cy, Cz, Ω, gamma,
                        F, W, α,
                        dynamic_period, center_spin_period,
                        periodic::Bool=false)

    # 随机场
    h = rand([-W, W], N)

    # 静态单体场 半步门 U_on(τ/2)
    field_half = [ exp(-1im*(τ/2)*( (-F*j + α*j^2/N^2 + h[j]) * op("Sz",s[j]) ))
                   for j in 1:N ]

    # 静态链内耦合：奇键与偶键 整步门 U_bond(τ)
    odd_bonds  = [ j for j in 1:2:(periodic ? N : N-1)
                    if !((j,j+1) in specified_sites || (j+1,j) in specified_sites) ]
    even_bonds = [ j for j in 2:2:((periodic && N>2) ? N : N-1)
                    if !((j,j+1) in specified_sites || (j+1,j) in specified_sites) ]

    odd_gates  = [ exp(-1im*τ*( Jx*op("Sx",s[j])*op("Sx",s[j+1]) +
                              Jy*op("Sy",s[j])*op("Sy",s[j+1]) +
                              Jz*op("Sz",s[j])*op("Sz",s[j+1]) )) for j in odd_bonds ]
    even_gates = [ exp(-1im*τ*( Jx*op("Sx",s[j])*op("Sx",s[j+1]) +
                               Jy*op("Sy",s[j])*op("Sy",s[j+1]) +
                               Jz*op("Sz",s[j])*op("Sz",s[j+1]) )) for j in even_bonds ]

    # 预分配结果
    num_steps = Int(round(ttotal/τ)) + 1
    t_array    = zeros(Float64, num_steps)
    Sz_array   = similar(t_array)
    Sy_array   = similar(t_array)
    Sx_array   = similar(t_array)
    Entropy    = similar(t_array)
    CEntropy   = similar(t_array)
    CSz_arr    = zeros(Float64, num_steps, N)
    Imbalance  = similar(t_array)
    zz_corrs   = Vector{Vector{Float64}}(undef, num_steps)
    sp_corrs   = Vector{Vector{Float64}}(undef, num_steps)
    sm_corrs   = Vector{Vector{Float64}}(undef, num_steps)
    tcorrs     = ComplexF64[]

    # 为关联函数准备 psi_z
    psi_z = deepcopy(psi)
    # 在中央自旋上施加 Sz
    A = op("Sz", s[end]) * psi_z[end]
    noprime!(A)
    psi_z[end] = A
    # 构造测量用 H = Sz_{central}
    os = OpSum(); os += "Sz", N+1
    Hcorr = MPO(os, s)

    step = 1
    @showprogress "Time Evolution" for t in 0:τ:ttotal
        # U1 半步 动态耦合
        psi   = apply_dynamic_coupling(psi,   s, specified_sites, t,
                                       gamma, Jx,Jy,Jz, Cx,Cy,Cz,
                                       Ω, τ/2, dynamic_period, center_spin_period)
        psi_z = apply_dynamic_coupling(psi_z, s, specified_sites, t,
                                       gamma, Jx,Jy,Jz, Cx,Cy,Cz,
                                       Ω, τ/2, dynamic_period, center_spin_period)

        # U2 二阶 TEBD 静态 H_2
        # 2.1 单体场 τ/2
        for g in field_half
            psi   = apply(g,   psi;   cutoff=1e-9)
            psi_z = apply(g,   psi_z; cutoff=1e-9)
        end
        # 2.2 奇键 τ
        for g in odd_gates
            psi   = apply(g,   psi;   cutoff=1e-9)
            psi_z = apply(g,   psi_z; cutoff=1e-9)
        end
        # 2.3 偶键 τ
        for g in even_gates
            psi   = apply(g,   psi;   cutoff=1e-9)
            psi_z = apply(g,   psi_z; cutoff=1e-9)
        end
        # 2.4 单体场 τ/2
        for g in field_half
            psi   = apply(g,   psi;   cutoff=1e-9)
            psi_z = apply(g,   psi_z; cutoff=1e-9)
        end

        # U1 半步 动态耦合 (t+τ/2)
        psi   = apply_dynamic_coupling(psi,   s, specified_sites, t+τ/2,
                                       gamma, Jx,Jy,Jz, Cx,Cy,Cz,
                                       Ω, τ/2, dynamic_period, center_spin_period)
        psi_z = apply_dynamic_coupling(psi_z, s, specified_sites, t+τ/2,
                                       gamma, Jx,Jy,Jz, Cx,Cy,Cz,
                                       Ω, τ/2, dynamic_period, center_spin_period)

        # — 测量 & 保存 —
        t_array[step]  = t
        Sz_array[step] = 2 * expect(complex(psi), "Sz"; sites=N+1)
        Sy_array[step] = 2 * expect(complex(psi), "Sy"; sites=N+1)
        Sx_array[step] = 2 * expect(complex(psi), "Sx"; sites=N+1)

        for j in 1:N
            szj = expect(complex(psi), "Sz"; sites=j)
            CSz_arr[step,j] = 2*szj
        end
        Imbalance[step] = sum(((-1)^(j-1)) * CSz_arr[step,j] for j in 1:N) / N

        Entropy[step]  = measure_EE(psi, div(N,2))
        CEntropy[step] = measure_EE(psi, N)

        zz_corrs[step] = real(correlation_matrix(psi, "Sz", "Sz")[N+1,:])
        sp_corrs[step] = real(correlation_matrix(psi, "S+", "S+")[N+1,:])
        sm_corrs[step] = real(correlation_matrix(psi, "S-", "S-")[N+1,:])

        # 计算 tcorr
        normalize!(psi)
        normalize!(psi_z)
        push!(tcorrs, inner(psi', Hcorr, psi_z))

        step += 1
    end

    # 将矢量化的 corr 转为矩阵
    zz_mat = hcat(zz_corrs...)'
    sp_mat = hcat(sp_corrs...)'
    sm_mat = hcat(sm_corrs...)'

    return t_array, Sz_array, Sy_array, Sx_array,
           Imbalance, CSz_arr,
           Entropy, CEntropy,
           zz_mat, sp_mat, sm_mat,
           tcorrs
end

#————————————————— 绘图 & 保存函数（保持不变） —————————————————#
function create_plots(t_array, Sz_array, Sy_array, Sx_array,
                      CSz_arr, Imbalance, Entropy, CEntropy)
    spin_plot = plot(t_array, Sz_array, label="⟨Sz⟩", xlabel="t", title="Central Spin")
    plot!(spin_plot, t_array, Sy_array, label="⟨Sy⟩")
    plot!(spin_plot, t_array, Sx_array, label="⟨Sx⟩")

    imb_plot = plot(t_array, Imbalance, xlabel="t", title="Imbalance")
    ent_plot = plot(t_array, Entropy, xlabel="t", title="Bath Entropy")
    cent_ent_plot = plot(t_array, CEntropy, xlabel="t", title="Central Entropy")

    display(plot(spin_plot, imb_plot, cent_ent_plot, ent_plot, layout=(2,2)))
end

function save_simulation_results(folder, t_array, Sz, Sy, Sx,
                                 Imbalance, CSz_arr, Entropy, CEntropy,
                                 zz_mat, sp_mat, sm_mat, tcorrs, params)
    mkpath(folder)
    h5open(joinpath(folder,"results.h5"),"w") do f
        write(f,"t",t_array)
        write(f,"Sz",Sz); write(f,"Sy",Sy); write(f,"Sx",Sx)
        write(f,"Imbalance",Imbalance)
        write(f,"CSz_arr",CSz_arr)
        write(f,"Entropy",Entropy); write(f,"CEntropy",CEntropy)
        write(f,"zz",zz_mat); write(f,"sp",sp_mat); write(f,"sm",sm_mat)
        write(f,"tcorr",tcorrs)
    end
    open(joinpath(folder,"params.txt"),"w") do io
        for (k,v) in params println(io,"$k = $v") end
    end
end

#————————————————— 运行函数 —————————————————#
function run_simulation(N, Jz, Jy, Jx,
                        Cz, Cy, Cx, Ω,
                        F, W, gamma,
                        τ, ttotal,
                        center_spin_initial_state,
                        chain_initial_state,
                        specified_sites,
                        dynamic_period,
                        center_spin_period,
                        periodic::Bool,
                        α)

    # 1) 初始化
    s = siteinds("S=1/2", N+1)
    init_fn = create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    psi = MPS(s, init_fn)

    # 2) 演化
    t_array, Sz, Sy, Sx,
    Imbalance, CSz_arr,
    Entropy, CEntropy,
    zz_mat, sp_mat, sm_mat,
    tcorrs = time_evolution(N, ttotal, τ, psi, s, specified_sites,
                             Jx, Jy, Jz, Cx, Cy, Cz, Ω, gamma,
                             F, W, α,
                             dynamic_period, center_spin_period,
                             periodic)

    # 3) 绘图 & 保存
    create_plots(t_array, Sz, Sy, Sx, CSz_arr, Imbalance, Entropy, CEntropy)

    params = OrderedDict(
        "N"=>N, "Jx"=>Jx, "Jy"=>Jy, "Jz"=>Jz,
        "Cx"=>Cx, "Cy"=>Cy, "Cz"=>Cz, "Ω"=>Ω,
        "F"=>F, "W"=>W, "γ"=>gamma,
        "τ"=>τ, "ttotal"=>ttotal,
        "center_state"=>center_spin_initial_state,
        "chain_state"=>chain_initial_state,
        "specified_sites"=>specified_sites,
        "dynamic_period"=>dynamic_period,
        "center_spin_period"=>center_spin_period,
        "periodic"=>periodic, "α"=>α
    )
    folder = joinpath("data", Dates.format(now(),"yyyy-mm-dd_HH-MM-SS"))
    save_simulation_results(folder, t_array, Sz, Sy, Sx,
                            Imbalance, CSz_arr, Entropy, CEntropy,
                            zz_mat, sp_mat, sm_mat, tcorrs, params)
end
