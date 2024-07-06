module Model

using ITensors

export initialize_mps, create_static_gates, measure_EE, coupling_strength, apply_dynamic_coupling

# 初始化MPS
function initialize_mps(N, center_spin_initial_state, chain_initial_state="Néel_state")
    s = siteinds("S=1/2", N + 1)
    
    # 根据 chain_initial_state 参数设置链的初态
    initial_state_function = n -> n == N + 1 ? center_spin_initial_state : (
        chain_initial_state == "Néel_state" ? (isodd(n) ? "Up" : "Dn") :
        chain_initial_state == "all_up" ? "Up" :
        chain_initial_state == "all_dn" ? "Dn" :
        chain_initial_state == "all_x+" ? "X+" :
        chain_initial_state == "all_x-" ? "X-" : error("Invalid chain initial state")
    )
    
    psi = MPS(s, initial_state_function)
    return psi, s
end

# 创建静态门操作
function create_static_gates(N, Jx, Jy, Jz, F, W, gamma, Cx, Cy, Cz, s, specified_sites, tau)
    gates = ITensor[]
    for j in 1:N
        s1, s2 = s[j], s[mod(j, N) + 1]  # 周期性边界条件
        if any(x -> (x == (j, mod(j, N) + 1) || x == (mod(j, N) + 1, j)), specified_sites)
            continue  # 跳过指定site之间的耦合
        end
        hj = Jx * op("Sx", s1) * op("Sx", s2) + Jy * op("Sy", s1) * op("Sy", s2) + Jz * op("Sz", s1) * op("Sz", s2)
        push!(gates, exp(-im * tau * hj))
    end

    h = rand([-W, W], N)  # 均匀分布
    for j in 1:N
        fz = F * j * op("Sz", s[j]) + h[j] * op("Sz", s[j])
        push!(gates, exp(-im * tau * fz))
    end

    center_spin = s[N + 1]
    for j in 1:N
        hjc = gamma / N * (Cx * op("Sx", s[j]) * op("Sx", center_spin) + Cy * op("Sy", s[j]) * op("Sy", center_spin) + Cz * op("Sz", s[j]) * op("Sz", center_spin))
        push!(gates, exp(-im * tau * hjc))
    end

    return gates
end

# 纠缠熵计算
function measure_EE(psi, n)
    orthogonalize!(psi, n)
    _, S = svd(psi[n], (linkind(psi, n-1), siteinds(psi,n)); cutoff=1e-18, alg="qr_iteration")
    SvN = -sum(p * log(p) for p in diag(S).^2)
    return SvN
end

# 周期性变化的函数，决定指定site之间的耦合强度
function coupling_strength(t, Jz, period)
    if mod(t, 2 * period) < period
        return Jz
    else
        return 0.0
    end
end

# 应用动态耦合
function apply_dynamic_coupling(psi, s, specified_sites, t, Jx, Jy, Jz, tau)
    for (j1, j2) in specified_sites
        J_dynamic = coupling_strength(t, Jz, 50.0)  # 假设周期为50
        s1, s2 = s[j1], s[j2]
        hj = Jx * op("Sx", s1) * op("Sx", s2) + Jy * op("Sy", s1) * op("Sy", s2) + J_dynamic * op("Sz", s1) * op("Sz", s2)
        G_dynamic = exp(-im * tau * hj)
        psi = apply(G_dynamic, psi; cutoff=1e-9)
    end
    return psi
end

end  # module Model
