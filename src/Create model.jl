using ITensors

"""
示例代码：
- center_spin_initial_state = "Up"  # 可选择 "Up", "Dn", "X+", "X-"
- chain_initial_state = "Néel_state"  # 可选择 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-"
- chain_initial_state = ["Up", "Dn", "X+", "X-", "Up", "Dn", "X+", "X-", "Up", "Dn"]  # 或者自定义
"""
function create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    return n -> n == N + 1 ? center_spin_initial_state : (  # 设置中心自旋状态
        isa(chain_initial_state, AbstractVector) && length(chain_initial_state) == N ? chain_initial_state[n] :
        chain_initial_state == "Néel_state" ? (isodd(n) ? "Up" : "Dn") :
        chain_initial_state == "all_up" ? "Up" :
        chain_initial_state == "all_dn" ? "Dn" :  # 设置自旋链状态
        chain_initial_state == "all_x+" ? "X+" :
        chain_initial_state == "all_x-" ? "X-" : error("Invalid chain initial state or array length")
    )
end

function create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, α, specified_sites, periodic::Bool)
    center_spin_site = N + 1  # 指定第N+1个site为中心自旋
    count = 0  # 计数关于中心自旋的需要特殊处理的site对，用于gate初始化
    # 指定site对,例如specified_sites = [(2, 1), (2, 3), (4, 3)]
    # 从specified_sites中剔除掉中心自旋，本函数只处理自旋链
    for (i, j) in specified_sites
        if i == center_spin_site || j == center_spin_site
            count += 1  #如果指定的site存在中心自旋，那么计数
        end
    end

    num_gates = (periodic ? N : N - 1) - length(specified_sites) + N + count  # 计算总门数，同时处理是否周期
    gates = Vector{ITensor}(undef, num_gates)  # 门序列初始化
    next_site = periodic ? [mod(j, N) + 1 for j in 1:N] : [j + 1 for j in 1:N-1]  # 生成下一个site序号，用于site对指定
    gate_idx = 1 

    # Create gate operations for spin chain coupling
    for j in 1:(periodic ? N : N-1)  # 对于开放边界条件跳过首尾相互作用
        s1 = s[j]
        s2 = s[next_site[j]]
        if (j, next_site[j]) in specified_sites || (next_site[j], j) in specified_sites
            continue  # 跳过指定的自旋链site对
        end
        # 定义相互作用方式
        hj = Jx * op("Sx", s1) * op("Sx", s2) + 
             Jy * op("Sy", s1) * op("Sy", s2) +
             Jz * op("Sz", s1) * op("Sz", s2)
        # 创建一阶Trotter分解
        gates[gate_idx] = exp(-im * tau * hj)
        gate_idx += 1
    end

    # Generate uniformly distributed random field
    h = rand([-W, W], N)  # MBL随机场
    L = N
    for j in 1:N
        Wj = -F * j + α * j^2 / (L-1)^2  # stark线性场
        fz = Wj * op("Sz", s[j]) + h[j] * op("Sz", s[j])
        # 创建一阶Trotter分解
        gates[gate_idx] = exp(-im * tau * fz)
        gate_idx += 1
    end

    return gates
end

function coupling_strength(t, Jx, Jy, Jz, period)
    if period == 0
        return Jx, Jy, Jz  #设置开关，取消周期，使用默认强度
    end 

    half_period = period / 2  
    in_first_half = mod(t, period) < half_period  # 如果在前半周期，强度为0；后半周期为指定
    
    Jx_dynamic = in_first_half ? 0.0 : Jx
    Jy_dynamic = in_first_half ? 0.0 : Jy
    Jz_dynamic = in_first_half ? 0.0 : Jz
    
    return Jx_dynamic, Jy_dynamic, Jz_dynamic
end

function center_spin_coupling_strength(t, Cx, Cy, Cz, center_spin_period)
    if center_spin_period == 0
        return Cx, Cy, Cz
        # return 0,0,0
    end

    half_period = center_spin_period / 2  
    in_first_half = mod(t, center_spin_period) < half_period  
    
    Cx_dynamic = in_first_half ? 0.0 : Cx
    Cy_dynamic = in_first_half ? 0.0 : Cy
    Cz_dynamic = in_first_half ? 0.0 : Cz
    
    return Cx_dynamic, Cy_dynamic, Cz_dynamic
end

function create_dynamic_coupling_gates(s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, Ω, tau, period, center_spin_period)
    # 计算动态耦合系数
    Jx_dynamic, Jy_dynamic, Jz_dynamic = coupling_strength(t, Jx, Jy, Jz, period)
    Cx_dynamic, Cy_dynamic, Cz_dynamic = center_spin_coupling_strength(t, Cx, Cy, Cz, center_spin_period)
    center_spin = s[end]  # 中心自旋的site索引

    gates = Vector{ITensor}(undef, length(specified_sites) + length(s))  # 初始化门的序列
    gate_idx = 1

    # 处理自旋链上的动态耦合
    for (j1, j2) in specified_sites
        if j1 != length(s) && j2 != length(s)
            s1, s2 = s[j1], s[j2]
            hj = Jx_dynamic * op("Sx", s1) * op("Sx", s2) + 
                 Jy_dynamic * op("Sy", s1) * op("Sy", s2) + 
                 Jz_dynamic * op("Sz", s1) * op("Sz", s2)
            G_dynamic = exp(-im * tau * hj)
            gates[gate_idx] = G_dynamic
            gate_idx += 1
        end
    end

    # 处理中心自旋的动态耦合
    gamma_over_len_s = gamma / (length(s) - 1)
    # 挑选出所有与中心自旋相关的自旋链site
    specified_center_spin_sites = Set((j1 == length(s) ? j2 : j1) for (j1, j2) in specified_sites if j1 == length(s) || j2 == length(s))
    
    for j in specified_center_spin_sites
        sj = s[j]
        hjc_dynamic = gamma_over_len_s * (Cx_dynamic * op("Sx", sj) * op("Sx", center_spin) +
                                          Cy_dynamic * op("Sy", sj) * op("Sy", center_spin) +
                                          Cz_dynamic * op("Sz", sj) * op("Sz", center_spin))
        Gjc_dynamic = exp(-im * tau * hjc_dynamic)
        gates[gate_idx] = Gjc_dynamic
        gate_idx += 1
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
        gates[gate_idx] = Gjc_default
        gate_idx += 1
    end

    # 中心自旋的额外项 \Omega * σz
    omega_z_term = 2 * Ω * op("Sz", center_spin)
    G_omega_z = exp(-im * tau * omega_z_term)
    gates[gate_idx] = G_omega_z

    return gates
end

function create_model(N, s, Jx, Jy, Jz, tau, F, W, α, specified_sites, periodic, t, gamma, Cx, Cy, Cz, Ω, period, center_spin_period, center_spin_initial_state, chain_initial_state)
    # 创建静态的自旋链 gates
    static_gates = create_spin_chain_gates(N, s, Jx, Jy, Jz, tau, F, W, α, specified_sites, periodic)
    
    # 创建动态耦合 gates
    dynamic_gates = create_dynamic_coupling_gates(s, specified_sites, t, gamma, Jx, Jy, Jz, Cx, Cy, Cz, Ω, tau, period, center_spin_period)
    
    # 合并 gates
    gates = vcat(static_gates, dynamic_gates)

    # 初始化 psi
    initial_state_function = create_initial_state_function(N, center_spin_initial_state, chain_initial_state)
    psi = MPS(s, initial_state_function)

    return psi, gates
end

