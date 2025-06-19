using CentralSpinChain
# using ITensors

# 参数设置
N = 100  # 粒子数
Jz = 1.00  # 自旋链交换系数
Jy = 1.00  # 自旋链交换系数
Jx = 1.00  # 自旋链交换系数
Cz = 1.00  # 耦合系数
Cy = 1.00  # 耦合系数
Cx = 1.00  # 耦合系数
Ω = 10.00   # Zeeman
F = 5.0  # 线性场系数
W = 0.00  # 随机系数
gamma = 1.00  # 中心自旋耦合系数
α = 1.0  # 非均匀场的非线性系数

tau = 0.01  # 时间步长
ttotal = 10.00  # 总时间

# 初始链状态和中心自旋状态
center_spin_initial_state = "X+"  # 可选择 "Up", "Dn", "X+", "X-"
chain_initial_state = "Néel_state"  # 可选择 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-",或者自定义如下
# chain_initial_state = ["Up", "Dn", "X+", "X-", "Up", "Dn", "X+", "X-", "Up", "Dn"]

# 动态耦合参数，中心自旋是第N+1个，指定时注意
specified_sites = []  # 指定site对,例如specified_sites = [(2, 1), (2, 3), (4, 3)]
dynamic_period = 0  # 动态耦合周期(设置为0取消Floquet)
center_spin_period = 0  # 中心自旋的动态耦合周期(设置为0取消Floquet)

# 是否需要中心自旋 true or false
do_you_need = false
if do_you_need == false
    center_spin_period = ttotal * 3
    center_spin_sites = [(N + 1, i) for i in 1:N]
    specified_sites = append!(specified_sites, center_spin_sites)
end

# 周期性边界条件设置
periodic = true  # true 表示使用周期性边界条件，false 表示使用开放边界

# 运行函数得到结果，并保存到data文件中
run_simulation(
    N,                      # 自旋链的长度
    Jz,                     # Z方向的交换系数
    Jy,                     # Y方向的交换系数
    Jx,                     # X方向的交换系数
    Cz,                     # 中心自旋与链上自旋Z方向的耦合系数
    Cy,                     # 中心自旋与链上自旋Y方向的耦合系数
    Cx,                     # 中心自旋与链上自旋X方向的耦合系数
    Ω,                      # Zeeman
    F,                      # 线性场系数
    W,                      # 随机场的范围
    gamma,                  # 中心自旋耦合系数
    tau,                    # 时间步长
    ttotal,                 # 总时间
    center_spin_initial_state,  # 中心自旋的初始状态
    chain_initial_state,        # 自旋链的初始状态
    specified_sites,        # 指定的site对
    dynamic_period,         # 动态耦合周期
    center_spin_period,     # 中心自旋的动态耦合周期
    periodic,               # 是否使用周期性边界条件
    α                       # 非均匀场的非线性系数
)

# sites = siteinds("S=1/2", N + 1)
# psi = randomMPS(sites)
# psi_x = deepcopy(psi)

# #apply Sx to psi_x 
# opp = op("Sx",sites[1])
# A_1_x = opp * psi_x[1]
# noprime!(A_1_x)
# psi_x[1] = A_1_x

# os = OpSum()
# os += "Sx",1
# H = MPO(os,sites)
# val = inner(psi',H,psi_x)
# println("<Sx_1(0)Sx_1(0)> = $val")