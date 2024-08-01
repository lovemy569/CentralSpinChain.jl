using CentralSpinChain

# 参数设置
N = 12  # 粒子数
Jz = 1.00  # 自旋链交换系数
Jy = 1.00  # 自旋链交换系数
Jx = 1.00  # 自旋链交换系数
Cz = 1.00  # 耦合系数
Cy = 1.00  # 耦合系数
Cx = 1.00  # 耦合系数
F = -0.0  # 线性场系数
W = 0.00  # 随机系数
gamma = 4.00  # 中心自旋耦合系数
α = 0.0  # 非均匀场的非线性系数

tau = 0.05  # 时间步长
ttotal = 50.00  # 总时间

# 初始链状态和中心自旋状态
center_spin_initial_state = "X+"  # 可选择 "Up", "Dn", "X+", "X-"
chain_initial_state = "Néel_state"  # 可选择 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-"

# 动态耦合参数
specified_sites = []  # 指定site对,example: specified_sites = [(2, 1),(2, 3),(4, 3)]
dynamic_period = 0  # 动态耦合周期(设置为0取消floquet)
center_spin_period = 0  # 中心自旋的动态耦合周期(设置为0取消floquet)

# 周期性边界条件设置
periodic = false  # true 表示使用周期性边界条件，false 表示使用开放边界

# 运行函数得到结果，并保存到data文件中
run_simulation(N, Jz, Jy, Jx, Cz, Cy, Cx, F, W, gamma, tau, ttotal, center_spin_initial_state, chain_initial_state, specified_sites, dynamic_period, center_spin_period, periodic, α)


