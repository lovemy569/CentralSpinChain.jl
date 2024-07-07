using CentralSpinChain

# 参数设置
N = 4  # 粒子数
Jz = 1.30  # 自旋链交换系数
Jy = 1.03  # 自旋链交换系数
Jx = 1.20  # 自旋链交换系数
Cz = 1.60  # 耦合系数
Cy = 1.08  # 耦合系数
Cx = 1.70  # 耦合系数
F = 0.88  # 线性场系数
W = 1.00  # 随机系数
gamma = 4.00  # 中心自旋耦合系数
α = 1.00  # 非均匀场的非线性系数

tau = 0.05  # 时间步长
ttotal = 5.00  # 总时间

# 初始链状态和中心自旋状态
center_spin_initial_state = "X+"  # 可选择 "Up", "Dn", "X+", "X-"
chain_initial_state = "all_up"  # 可选择 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-"

# 动态耦合参数
specified_sites = [(2, 3)]  # 指定site对,example: specified_sites = [(2, 1),(2, 3),(4, 3)]
dynamic_period = 1.0  # 动态耦合周期(设置为0取消floquet)
center_spin_period = 1.0  # 中心自旋的动态耦合周期(设置为0取消floquet)

# 周期性边界条件设置
periodic = true  # true 表示使用周期性边界条件，false 表示使用开放边界

# 运行函数得到结果，并保存到本目录文件下
run_simulation(N, Jz, Jy, Jx, Cz, Cy, Cx, F, W, gamma, tau, ttotal, center_spin_initial_state, chain_initial_state, specified_sites, dynamic_period, center_spin_period, periodic, α)


