using CentralSpinChain

# 参数设置
N = 4  # 粒子数
Jz = 1.00  # 自旋链交换系数
Jy = 1.00  # 自旋链交换系数
Jx = 1.00  # 自旋链交换系数
Cz = 1.00  # 耦合系数
Cy = 1.00  # 耦合系数
Cx = 1.00  # 耦合系数
F = 0.00   # 线性场系数
W = 0.00   # 随机系数
gamma = 5.00  # 中心自旋耦合系数

tau = 0.05  # 时间步长
ttotal = 5.00  # 总时间
cutoff = 1E-9  # 截断精度

# 初始链状态和中心自旋状态
center_spin_initial_state = "X+"
chain_initial_state = "all_up"  # 可选择 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-"

# 动态耦合参数
# specified_sites = [(2, 3), (5, 6)]  # 指定site之间的耦合对
specified_sites = [(2, 3)]
dynamic_period = 0.2  # 动态耦合周期
center_spin_period = 4  # 中心自旋的动态耦合周期

# 运行函数得到结果的图像
run_simulation(N, Jz, Jy, Jx, Cz, Cy, Cx, F, W, gamma, tau, ttotal, cutoff, center_spin_initial_state, chain_initial_state, specified_sites, dynamic_period, center_spin_period)


