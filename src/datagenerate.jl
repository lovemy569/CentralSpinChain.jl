module DataGenerate

using CentralSpinChain.Model
using ITensors
using Plots
using HDF5
using ProgressMeter

export run_simulation, plot_results

function run_simulation()
    N = 12  # 粒子数
    Jx, Jy, Jz = 1.00, 1.00, 1.00  # 自旋链交换系数
    F, W = 0.00, 0.20  # 线性场系数和随机系数
    gamma = 5.00  # 中心自旋耦合系数
    Cx, Cy, Cz = 1.00, 1.00, 1.00  # 耦合系数
    tau = 0.05  # 时间步长
    ttotal = 500.00  # 总时间

    center_spin_initial_state = "X+"
    chain_initial_state = "all_up"  # 链的初态，可选 "Néel_state", "all_up", "all_dn", "all_x+", "all_x-"
    specified_sites = [(2, 3), (5, 6)]  # 指定site之间的耦合对

    psi, s = initialize_mps(N, center_spin_initial_state, chain_initial_state)
    gates = create_static_gates(N, Jx, Jy, Jz, F, W, gamma, Cx, Cy, Cz, s, specified_sites, tau)

    Sz_array, Sy_array, Sx_array, Entropy_array, t_array = Float64[], Float64[], Float64[], Float64[], Float64[]
    num_steps = Int(ttotal / tau) + 1
    CSz_array = zeros(Float64, num_steps, N)
    progress = Progress(num_steps, desc="Time Evolution", barlen=100)  # 创建进度条

    step_index = 1
    for t in 0:tau:ttotal
        Sz = expect(psi, "Sz"; sites = N + 1) * 2
        Sy = expect(psi, "Sy"; sites = N + 1) * 2
        Sx = expect(psi, "Sx"; sites = N + 1) * 2
        push!(t_array, t)
        push!(Sz_array, Sz)
        push!(Sy_array, Sy)
        push!(Sx_array, Sx)

        for j in 1:N
            CSz_array[step_index, j] = (-1) ^ (j + 1) * (expect(psi, "Sz"; sites = j) + 0.5)
        end
        SvN = measure_EE(psi, N)
        push!(Entropy_array, SvN)

        # 动态调整指定site之间的耦合强度
        psi = apply_dynamic_coupling(psi, s, specified_sites, t, Jx, Jy, Jz, tau)

        # 应用静态门操作进行时间演化
        psi = apply(gates, psi; cutoff=1e-9)
        normalize!(psi)

        next!(progress)  # 更新进度条
        step_index += 1
    end

    save_results(t_array, Sz_array, Sy_array, Sx_array, Entropy_array, CSz_array, CSz)
end

function save_results(t_array, Sz_array, Sy_array, Sx_array, Entropy_array, CSz_array, CSz)
    # 动态生成文件名
    filename = string("central", "_chain", N, "_F", F, "_W", W, "_J", gamma, ".h5")

    # 打开HDF5文件，准备写入
    f = h5open(filename, "w")
    write(f, "t", t_array)
    write(f, "Sz", Sz_array)
    write(f, "Sy", Sy_array)
    write(f, "Sx", Sx_array)
    write(f, "CSz", CSz)
    # 关闭HDF5文件
    close(f)
end

function plot_results()
    # 从HDF5文件加载数据
    filename = string("central", "_chain", N, "_F", F, "_W", W, "_J", gamma, ".h5")
    t_array = h5read(filename, "t")
    Sz_array = h5read(filename, "Sz")
    Sy_array = h5read(filename, "Sy")
    Sx_array = h5read(filename, "Sx")
    CSz = h5read(filename, "CSz")
    Entropy_array = h5read(filename, "Entropy_array")

    # 创建中心自旋的图表
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

    # 创建不平衡度的图表
    imbalance = sum(CSz_array, dims=2) ./ (sum(CSz./ 2 .+ 0.5, dims=2) )
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

    # 创建纠缠熵的图表
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

    # 使用 layout 创建子图布局
    combined_plot = plot(spin_plot, imbalance_plot, entropy_plot, layout=(1, 3), size=(1080, 320))

    # 保存图像
    savefig(combined_plot, "combined_plot.png")
end

end  # module DataGenerate
