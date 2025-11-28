from Potts_Model_2D import mcmc_without_external_field
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


def estimate_constant(T_values, y_values, T_star, arg_name):
    
    logT = np.log(np.abs(1 - T_values / T_star))
    logy = np.log(y_values)
    
    # 移除无穷大或NaN值
    valid_indices = np.isfinite(logT) & np.isfinite(logy)
    logT = logT[valid_indices]
    logy = logy[valid_indices]
    T_values = T_values[valid_indices]
    y_values = y_values[valid_indices]
    
    # 线性回归
    slope, intercept, r_value, p_value, _ = stats.linregress(logT, logy)
    r_squared = r_value ** 2
    
    # 绘制图形
    plt.figure(figsize=(10, 4))
    
    # 子图1: y 关于 T 的图像
    plt.subplot(1, 2, 1)
    plt.plot(T_values, y_values, 'b-', linewidth=2)
    plt.axvline(x=T_star, color='r', linestyle='--', alpha=0.7, label=f'T_star = {T_star}')
    plt.xlabel('T')
    plt.ylabel(arg_name)
    
    # 子图2: ln y 关于 ln |1-T/T_star| 的散点图和回归直线
    plt.subplot(1, 2, 2)
    plt.scatter(logT, logy, alpha=0.7)
    plt.xlabel('ln epsilon')
    plt.ylabel('ln '+arg_name)
    x_fit = np.linspace(min(logT), max(logT), 100)
    y_fit = intercept + slope * x_fit
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'y = {intercept:.3f} {slope:+.3f}x')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f"results/epsilon+{arg_name}.pdf")
    plt.show()
    
    print("回归结果:")
    print(f"gamma 估计值: {-slope:.6f}")
    print(f"回归R²值: {r_squared:.6f}, p值: {p_value}")
    
    return -slope

# 运行分析
if __name__ == "__main__":

    q = 3
    N = 30
    
    temperatures = [0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1]
    # temperatures = [0.9, 0.95, 1.0, 1.05, 1.1]

    corr_k = np.arange(5, 15, 1)
    corr_length = []
    r_square = []
    # _, _, lattice = mcmc_without_external_field(N, q, np.min(temperatures), n_tempering=200, n_measure=10000, RATE=1, mes_energy=False, get_energy=True)
    # np.save('lattice/N=100,T=0.50.npy', lattice)
    # lattice =  np.load('lattice/N=100,T=0.50.npy')

    for T in temperatures:
        results,_,lattice = mcmc_without_external_field(N, q, T, n_tempering=0 , n_measure=200, RATE=5, n_step=10, mes_energy=False, corr_k=corr_k, get_energy=True)
        corr_gamma = np.maximum(results["corr_gamma"], 1e-8)
        reg = stats.linregress(corr_k, np.log(corr_gamma))
        print(corr_gamma)
        print(f"回归R²值: {reg.rvalue**2}, p值: {reg.pvalue}. 相干长度 = {-1/reg.slope}")
        corr_length.append(-1/reg.slope)
        r_square.append(reg.rvalue**2)

    temperatures = np.array(temperatures)
    r_square = np.array(r_square)
    corr_length = np.array(corr_length)
    indice =  np.where(r_square > 0.75)[0]

    data = np.column_stack((temperatures[indice], corr_length[indice], r_square[indice]))
    np.savetxt('results/q=3,corr_len.txt', data, header='Temperature Correlation_Length Reg_R_Square', fmt='%.6f', delimiter=' ')
    
    