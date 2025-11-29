from Potts_Model_2D import mcmc_without_external_field
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def plot_correlation_length(temperatures, corr_length, fig_name):
    
    if len(temperatures) != len(corr_length):
        raise ValueError("temperatures和corr_length的长度必须相同")
    
    sorted_indices = np.argsort(temperatures)
    
    sorted_temperatures = temperatures[sorted_indices]
    sorted_corr_length = corr_length[sorted_indices]
    
    plt.figure(figsize=(10, 6))
    
    plt.plot(sorted_temperatures, sorted_corr_length, 'o-', linewidth=2, markersize=6)
    plt.xlabel('Temperature', fontsize=14)
    plt.ylabel('Correlation Length', fontsize=14)
    plt.tight_layout()

    plt.savefig(fig_name + '.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    return sorted_temperatures, sorted_corr_length

# 运行分析
if __name__ == "__main__":

    q = 3
    N = 100
    
    corr_k = np.arange(10, 41, 5)
    temperatures = []
    corr_length = []
    r_square = []

    # _, _, lattice = mcmc_without_external_field(N, q, 0.995, n_tempering=200, n_measure=10000, RATE=1, mes_energy=False, get_energy=True)
    # np.save('lattice/N=100,q=3,T=0.995.npy', lattice)

    # lattice =  np.load('lattice/N=100,q=3,T=0.995.npy')
    # for T in [0.995, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4]:
    #     corr_gamma = np.zeros(len(corr_k))
    #     for i in range(5):
    #         results = mcmc_without_external_field(N, q, T, n_tempering=0, n_measure=1000, RATE=4, n_step=2, mes_energy=False, corr_k=corr_k)
    #         corr_gamma += np.maximum(results["corr_gamma"], 0)
    #     corr_gamma /= 5
    #     print(corr_gamma)
    #     reg = stats.linregress(corr_k, np.log(corr_gamma))
    #     print(f"回归R²值: {reg.rvalue**2}, p值: {reg.pvalue}. 相干长度 = {-1/reg.slope}")
    #     temperatures.append(T)
    #     corr_length.append(-1/reg.slope)
    #     r_square.append(reg.rvalue**2)
    A = []
    lattice_0 =  np.load('lattice/N=100,q=3,T=0.995.npy')
    for T in [0.98, 0.95, 0.90, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4]:
        corr_gamma = np.zeros(len(corr_k))
        for i in range(5):
            results,_,lattice = mcmc_without_external_field(N, q, T, n_tempering=0, n_measure=1000, RATE=2, n_step=1, mes_energy=False, corr_k=corr_k, get_energy=True, lattice=lattice_0)
            corr_gamma += np.maximum(results["corr_gamma"], 0)
        lattice_0 = lattice
        corr_gamma /= 5
        print(corr_gamma)
        reg = stats.linregress(corr_k, np.log(corr_gamma))
        print(f"回归R²值: {reg.rvalue**2}, p值: {reg.pvalue}. 相干长度 = {-1/reg.slope}")
        temperatures.append(T)
        corr_length.append(-1/reg.slope)
        r_square.append(reg.rvalue**2)
        A.append(corr_gamma)

    temperatures = np.array(temperatures)
    r_square = np.array(r_square)
    corr_length = np.array(corr_length)
    # indice =  np.where(r_square > 0.75)[0]

    data = np.column_stack((temperatures, corr_length, r_square))
    np.savetxt('results/q=3,corr_len.txt', data, header='Temperature Correlation_Length Reg_R_Square', fmt='%.6f', delimiter=' ')
    
    np.savetxt("corr_gamma_list", A)