from Potts_Model_2D import mcmc_without_external_field
import matplotlib.pyplot as plt
import numpy as np
import os
def plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name):
    """绘制结果"""
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))  
    
    # internal energy vs temperature
    axes[0].plot(temperatures, internal_energies)
    axes[0].set_xlabel("Temperature")
    axes[0].set_ylabel("Internal Energy")
    
    # specific heat vs temperature  
    axes[1].plot(temperatures, specific_heats)
    axes[1].set_xlabel("Temperature")
    axes[1].set_ylabel("Specific Heat")
    
    os.makedirs('results', exist_ok=True)
    plt.tight_layout()
    plt.savefig(f"results/{fig_name}.pdf")
    plt.show()
    


if __name__ == "__main__":
    
    N = 500
    q = 3
    ###  q=3, T以细尺度变化, 加密到N=500
    temperatures = np.linspace(0.950, 1.005, 12)
    internal_energies, specific_heats = [], []
    lattice = np.ones((N, N), dtype=int)
    _, _, lattice = mcmc_without_external_field(500, q, 0.95, lattice = lattice, n_tempering=0, n_measure=200, RATE=1, get_energy=True, mes_energy=False)
    for T in temperatures:
        results, _, lattice = mcmc_without_external_field(500, q, T, lattice = lattice, n_tempering=0, n_measure=200, RATE=2, get_energy=True)
        internal_energies.append(results['internal_energy'])
        specific_heats.append(results['specific_heat'])
        print(f"u={results['internal_energy']}, c={results['specific_heat']}" )
    
    fig_name = "uc-T,q=3,N=500"
    plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name)
    
    # 保存为txt文件
    data = np.column_stack((temperatures, internal_energies, specific_heats))
    np.savetxt('results/q=3,N=500.txt', data, header='Temperature Internal_Energy Specific_Heat', 
           fmt='%.6f', delimiter=' ')
   
