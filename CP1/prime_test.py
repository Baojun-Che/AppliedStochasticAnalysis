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
    
    N = 100  
    h = 0
    q = 10
    ## q=10, T以粗尺度变化
    # temperatures = np.linspace(0.1, 5.0, 50)
    # internal_energies, specific_heats = [], []
    # _, _, lattice = mcmc_without_external_field(N, q, 0.1, n_tempering=500, n_measure=1000, get_energy=True)
    # for T in temperatures:
    #     results, _, lattice = mcmc_without_external_field(N, q, T, lattice=lattice, n_tempering=0, n_measure=2000, get_energy=True)
    #     internal_energies.append(results['internal_energy'])
    #     specific_heats.append(results['specific_heat'])
    #     print(f"u={results['internal_energy']}, c={results['specific_heat']}" )
    
    # fig_name = "uc-T,q=10,large"
    # plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name)


    ### q=10, T以中等尺度变化
    # temperatures = np.linspace(0.60, 0.90, 16)
    # internal_energies, specific_heats = [], []
    # _, _, lattice = mcmc_without_external_field(N, q, 0.60, n_tempering=0, n_measure=4000, get_energy=True)
    # for T in temperatures:
    #     results, _, lattice = mcmc_without_external_field(N, q, T, lattice=lattice, n_tempering=0, n_measure=2000, get_energy=True)
    #     internal_energies.append(results['internal_energy'])
    #     specific_heats.append(results['specific_heat'])
    #     print(f"u={results['internal_energy']}, c={results['specific_heat']}" )
    
    # fig_name = "uc-T,q=10,middle"
    # plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name)

    ## q=10, T以细尺度变化
    # temperatures = np.linspace(0.660, 0.720, 13)
    # internal_energies, specific_heats = [], []
    # _, _, lattice = mcmc_without_external_field(N, q, 0.680, n_tempering=0, n_measure=8000, get_energy=True)
    # for T in temperatures:
    #     results, _, lattice = mcmc_without_external_field(N, q, T, lattice=lattice, n_tempering=0, n_measure=2000, get_energy=True)
    #     internal_energies.append(results['internal_energy'])
    #     specific_heats.append(results['specific_heat'])
    #     print(f"u={results['internal_energy']}, c={results['specific_heat']}" )
    
    # fig_name = "uc-T,q=10,small"
    # plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name)
    
    # # 保存为txt文件
    # data = np.column_stack((temperatures, internal_energies, specific_heats))
    # np.savetxt('results/q=10,N=100.txt', data, header='Temperature Internal_Energy Specific_Heat', 
    #        fmt='%.6f', delimiter=' ')
    

    ##  q=10, T以细尺度变化, 加密到N=500
    temperatures = np.linspace(0.660, 0.720, 13)
    internal_energies, specific_heats = [], []
    _, _, lattice = mcmc_without_external_field(500, q, 0.660, n_tempering=0, n_measure=400000, RATE=1, get_energy=True, mes_energy=False)
    for T in temperatures:
        results, _, lattice = mcmc_without_external_field(500, q, T, lattice = lattice, n_tempering=0, n_measure=4000, RATE=5, get_energy=True)
        internal_energies.append(results['internal_energy'])
        specific_heats.append(results['specific_heat'])
        print(f"u={results['internal_energy']}, c={results['specific_heat']}" )
    
    fig_name = "uc-T,q=10,N=500"
    plot_temperature_enengy(temperatures, internal_energies, specific_heats, fig_name)
    
    # 保存为txt文件
    data = np.column_stack((temperatures, internal_energies, specific_heats))
    np.savetxt('results/q=10,N=500.txt', data, header='Temperature Internal_Energy Specific_Heat', 
           fmt='%.6f', delimiter=' ')
