from Potts_Model_2D import MCMC_simulation
import matplotlib.pyplot as plt
import numpy as np
def plot_temperature_enengy(temperatures, internal_energies, specific_heats):
    """绘制结果"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    

# 示例使用
if __name__ == "__main__":
    # 模拟参数
    N = 100  # 晶格尺寸
    temperatures = np.linspace(0.1, 5, 20)
    internal_energies, specific_heats = [], []
    for T in temperatures:
        q = 3
        results = MCMC_simulation(N, q, T, n_tempering=200, n_measure=1000, n_step=5)
        internal_energies.append(results['internal_energy'])
        specific_heats.append(results['specific_heat'])
    
    plot_temperature_enengy(temperatures, internal_energies, specific_heats)
    
    
