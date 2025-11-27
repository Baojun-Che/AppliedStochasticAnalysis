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
    data = np.loadtxt('results\q=3,N=100.txt', skiprows=1)
    fig_name = "uc-T,q=3,small"
    plot_temperature_enengy(data[:, 0], data[:, 1], data[:, 2], fig_name)
