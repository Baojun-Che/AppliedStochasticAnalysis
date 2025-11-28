from Potts_Model_2D import mcmc_without_external_field
import matplotlib.pyplot as plt
import numpy as np
import os

if __name__ == "__main__":
    
    N = 100
    q = 3
    ###  q=3, T以细尺度变化, 加密到N=500
    T = 0.5
    _, energies, lattice = mcmc_without_external_field(N, q, T, n_tempering=200, n_measure=200, RATE=1, n_step=50, get_energy=True)   
    np.save('lattice/N=100,T=0.50.npy', lattice)