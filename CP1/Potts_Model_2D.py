import numpy as np
# from numba import jit
import time

class PottsModel2D:
    def __init__(self, N, q, J=1.0, h=0.0):
        self.N = N
        self.q = q
        self.J = J
        self.h = h
        self.lattice = np.random.randint(1, q+1, size=(N, N))
        self.k_beta = 1.0 
        self.beta = 1.0
        
    def set_beta(self, T):
        self.beta = 1.0 /(self.k_beta * T)
    
    def wolff_flip(self):

        N, q, beta, J = self.N, self.q, self.beta, self.J
        
        i, j = np.random.randint(0, N, 2)
        old_spin = self.lattice[i, j]
        
        current_set = set()
        boundary = set([(i, j)])
        current_set.add((i, j))
        
        while boundary:
            current_point = boundary.pop()
            i, j = current_point
            
            neighbors = [((i+1)%N, j), ((i-1)%N, j), (i, (j+1)%N), (i, (j-1)%N)]
            
            for nb in neighbors:
                ni, nj = nb
                if nb not in current_set and self.lattice[ni, nj] == old_spin:
                    if np.random.random() <  1.0 - np.exp(-beta * J):
                        current_set.add(nb)
                        boundary.add(nb)
        
        new_spin = np.random.randint(1, q+1)
        while new_spin == old_spin:
            new_spin = np.random.randint(1, q+1)
        
        for site in current_set:
            i, j = site
            self.lattice[i, j] = new_spin
    
    def metropolis_flip(self):
        N = self.N
        for _ in range(N*N):
            i, j = np.random.randint(0, N, 2)
            old_spin = self.lattice[i, j]
            new_spin = np.random.randint(1, self.q+1)
            
            delta_H = 0
            
            neighbors = [((i+1)%N, j), ((i-1)%N, j), (i, (j+1)%N), (i, (j-1)%N)]
            for nb in neighbors:
                ni, nj = nb
                if self.lattice[ni, nj] == old_spin:
                    delta_H += self.J
                if self.lattice[ni, nj] == new_spin:
                    delta_H -= self.J
            
            delta_H += self.h * (old_spin - new_spin)
            
            if delta_H <= 0 or np.random.random() < np.exp(-self.beta * delta_H):
                self.lattice[i, j] = new_spin
    
    def func_H(self):
        """Hamiltonian量"""
        N = self.N
        energy = 0.0
        
        for i in range(N):
            for j in range(N):
                neighbors = [((i+1)%N, j), (i, (j+1)%N)]
                for nb in neighbors:
                    ni, nj = nb
                    if self.lattice[i, j] == self.lattice[ni, nj]:
                        energy -= self.J
                energy -= self.h * self.lattice[i, j]
        
        return energy
    
    def spin_sum(self):
        return np.sum(self.lattice)
    
    def correlation_compute(self, corr_k):
        corr_gamma_ = np.zeros(len(corr_k))
        N = self.N
        for id in range(0,len(corr_k)):
            k = corr_k[id]
            for i in range(N):
                for j in range(N):
                    corr_gamma_[id] += self.lattice[i,j] * (self.lattice[(i+k)%N, j] + self.lattice[i, (j+k)%N])
        return corr_gamma_



def temperature_scheduler(iter, n_iter, T_min, T_max, n_decay):
    """stable cos decay"""
    if iter <= n_iter - n_decay:
        return  T_max
    else:
        return T_min + 0.5 * (T_max - T_min) * (1 + np.cos(np.pi * (iter - (n_iter - n_decay)) / n_decay))


def MCMC_simulation(N, q, T, h, n_tempering, n_measure, n_step=5,
                            mes_energy = True, mes_manetization = False, corr_k = []):

    results = {
        'temperature': T,
        'internal_energy': 0,
        'specific_heat': 0,
        'magnetization': 0,
        'corr_k': corr_k,
        'corr_gamma': np.zeros(len(corr_k))
    }
    
    print(f"MCMC模拟, q={q}, T={T:.3f}, 测量共{n_measure}次")
    
    model = PottsModel2D(N, q, h = h)
    
    # 热化过程
    n_decay = round(n_tempering * 0.2) 
    T_min, T_max = T, max(2*T, 5)
    for step in range(n_tempering):

        T_tempering = temperature_scheduler(step, n_tempering, T_min, T_max, n_decay)
        model.set_temperature(T_tempering)

        model.wolff_step()
        if model.h != 0:
            model.metropolis_step()
    
    
    # 测量过程
    model.set_temperature(T)

    energies = np.zeros(n_measure)
    manetizations = 0
    if len(corr_k)>0:
        lattice_sum = np.zeros((N, N))
        multiple_dis_k = np.zeros(len(corr_k))

    
    for step in range(n_measure):
        
        for _ in range(n_step):
            model.wolff_flip()
            if model.h != 0:
                model.metropolis_flip()
        
        # 测量物理量
        energy = 0
        if mes_energy:
            energy = model.func_H()
            energies[step] = energy
        if mes_manetization:
            manetizations += model.spin_sum()
        if len(corr_k)>0:
            multiple_dis_k += model.correlation_compute(corr_k)
            lattice_sum += model.lattice

        if step % (n_measure/10) == 0:
            if mes_energy:
                print(f"测量过程运行至第{step}步, Hamilton量={energy}")
            else:
                print(f"测量过程运行至第{step}步")
        


    # 计算统计量
    internal_energy = np.mean(energies) / (N**2)
    specific_heat = np.cov(energies) * model.k_beta * model.beta**2 / (N**2)
    manetizations /= n_measure * N**2
    if len(corr_k)>0:
        multiple_dis_k /= 2 *n_measure * N**2
        lattice_sum /= 2 *n_measure**2 * N**2
        for id in range(0,len(corr_k)):
            k = corr_k[id]
            for i in range(N):
                for j in range(N):
                    multiple_dis_k[id] -= lattice_sum[i,j] * (lattice_sum[(i+k)%N, j] + lattice_sum[i, (j+k)%N])
    
    results['internal_energy'] = internal_energy
    results['specific_heat'] = specific_heat
    results['magnetization'] = manetizations
    results['corr_gamma'] = multiple_dis_k
    
    return results

