# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:49:39 2023

@author: gopalsamy
"""



## program writen to see if for equilent gkv and gm parameters, 
## same local solution is produced  (1d bar - homogeneous solution)


import numpy as np
import matplotlib.pyplot as plt


def bar_analytical_sol_gkv(E, tau,eps_dot, ni,eps_max,dt=1e-2):
    
    ## list to store epsilon_i  
    ## epsilon_i[j][k] = epsilon_k for time step j
    epsilon_i = [[0]*ni]     ## internal strain  epsilon_0, epsilon_1,,...
    epsilon = [0]          ## total strain
    sigma =[0]             ## sigma (constant along the bar)
    t = [0]                ## time 
    
    while epsilon[-1] < eps_max:
        
        k1 = 0; k2=0
        for j in range(ni-1):
            k1 += tau[j] * epsilon_i[-1][j+1]/(tau[j] + dt)
            k2 += dt/((tau[j]+dt)*E[j+1])
        sigma_tdt =  ((eps_dot*(t[-1]+dt)) -(k1))/(1/E[0] + k2)
        epsilon_tdt = eps_dot*(t[-1]+dt)
        ## temporary list to store the internal variable epsilon_i at time step i+1
        tmp_eps_i = []
        tmp_eps_i.append(sigma_tdt/E[0])
        for j in range(ni-1):
            k3 = (tau[j]*epsilon_i[-1][j+1]/(tau[j] + dt))  + (dt*sigma_tdt/((tau[j] + dt)*E[j+1]))
            tmp_eps_i.append(k3)
            
        epsilon.append(epsilon_tdt)  
        t.append(t[-1]+dt)
        sigma.append(sigma_tdt)    
        epsilon_i.append(tmp_eps_i)
        
    return epsilon, sigma

def bar_analytical_sol_gm(E, tau, eps_dot, ni, eps_max, dt = 1e-3):
    sigma = [0]
    epsilon = [0]
    sigma_t = 0
    t = 0
    while epsilon[-1] < eps_max:
        E_t = E[0] 
        for i in range(ni-1):
            E_t += E[i+1]*np.exp(-t/tau[i])
        sigma_t += eps_dot *  E_t * dt
        
        t += dt
        sigma.append(sigma_t)
        
        epsilon.append(eps_dot*t)
        
    return epsilon,sigma




E_gmm= [1.33e2, 3.81e2,4.61e2,9.38e2,1.22e3,
        9.54e2,5.07e2,2.70e2,2.69e3, 3.28e3,
        2.44e3,2.94e3,2.84e3,1.59e3,2.53e3,1.63e3]

tau_gmm = [3.02e4,3.29e3,2e3,4.12e2,1.69e2,
           1.28e2,1.2e2,3.17e1,4.74,7.88e-1,1.52e-1,
           1.23e-2,1.09e-3,5.33e-5,2.25e-7]


E_gmm = [1.07e2, 2.48e2, 5.54e2,6.67e2,1.19e3,
       3.81e2,1.57e2,2.74e3,3.74e2,2.26e3,1.88e3,
       2.93e3,2.99e3,3.04e3,2.49e3,2.20e3]                      ##  MPa

#E_gmm = [i*1e6 for i in E_gmm]    ## Pa

tau_gmm = [8.21e1,1.19e1,3.03,9.35e-1,6.02e-1,
           5.77e-1,1.39e-1,2.34e-2,2.29e-2,6.35e-3,
           1.33e-3,1.65e-4,1.33e-5,5.16e-7,4.13e-9]

E_gmm = [ 275.85407523, 2152.22939477, 2243.32000377, 2384.95964582,
       2711.76116699, 2673.24436825, 2769.93083808, 2100.50288175,
       1224.30984935, 2019.13331597]

tau_gmm =[1.57627333e-06, 1.74634962e-05, 1.13884723e-01, 1.57316509e-02,
       1.96913752e-04, 2.06974102e-03, 7.16389390e-01, 9.08833781e+01,
       5.90861097e+00]

n_gmm = len(E_gmm)

E_kv =[2.48039994e+04, 1.44884216e+03, 1.75235504e+02, 5.59975404e+03,
       3.53003103e+05, 8.07536029e+04, 1.90960402e+05, 1.12055417e+05,
       1.48894840e+04, 2.52996598e+05, 6.35403845e+04, 2.80123475e+04]

tau_kv = [6.08095858e+03, 1.00000000e+05, 4.75082423e+02, 2.40819595e-07,
       1.85547684e-01, 5.98048026e-05, 1.44457131e-02, 5.63167665e+01,
       1.18082841e-03, 9.82333948e-01, 6.93728393e+00]

E_kv = [2.39845097e+04, 9.09466809e+03, 2.32046243e+04, 5.60390911e+04,
       2.42834974e+05, 7.90630834e+04, 3.53528269e+03, 5.37190564e+04,
       1.27094674e+02, 1.82996182e+05, 1.10104824e+05, 1.15282316e+03]

tau_kv = [2.77109889e-01, 3.51831534e-02, 8.32942868e-03, 1.00000000e-08,
       2.03526646e-04, 2.31101743e+00, 1.72112000e-03, 3.06892337e+02,
       6.67101723e-07, 1.61628261e-05, 1.86580281e+01]

E_kv = [ 20555.40787441, 135374.72188481,  84505.60095753,  11142.58265657,
        56337.78789284,   3195.23653526, 181136.00646417,    327.13513483,
        21958.9894283 ,  34646.72374612]
tau_kv = [1.98638904e-05, 2.35339793e-04, 1.12355128e+00, 2.58717530e-03,
       1.34975082e+01, 1.75757042e-06, 5.36062114e+02, 1.59580031e-01,
       2.08841625e-02]

E_kv = [i*1e6 for i in E_kv]
E_gmm = [i*1e6 for i in E_gmm]


n_kv = len(E_kv)

eps_dot = [1e-1, 1e-3]
eps_dot = [1e-1]

for i in eps_dot:
    e1,s1 = bar_analytical_sol_gkv(E_kv, tau_kv, i, n_kv, .1)
    e2,s2 = bar_analytical_sol_gm(E_gmm, tau_gmm, i, n_gmm, .1)
    
    plt.plot(e1,s1)
    plt.plot(e2,s2,'--')
    plt.legend(['GKV', 'GM'])
    plt.xlabel('strain')
    plt.ylabel('stress (Pa)')