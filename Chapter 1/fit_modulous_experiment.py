
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 22:41:09 2022

@author: gopalsamy
"""

## Tuning of initital guess and other parameters required for better fit

######################################################################################################
########## PROGRAM TO FIT DYNAMIC MODULOUS  WITH GKV  using Least Squares method #############
######################################################################################################

## Data for dynamic modulous obtained from Ibishola's thesis for sand bitumen 0/2
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit, least_squares
import pandas as pd

import os
import sys

import matplotlib


font = {'size'   : 16}

matplotlib.rc('font', **font)


##################################################
#### PART 1: Load the data for frequency and temperatures ######
##################################################

"""
file_path =  r'.\theta_frequency.csv'

theta =[]   ## temperature
w =[]    ## frequency
with open(file_path) as f:
    for row in f: # read a row as {column1: value1, column2: value2,...}
            
            a,b = row.split(';')                # based on column name k
            try:
                theta.append(float(a))
                w.append(float(b))
                
            except:
                print('warning:could not convert string to float: '+ a+' , '+b)

w = np.array(w)
theta = np.array(theta)
"""


data = pd.read_excel(r'.\1.xlsx') 

### experimetnal data
f_ex = np.array(list(data['actual frequency']))  
f_r_ex  = np.array(list(data['reduced f']))
T_ex = np.array(list(data['T']))
E_ex = np.array(list(data['E*']))

w = f_ex
theta = T_ex

#T_ex = np.array(T_ex)

n_th = len(theta)

data1 = pd.read_excel(r'.\2.xlsx') 


E1 = list(data1['E1'])
E2 = list(data1['E2'])
E1 = [i/1e6 for i in E1]
E2 = [i/1e6 for i in E2]




###################################################
#### PART 2: Comlpex modulus #############
#########################################################

def HS(T, w):
    ## Complex moudlous for HS model ; parameters obtained from Ibishola's thesis (provided excel sheets )
    a0 = 3.285 ; a1 = -0.40385; a2 = 0.0016562;
    E0 = 52.1947; E_inf = 25067.7;
    delta = 2.8947;
    h = 0.52755;
    k = 0.1817;
    
    def alpha(T):
        return np.exp(a0+a1*T+a2*T**2)
    
    den = 1+ delta*(1j*w*alpha(T))**(-k) +(1j*w*alpha(T))**(-h)  
    E_HS = ((E_inf - E0)/(den)) + E0
    
    return E_HS



def E_star_KV(E_kv,tau_kv,T,w):
    ## Complex moudlous for GKV model  
    
    def shift_factor(T):
        C1 = 32.6331; C2 = 214.5203;
        T_ref = 15
        al = 10**(-C1*(T-T_ref)/(C2+T-T_ref))
        return al
    
    
    
    sf = shift_factor(T)
    
    E_s = 1/E_kv[0]
    for i in range(len(tau_kv)):
        E_s += 1/(E_kv[i+1]*(1+1j*tau_kv[i]*sf*w))
    return 1/E_s


def E_star_GM(E,tau,T,w):
    ## Complex moudlous for GM model  
    
    def shift_factor(T):
        C1 = 32.6331; C2 = 214.5203;
        T_ref = 15
        al = 10**(-C1*(T-T_ref)/(C2+T-T_ref))
        #a0 = 3.285 ; a1 = -0.40385; a2 = 0.0016562;
        #return np.exp(a0+a1*T+a2*T**2)
        return al
    
    
    
    sf = shift_factor(T)
    
    E_s = E[0]
    for i in range(len(tau_kv)):
        E_s += (E[i+1]*(sf*tau[i]*w)**2 + 1j*(sf*tau[i]*w*E[i+1]))/(1+(sf*tau[i]*w)**2)
    return E_s


def obj_func(x, n_kv, T, w,  ind = None , penalty = 10):
    ## (function that provides the residual)
    ## objective function to minimize to fit the paramters for the GKV/GM model
    ### n_kv : number of springs in GKV/GM model
    ### n_kv - number of units to fit
    ### w - numpy array of frequency ranges to fit
    ## y_data - existing data of complex modulous to fit
    ### x[:n_kv] - spring constants  E_i
    ### x[n_kv:] - retardation times of dashpots     (  x[n_kv:2*n_kv-1] )
    weights = 1
    weights1 = 2 ## weights for imaginary part 
    weights3 = 2 ## weights for modulous |E*|
    #hs = HS(T,w)
    #y_data = np.concatenate((hs.real, hs.imag))
    y_data = np.concatenate((np.array(E1), np.array(E2)))
    
    #n= len(hs)
    
    E_kv = x[:n_kv]
    tau_kv = x[n_kv:]
    
    
    ## 
    E_s = E_star_GM(E_kv, tau_kv, T,w)
    n = len(E_s)
    y_kv = np.concatenate((E_s.real, E_s.imag))
    
    
    ## calibration from experimetnal
    #E_s_cal = E_star_KV(E_kv, tau_kv, T_ex,f_ex)
    
    
    
    if ind is not None:
        y_kv_ind = np.zeros_like(y_kv);
        y_data_ind = np.zeros_like(y_data);
        y_kv_ind[ind] = y_kv[ind]
        y_data_ind[ind] = y_data[ind]
        ob_fun =   (y_kv -  y_data)**2  + penalty*((y_kv_ind - y_data_ind))**2
    else:
        #ob_fun =   weights*(y_kv[:n] -  y_data[:n])**2 +  weights1*(y_kv[n:] -  y_data[n:])**2
        ob_fun = np.concatenate(( weights*(y_kv[:n] -  y_data[:n]),weights1*(y_kv[n:] -  y_data[n:])),axis=0)
        ob_fun = np.concatenate((ob_fun,weights3*(np.abs(E_s)-E_ex)))
        #ob_fun = np.abs(E_s_cal)-E_ex
        #ob_fun[-3:]*=100
        #ob_fun[n:5*n]*=2
        
        #mod_ob_fun = np.abs(E_s) - np.abs(hs)
        #mod_ob_fun[-6:]*= weights3
        #mod_ob_fun[-8:-6]*= weights3*2
        #ob_fun = np.concatenate((ob_fun, mod_ob_fun),axis=0)
    return ob_fun

##########################################
#######  FIT PARAMTERS OF GKV MODEL ######
##########################################





## number of units to fit including free spring
n_kv = 10

## initial guess  for stiffnesses and retardation/relaxation times for GKV/GM model
E_kv =[1e3]*n_kv
tau_kv = [1e-5]*(n_kv-1)
## unknowns in an array
x0 = np.array(E_kv+tau_kv)



## define bounds for the paramters  (E_kv and tau_kv)

lb=  [1]*n_kv + [1e-10]*(n_kv-1)
ub = [1e12]*n_kv + [1e5]*(n_kv-1)


## Choice 1 (without penalty)

#popt = least_squares(obj_func, x0, args = (n_kv,theta,w), bounds=((lb,ub)))



## Choice 2 (with penalty to force the residaul to be minimum at selected points given by ind)
ind = [0]+[len(w)-i for i in range(2)]    ## indices where the residual is minimized with additionla effort
penalty = 10          ## weight/penalty factor for minimization of residaul at ind
#popt = least_squares(obj_func, x0, args = (n_kv,theta,w, ind, penalty), bounds=((lb,ub)))





## Choice 2.1 (same as 2 with some strict error tolerance) (doesnt have notcible effect)
ftol =1e-15
xtol = 1e-15
gtol = 1e-15
#ind = [0,1]    ## indices where the residual is minimized with additionla effort
#penalty = 10          ## weight/penalty factor for minimization of residaul at ind
popt = least_squares(obj_func, x0, args = (n_kv,theta,w), bounds=((lb,ub)), ftol=ftol,
                     xtol =xtol, gtol = gtol)


x= popt.x

#x[0] = y_data[-1]


E_kv = x[:n_kv]
tau_kv = x[n_kv:]

##################################################################################
################## COMPARE THE FIT WITH INITAIL COMPLEX MODULOUS #################
##################################################################################

### IMP : Multiply the final values of E_kv by 1e6 (in the considerd case as y_data is in MPa)

fit_mod = (E_star_GM(E_kv, tau_kv, theta,w))
#fit_hs = HS(theta,w)

plt.figure(1)
#plt.plot(fit_hs.real, fit_hs.imag,'*')
plt.plot(fit_mod.real,fit_mod.imag,'ok', markersize=3)
try:
    for i in range(n_th):
        plt.plot(E1[6*i:6*(i+1)], E2[6*i:6*(i+1)])
except:
    print()
#plt.loglog(w,E_star_GM(w)/1e6,'--')
#plt.semilogx(w, np.abs(poisson_star(w)))
#plt.legend(['Huet-Sayegh', 'GKV fit'])
plt.legend(['fit',r'$-10^oC$',r'$0^oC$',r'$10^oC$',r'$15^oC$',r'$20^oC$',r'$30^oC$',r'$40^oC$'])
plt.xlabel('E1 (MPa)')
plt.ylabel('E2 (MPa)')

"""
for i in range(int(len(w)/6)):
    plt.loglog(w[6*i:6*(i+1)], np.abs(fit_hs[6*i:6*(i+1)]),'*')
    plt.loglog(w[6*i:6*(i+1)], np.abs(fit_mod[6*i:6*(i+1)]))
"""


## perform DMA (Dynamic mexhanical analysis to check the validate the parameters)
#E_kv = [1e6*i for i in E_kv]

def shift_factor(T):
        C1 = 32.6331; C2 = 214.5203;
        T_ref = 15
        al = 10**(-C1*(T-T_ref)/(C2+T-T_ref))
        return al
    
def find_sigma(eps,theta, dt, eps_i_n):
    
    
    sf = shift_factor(theta)
    tau = [sf*i for i in tau_kv]
    
    H = 1
    eps_int = 0
    for i in range(len(tau_kv)):
        H += (dt/(dt+tau[i]))*E_kv[0]*(E_kv[i+1])**-1
        eps_int += (tau[i]/(tau[i]+dt))*eps_i_n[i]
    H = (H**-1)*E_kv[0]
    sigma = H*(eps-eps_int)
    return sigma

def find_eps_i(sigma,theta,eps_i_n,dt):
    sf = shift_factor(theta)
    tau = [sf*i for i in tau_kv]
    eps_i = []
    for i in range(len(tau)):
        tmp = sigma +((E_kv[i+1]*tau[i]/dt)*eps_i_n[i])
        tmp = tmp * dt/(E_kv[i+1]*(dt) + E_kv[i+1]*tau[i])
        eps_i.append(tmp)
    return eps_i
        
def DMA(theta,w,a=10, e0= 1):
    ## a: a time steps per cycle
    dt = 1/(a*w)
    T = 1  ## final time
    nt = int(T/dt)
    t = np.linspace(0,T,nt)
    eps = e0*np.sin(w*t)
    
    eps_i_n= [0]*len(tau_kv)
    sigma = [0]
    for i in range(len(t[1:])):
        sigma.append(find_sigma(eps[i+1], theta, dt, eps_i_n))
        eps_i_n = find_eps_i(sigma[-1], theta, eps_i_n, dt)
    
    return max(sigma)/e0

"""

## numerical Dynamic Mechanical Analysis (DMA) test for GKV model to validate the fitting process


E_dma = []
fit_hs = HS(theta,w)
for i in range(len(w)):
    print("performing for i ="+ str(i))
    
    E_dma.append(DMA(theta[i],w[i]))
    
a_the=shift_factor(theta)

w_r = a_the*w

plt.figure(2)
plt.plot(np.log(w_r),np.abs(fit_hs),'*')
plt.plot(np.log(w_r),E_dma,'o')    
"""

a_the=shift_factor(theta)

w_r = a_the*w

plt.figure(2)
#plt.semilogx(w_r,np.abs(fit_hs),'*')
plt.semilogx(w_r,np.abs(fit_mod),'ok',markersize=3)  
try:
    for i in range(n_th):
        plt.semilogx(f_r_ex[6*i:6*(i+1)], E_ex[6*i:6*(i+1)])
except:
    print()
plt.ylabel('|E*| (MPa)')
plt.xlabel('frequency (Hz)')
plt.legend(['fit',r'$-10^oC$',r'$0^oC$',r'$10^oC$',r'$15^oC$',r'$20^oC$',r'$30^oC$',r'$40^oC$'])
