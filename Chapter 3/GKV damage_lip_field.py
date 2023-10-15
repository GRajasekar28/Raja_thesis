# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 19:42:02 2023

@author: gopalsamy
"""

## 1D FE code for damage in a homogenous bar using the softening GKV model
## use of an implicit time discretisation

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds
import scipy



from abc import ABC, abstractmethod



class HPoly():
    def __init__(self, coef=[0,2.], name='hpoly'):
        self._poly = np.polynomial.Polynomial(coef)
        self._deri1, self._deri2  = (self._poly.deriv(),self._poly.deriv(2))
    def __call__(self,d): return self._poly(d)
    def jac(self, d):     return self._deri1(d)
    def hess(self,d):     return self._deri2(d)

class HRPoly():
    def __init__(self, coefn, coefd, name='hrpoly'):
        self._num, self._den =  (np.polynomial.Polynomial(coefn), np.polynomial.Polynomial(coefd))
        self._deriv1_num, self._deriv2_num = (self._num.deriv(),self._num.deriv(2))
        self._deriv1_den, self._deriv2_den = (self._den.deriv(),self._den.deriv(2))
    def __call__(self,d) : 
        return self._num(d)/self._den(d)
    def jac(self, d) : 
        f, df = (self._num(d), self._deriv1_num(d))
        g, dg = (self._den(d), self._deriv1_den(d))
        return (df*g-f*dg)/g**2
    def hess(self, d) : 
        f, df, d2f = (self._num(d), self._deriv1_num(d), self._deriv2_num(d) )
        g, dg, d2g = (self._den(d), self._deriv1_den(d), self._deriv2_den(d))
        return ((d2f*g-f*d2g)*g -2.*(df*g-f*dg)*dg)/g**3    

def HLinear():      return HPoly(coef=[0.,2.],name='h(d) = 2d') # Gc = 2YcLc
def HQuadratic():   return HPoly(coef=[0,2,3],name='h(d) = 2d+3d²') #Gc = 4yc lc ??
def HQuadratic2():  return HPoly(coef=[0,0,2],name='h(d) = 2d²')   #Gc = 2/3 Yclc damage imediat
# to get a convex function we need, 0<lm<0.5
# for a given Gc, lc, lm = 2.Yc*lc/Gc 
# for example, Yc = 1, Gc =1, 0.1<lc < 0.25 
# -> lc =0.2 to be safe for the coarser mesh. then lm = 0.4
#    now, refine h<- h/2, we can take lc = 0.1 and  lm = 0.2
def HCohesive(lm) : return HRPoly([0.,2.,-1.],[1.,-2,(1.+2*lm),-2.*lm, lm**2], name = 'cohesive l='+str(lm)) 

def GLinear():      return HPoly(coef=[1.,-1.],name= 'g(d) = 1-d')
def GQuadratic():   return HPoly(coef=[1.,-2.,1.],name= 'g(d) = (1-d)^2') #1.-2*d+d^2
def G_Const_1(): return HPoly(coef=[1.],name='g(d)=1')
#typical value for eta = 0.5h/l
def GO3Eta(eta):   return HPoly(coef=[1.,-2., 1.+eta, -eta],name= '(1.-d)^2 +eta*(1.-d)*d**2') #1.-2d+d^2 +eta*d^2 - eta*d^3
def GO4Eta(eta):   return HPoly(coef=[1.,-2., 1.,      eta, -eta],name= 'O4_LE') 





eta = 0.
g1d = GQuadratic()
g2d = G_Const_1()
hd = HQuadratic()



## properties

    ## Material
E = [ 20555.40787441, 135374.72188481,  84505.60095753,  11142.58265657,
        56337.78789284,   3195.23653526, 181136.00646417,    327.13513483,
        21958.9894283 ,  34646.72374612]
tau = [1.98638904e-05, 2.35339793e-04, 1.12355128e+00, 2.58717530e-03,
       1.34975082e+01, 1.75757042e-06, 5.36062114e+02, 1.59580031e-01,
       2.08841625e-02]
E = [a*1e6 for a in E]
n = len(tau)
"""
n=1
E= [14000, 5600]
tau = [2.59e-1]
E = [a*1e6 for a in E]
"""

Y_c = 500
    ## Geometric and time domain 
    
L = 1.  # m
T = 400 # s

l = .5

### Applied montonous displacement rate at right edge
U_dot = 1e-3 ## m/s

## Discretisation

## number of nodes
N = 40
x = np.linspace(0,L,N)
dx = L/(N-1)

## number of time steps
dt = 1e-3
n_t = int(T/dt +1)
t = np.linspace(0,T,n_t)

## matrices to store variables

## row 'i' correpsond to time step 'i'
## coloumns correspond to nodal/elemental values
u = np.zeros([n_t, N])
eps = np.zeros([n_t, N-1])
sigma = np.zeros([n_t,N-1])
d = np.zeros([n_t, N-1])


## matrices to store internal strain variables (IMP:eps_0 not stored;; eps_i := eps_{i})

epsilon_i =[]
for n_u in range(1,n+1):
    epsilon_i.append(np.zeros([n_t,N-1]))


## Initial conditions imposed already in u field by np.zeros
## applying displacement BCDs at right end of bar (left end fixed)
u[:,-1] = U_dot*t



## totoal energies at each time step
fe = [0]
ve = [0]
de = [0]
wi = [0]


def stiffness(E,tau,n,N,dx,dt,d):
    K = np.zeros([N,N])
    for i in range(N-1):
        k = modified_stifness(E, tau,n, dt, d[i])
        K[i,i] += (k/dx)
        K[i+1,i+1] += (k/dx)
        K[i,i+1] = -(k/dx)
        K[i+1,i] = -(k/dx)
    return K



def modified_stifness(E,tau,n,dt,d_i):
    g1 = g1d(d_i)
    g2 = g2d(d_i)
    k = 1/(E[0])
    for i in range(1,n+1):
        k += g1*dt/(g1*dt*E[i]+g2*tau[i-1]*E[i])
    
    return g1/k
    


###IMP: BCD yet to be enforced in rhs_1
def rhs_1(eps_i, E,tau,dx, dt, t_np1,n,d):
    ## t_n - previous time step index
    ## eps_i - list of matrices of internal variables (eps_1, eps_2.. excluding eps_0)
    g1 = g1d(d)
    g2 = g2d(d)
    f= np.zeros(N)
    for i in range(N-1):
        k = modified_stifness(E, tau,n, dt, d[i])
        for j in range(1,n+1):
            temp = g2[i]*tau[j-1]/(g1[i]*dt+g2[i]*tau[j-1])
            f[i] += k*temp*eps_i[j-1][t_np1-1,i]*(-1)
            f[i+1] += k*temp*eps_i[j-1][t_np1-1,i]*(1)
    return f


def psi_0(E, epsilon, epsilon_i, t_np1, dx):
    ## function to calculate undamaged free energies at each element for a 
    ## given time step index t_np1
    psi_e = np.zeros(N-1)  ## free enegy of element
    ## finding eps_0 = epsilon - sum(eps_i) ;; i =1 to n;    i != 0
    eps_0 = np.copy(epsilon[t_np1,:])
    for i in range(1,n+1):
        eps_0 -= epsilon_i[i-1][t_np1,:]
        psi_e +=  .5*E[i]*epsilon_i[i-1][t_np1,:]**2
    psi_e += .5*E[0]*eps_0**2  ##not multiplied  dx as the element size are same
    return psi_e

## function to calulcate undamaged viscous dissipation potential*dt at each element 
def func_psi_v(N,n,E,tau,epsilon_i,dt, t_np1):
    psi_v = np.zeros(N-1)
    for i in range(1,n+1): 
        psi_v += .5*tau[i-1]*E[i]*((epsilon_i[i-1][t_np1,:]-epsilon_i[i-1][t_np1-1,:])/dt)**2
        ## not multiplied by dx as element size are same
    return psi_v*dt


## no regularisation 
## incremental potential to be minimized for finding damage 
## entire expression multiplied by dx for numerically approximated f;
def f_d_1(d,psi_0,psi_v,Y_c):
    g1 = g1d(d)
    g2 = g2d(d)
    h = hd(d)
    return sum(g1[:] * psi_0[:] + g2[:] * psi_v[:] + Y_c * h[:])





k = modified_stifness(E, tau, n, dt, 0)

## stiffness tensor
K = stiffness(E, tau, n, N, dx, dt, np.zeros(N-1))

## index to store the indices correpsonding to beginning of damage phase
ind_d = 1
i = 0
## time iteration
while 1:

    i+=1

    
    f = rhs_1(epsilon_i, E, tau, dx, dt, i, n, np.zeros(N-1))
    ## enforcing BCD's (0 on left and U_dot on right)
    f[1] -= K[1,0]*u[i,0];    f[-2] -= K[-2,-1]*u[i,-1]
    
    ## finding u at i th time step (leaving out 1st and last rows and coloumns due to imposed BCD's)
    U = np.linalg.solve(K[1:-1,1:-1], f[1:-1])
    
    u[i,1:-1] = np.copy(U)
    
    eps[i,:] = np.diff(u[i,:])/dx
    
    
    
    ## find sigma at i th time step
    for i1 in range(1,n+1):
        temp = tau[i1-1]/(tau[i1-1]+dt)
        for i2 in range(N-1):
            sigma[i,i2] += temp*epsilon_i[i1-1][i-1,i2]
    sigma[i,:] = k*(eps[i,:] - sigma[i,:])
    
    
    ## sigma / E[0]  provides eps_0 
    ## update internal strain variables eps_1, eps_2, .. excluding eps_0
    for i1 in range(1,n+1):
        epsilon_i[i1-1][i,:] = (1/(tau[i1-1]+dt))*(dt*sigma[i,:]/E[i1]  + tau[i1-1]*epsilon_i[i1-1][i-1,:])
        
    psi_e = psi_0(E, eps, epsilon_i, i, dx)
    psi_v = func_psi_v(N, n, E, tau, epsilon_i, dt, i)
    
    if (i%100 ==0):
        print(max(psi_e))
    
    if max(psi_e) > Y_c:
        ### removing all changes made in the current time step
        u[i,1:-1] = np.zeros(N-2)
        eps[i,:] = np.zeros(N-1)
        sigma[i,:] = np.zeros(N-1)
        for i1 in range(1,n+1):
            epsilon_i[i1-1][i,:] = np.zeros(N-1)
        ind_d = i
        print('Damage started')
        break
    else:
        ## store respective enrgies
        temp_visc_dissip = np.zeros(N-1)
        for i1 in range(1,n+1):        ## n since only viscous components are accountes
            ## array of viscous strain rates
            temp_eps_i_dot = (epsilon_i[i1-1][i,:]-epsilon_i[i1-1][i-1,:])/dt
            temp_visc_dissip += dt*(tau[i1-1]*E[i1]*temp_eps_i_dot**2)
        
        fe.append(np.sum(psi_e)*dx)
        ve.append(ve[-1]+np.sum(temp_visc_dissip)*dx)
        de.append(0)
        wi.append(wi[-1]+np.mean(sigma[i,:])*U_dot*dt)


"""Damage phase"""
## starting from ind_d

max_iter = 1000

i = ind_d-1
#d[i, int(N//2)]  = .5
## time iteration
while 1:
    
    i += 1
    ## Initial values of non-linear iteration for coupling b/w visco-elasticity 
    ## and damage problem
    d_m = np.copy(d[i-1,:]);  u_m = np.copy(u[i-1,:]);
    ## inital trigger for localsiation
    
    d_m[int(N//2)] += 1e-3
    ###variable to store values at  m+1 th non-linear iteration
    u_mp1 = np.zeros(N)
    ##BCD's on u_m and u_mp1
    u_m[0] = u[i,0]; u_m[-1] = u[i,-1];  
    u_mp1[0] = u[i,0]; u_mp1[-1] = u[i,-1]; 
    bounds = Bounds(d[i-1,:], np.ones(N-1))   
    for j in range(max_iter+1):
        K = stiffness(E, tau, n, N, dx, dt, d_m)
        f = rhs_1(epsilon_i, E, tau, dx, dt, i, n, d_m)
        ## enforcing BCD's (0 on left and U_dot on right)
        f[1] -= K[1,0]*u_m[0];    f[-2] -= K[-2,-1]*u_m[-1]
        
        ## finding u at i th time step (leaving out 1st and last rows and coloumns due to imposed BCD's)
        U = np.linalg.solve(K[1:-1,1:-1], f[1:-1])
        
        u_mp1[1:-1] = np.copy(U)
        
        eps[i,:] = np.diff(u_mp1)/dx
        
        g1 = g1d(d_m) ; g2 = g2d(d_m)
        
        ##  sigma reinitalized to zero within each non-linar iteration due to += operator invloved with it
        sigma[i,:] = np.zeros(N-1)
        for i1 in range(N-1):
            k = modified_stifness(E, tau, n, dt, d_m[i1])
            for i2 in range(1,n+1):
                temp = g2[i1]*tau[i2-1]/(g2[i1]*tau[i2-1]+g1[i1]*dt)
                sigma[i,i1] += temp*epsilon_i[i2-1][i-1,i1]
            sigma[i,i1] = k*(eps[i,i1] - sigma[i,i1])
        
        
        for i1 in range(1,n+1):
            epsilon_i[i1-1][i,:] = (1/(g2*tau[i1-1]+g1*dt))*(dt*sigma[i,:]/(E[i1])  + g2*tau[i1-1]*epsilon_i[i1-1][i-1,:])
        
        
        psi_e = psi_0(E, eps, epsilon_i, i, dx)
        psi_v = func_psi_v(N, n, E, tau, epsilon_i, dt, i)
        
        
        Aslope=scipy.sparse.eye(N-2,N-1) - scipy.sparse.eye(N-2,N-1,1)
        ## for Lipschitz d
        slopeconstrain = scipy.optimize.LinearConstraint(Aslope, -dx/l *np.ones(N-2), dx/l*np.ones(N-2) )
        ## for homogeneous d
        slopeconstrain = scipy.optimize.LinearConstraint(Aslope, -1e-10 *np.ones(N-2), 1e-10*np.ones(N-2) )
        lipschitz_constrain = slopeconstrain
        
        
        res = minimize(f_d_1,d_m,args=(psi_e,psi_v,Y_c) ,options={'ftol':1e-14, 'disp':False},bounds=bounds,constraints= lipschitz_constrain)
        #if not res.success:
        #    print("Minimization failed")
        d_mp1 = np.copy(res.x)
        
        
        
        err_u = dx*np.linalg.norm(u_mp1-u_m);
        err_d = dx*np.linalg.norm(d_mp1-d_m);
        #print(i,j, err_u,err_d, max(psi_e))
        
        if ((err_u <= 1e-11) and (err_d <= 1e-11)):
            print(i,j, err_u,err_d, min(d_mp1), max(d_mp1))
            u[i,:] = np.copy(u_mp1);
            d[i,:] = np.copy(d_mp1)
            temp_visc_dissip = np.zeros(N-1)
            for i1 in range(1,n+1):        ## n since only viscous components are accountes
                ## array of viscous strain rates
                temp_eps_i_dot = (epsilon_i[i1-1][i,:]-epsilon_i[i1-1][i-1,:])/dt
                temp_visc_dissip += dt*(tau[i1-1]*E[i1]*temp_eps_i_dot**2)
            
            fe.append(np.sum(g1d(d_mp1)*psi_e)*dx)
            ve.append(ve[-1]+np.sum(g2d(d_mp1)*temp_visc_dissip)*dx)
            de.append(Y_c*np.sum(hd(d_mp1))*dx)
            wi.append(wi[-1]+np.mean(sigma[i,:])*U_dot*dt)
            break
        elif j == max_iter:
            raise Exception('Maximum iteration reached in the non-linear iteration')
        else:
            u_m = np.copy(u_mp1)
            d_m = np.copy(d_mp1)
    if max(d[i,:]) > .99:
        raise Exception('Maximum set damage value reached')
        

## for stress-strain plot
plt.plot((U_dot/L)*t[:i], [np.mean(sigma[j,:]) for j in range(len(t[:i]))])
plt.xlabel('strain')
plt.ylabel('stress (Pa)')



## for energy plots

plt.figure(2)
plt.plot((U_dot/L)*t[:i],wi)
plt.plot((U_dot/L)*t[:i],fe)
plt.plot((U_dot/L)*t[:i],ve)
plt.plot((U_dot/L)*t[:i],de)
plt.plot((U_dot/L)*t[:i],np.array(fe)+np.array(ve)+np.array(de))



plt.legend(['Work input','fe','vd','de','fe+vd+de'])
plt.xlabel('time step index')
plt.ylabel(' energy  E (J/m2)')
            
            
    
        
