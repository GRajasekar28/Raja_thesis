# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 08:21:17 2023

@author: gopalsamy
"""

## program to obtain local solution for GKV/GM model with damage
## use of an explicit time scheme.. hence a smaller dt is better

### This program is used to validate the 1D FE code for a bar with homogenous damage (local solution)

import numpy as np
import pylab as plt
from abc import ABC, abstractmethod
import scipy
from scipy import optimize
import copy


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


class H_TLS_1:
    def __init__(self, alpha =3.5, beta= .85, c=2, name = 'TLS_eq_soft_func'):
        ## c=2 for GQuadrati and c=1 for GLinear
        self.name = name+' alpha ='+str(alpha)+' beta ='+str(beta)+' c ='+str(c)
        self.a = alpha
        self.b = beta
        self.c = c
        
    def __call__(self,d):
        a = self.a
        b= self.b
        c= self.c
        d1 = copy.copy(d)
        try:
            d1[d1>.999] = .999
        except:
            if d1> .999:
                d1 = .999
        return ((c*(1-b*d1)**(1-a))/(b*(a-1))) -c/(b*(a-1))
    def jac(self,d):
        d1 = copy.copy(d)
        c = self.c
        try:
            d1[d1>.999] = .999
        except:
            if d1> .999:
                d1 = .999
        return c*(1-self.b*d1)**(-self.a)
    def hess(self,d):
        d1 = copy.copy(d)
        c = self.c
        try:
            d1[d1>.999] = .999
        except:
            if d1> .999:
                d1 = .999
        return c*self.a *self.b *(1-self.b*d1)**(-self.a-1.)



def GLinear():      return HPoly(coef=[1.,-1.],name= 'g(d) = 1-d')
def GQuadratic():   return HPoly(coef=[1.,-2.,1.],name= 'g(d) = (1-d)^2') #1.-2*d+d^2
def G_Const_1(): return HPoly(coef=[1.],name='g(d)=1')
#typical value for eta = 0.5h/l
def GO3Eta(eta):   return HPoly(coef=[1.,-2., 1.+eta, -eta],name= '(1.-d)^2 +eta*(1.-d)*d**2') #1.-2d+d^2 +eta*d^2 - eta*d^3
def GO4Eta(eta):   return HPoly(coef=[1.,-2., 1.,      eta, -eta],name= 'O4_LE') 





eta = 0.
g1 = GQuadratic()
g2 = G_Const_1()
h = HQuadratic()

#h = H_TLS_1(1.8,.99, 2)




def bar_analytical_sol_gm(E, tau, eps_dot, ni, eps_max, Yc, dt = 1e-1):
    ## ni = len(tau)
    
    def Y(d,eps,eps_i):
        phie = .5*E[0]*eps**2
        phiv = 0
        for i in range(1,ni+1):
            phie += .5*E[i]*(eps-eps_i[i-1][-1])**2
            phiv += .5*dt*tau[i-1]*E[i]*((eps_i[i-1][-1]-eps_i[i-1][-2])/dt)**2
        return -g1.jac(d)*phie  -g2.jac(d)*phiv  - Yc * h.jac(d)


    sigma = [0]
    epsilon = [0]
    eps_i = [[0] for i in range(ni)]   ## eps_0 not stored  in eps_i (since eps_0 = eps)
    d = [0]
    sigma_t = 0
    t = 0
    while epsilon[-1] < eps_max:
        
        ## find eps_i with d_n
        for i in range(ni):
            temp =(dt/(g1(d[-1])*dt + g2(d[-1])*tau[i])) *(g1(d[-1])*epsilon[-1] + g2(d[-1])*tau[i]*eps_i[i][-1]/dt)
            eps_i[i].append(temp)
        
        ## find sigma with d_n
        sigma_t = g1(d[-1])*E[0]*epsilon[-1]
        for i in range(1,ni+1):
            sigma_t += g1(d[-1])* E[i]*(epsilon[-1]- eps_i[i-1][-1])
        sigma.append(sigma_t)
        
        fun = lambda x: Y(x, epsilon[-1],eps_i)
        ## find d_n+1
        if fun(1.)>0.: 
            d.append(1.)
        elif Y(d[-1], epsilon[-1], eps_i) >0.:    
            d.append( scipy.optimize.brentq(fun, d[-1], 1.))
        else:
            d.append(d[-1])
        ##  constraints on d
        ##irreversibility constraint
        if d[-1] < d[-2]:
            d[-1] = d[-2]
        if d[-1]> 1.:
            d[-1] = 1.
            
        if d[-1] > .9:
            break
        
        
        t += dt
        epsilon.append(eps_dot*t)
        
        
        
        
    return epsilon,sigma, d


def bar_analytical_sol_gkv(E, tau, eps_dot, ni, eps_max, Yc, dt = 1e-1):
    ## ni = len(tau)
    
    def Y(d,eps,eps_i):
        phie = 0
        phiv = 0
        eps_0 = eps
        for i in range(1,ni+1):
            phie += .5*E[i]*(eps_i[i-1][-1])**2
            phiv += .5*dt*tau[i-1]*E[i]*((eps_i[i-1][-1]-eps_i[i-1][-2])/dt)**2
            eps_0-= eps_i[i-1][-1]
        phie += .5*E[0]*eps_0**2
        return -g1.jac(d)*phie  -g2.jac(d)*phiv  - Yc * h.jac(d)


    sigma = [0]
    epsilon = [0]
    eps_i = [[0] for i in range(ni)]   ## eps_0 not stored  in eps_i (since eps_0 = eps- sum(eps_i))
    d = [0]
    sigma_t = 0
    t = 0
    while epsilon[-1] < eps_max:
        
        
        
        ## find eps_i with d_n
        for i in range(ni):
            temp =(dt/(g1(d[-1])*dt + g2(d[-1])*tau[i])) *(sigma[-1]/E[i+1] + g2(d[-1])*tau[i]*eps_i[i][-1]/dt)
            eps_i[i].append(temp)
        
        ## find sigma with d_n
        k = 1
        for i in range(ni):
            k += g1(d[-1]) * E[0]*dt/(E[i+1]*(g1(d[-1])*dt+g2(d[-1])*tau[i]))
        sigma_t = g1(d[-1])*E[0]*(eps_dot*(t+dt))
        for i in range(1,ni+1):
            sigma_t -= g1(d[-1])* g2(d[-1])*E[0] *tau[i-1] * eps_i[i-1][-2]/((g1(d[-1])*dt + g2(d[-1])*tau[i-1]))
        sigma.append(sigma_t/k)
        
        fun = lambda x: Y(x, epsilon[-1],eps_i)
        ## find d_n+1
        if fun(1.)>0.: 
            d.append(1.)
        elif Y(d[-1], epsilon[-1], eps_i) >0.:    
            d.append( scipy.optimize.brentq(fun, d[-1], 1.))
        else:
            d.append(d[-1])
        ##  constraints on d
        ##irreversibility constraint
        if d[-1] < d[-2]:
            d[-1] = d[-2]
        if d[-1]> 1.:
            d[-1] = 1.
            
        if d[-1] > .9:
            break
        
        
        t += dt
        epsilon.append(eps_dot*t)
        
        
        
        
    return epsilon,sigma, d

E_gmm = [ 275.85407523, 2152.22939477, 2243.32000377, 2384.95964582,
       2711.76116699, 2673.24436825, 2769.93083808, 2100.50288175,
       1224.30984935, 2019.13331597]

tau_gmm =[1.57627333e-06, 1.74634962e-05, 1.13884723e-01, 1.57316509e-02,
       1.96913752e-04, 2.06974102e-03, 7.16389390e-01, 9.08833781e+01,
       5.90861097e+00]
E_gmm = [i*1e6 for i in E_gmm]    ## Pa


n_gm = len(tau_gmm)



E_kv = [ 20555.40787441, 135374.72188481,  84505.60095753,  11142.58265657,
        56337.78789284,   3195.23653526, 181136.00646417,    327.13513483,
        21958.9894283 ,  34646.72374612]
tau_kv = [1.98638904e-05, 2.35339793e-04, 1.12355128e+00, 2.58717530e-03,
       1.34975082e+01, 1.75757042e-06, 5.36062114e+02, 1.59580031e-01,
       2.08841625e-02]
E_kv = [a*1e6 for a in E_kv]


"""
E_kv = [2.27308943e+08, 5.91993268e+06, 2.32641743e+07]
tau_kv = [105.87820062,   6.64718719]
"""
"""
tau_kv = [0]
E_kv = [2,2]
"""

n_kv = len(tau_kv)





eps_dot  = 1e-2
Yc = 500
dt = 1e-4
e1,s1,d1 = bar_analytical_sol_gm(E_gmm, tau_gmm, eps_dot, n_gm, 100, Yc,dt)
e2,s2,d2 = bar_analytical_sol_gkv(E_kv, tau_kv, eps_dot, n_kv, 100, Yc,dt)
plt.plot(e1,s1[:len(e1)])
plt.plot(e2,s2[:len(e2)],'--')
plt.legend(['GKV', 'GM'])
plt.xlabel('strain')
plt.ylabel('stress (Pa)')

