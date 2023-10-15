# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 17:25:56 2023

@author: gopalsamy
"""

"""1D bar validation with damage (localisation triggered at center)
## make sure the trigger is larger to avoid mesh irregularities affecting the localization location"""


#######################################################################
#########  EXAMPLE 1: HOMOGENEOUS BAR  ###########################

### VALIDATION OF GKV IMPLEMENTATION WITH DAMAGE WITH (implicit) ANALYTIC SOLUTION ##########
###  can VALIDATE FOR BOTH LINEAR AND NON-LINEAR CASE (non-linear -> asymmetric case)
#######################################################################

import os
import sys
sys.path.append('../lib')
#sys.path.append('./.')

import viscoelasticity_law2D as visclaw
import viscoelasticymech2d   as mech2d
from mesh import simplexMesh, dualLipMeshTriangle
import liplog
import logging

import numpy as np
import matplotlib.pylab as plt
import cProfile
import time as time_mod
import lipdamage as lipprojector

import matplotlib


##############################################################################
#########################  DETAILs OF THE MESH ###############################
##############################################################################

############# load the mesh ##############
basemeshname = 'bar2_2d_new';   ## name of the mesh
meshid = ''
meshfilename = '../msh/'+basemeshname + meshid +'.msh'
mesh_in_mm = False       ## True when mesh details are in 'mm' else False when in 'm'

########### Dimensions of the geometry ######################
## in 'm'
thickness = 1; 
b= thickness;  ## thickness in z-direction 
L = 1 
h = .25



##############################################################################
##########################  FOR WRITING RESULTS TO FILE ######################
##############################################################################



### name for outputfiles 

simulation_id = str(2)     ## give new id for performing simulation on same mesh with different setting
expname = basemeshname+ simulation_id

### for outputfiles 
respath = '../tmp'
respath1 = os.path.join(respath,expname)
os.mkdir(respath1)  ## make sure a directory doesnt exist in the same name as expname


respath2 = os.path.join(respath1,'results')
os.mkdir(respath2)




###  initiate logger to log the actions
log_file_name = expname+'_log'
logger =  liplog.setLogger(respath = respath1, basefilename = log_file_name)





##################################################################
################# CONSTITUTIVE BEHAVIOUR #########################
##################################################################



#### Material properties  #######


## Material parameters
# E0=14000e6
# E1=5600e6
# nu = 0.0
# tau = 2.59e-2
# Yc = 30
 




E = [31770, 87398, 123414, 65830, 62457, 62661, 7305, 12500, 418, 1743, 79, 39]    ## MPa
E = [i*1e6 for i in E]     ## Pa
tau = [1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 5e2, 1e3]   ## s
Yc = 23000    ## J/m^3
ni = 11+1     ## number of variables (excluding damage variable)  (or) number of units in GKV 

"""
E = E[:3]
tau = tau[:2]
ni = 2+1
"""


E = [ 20555.40787441, 135374.72188481,  84505.60095753,  11142.58265657,
        56337.78789284,   3195.23653526, 181136.00646417,    327.13513483,
        21958.9894283 ,  34646.72374612]
tau = [1.98638904e-05, 2.35339793e-04, 1.12355128e+00, 2.58717530e-03,
       1.34975082e+01, 1.75757042e-06, 5.36062114e+02, 1.59580031e-01,
       2.08841625e-02]
E = [a*1e6 for a in E]
ni = len(E)

nu=0.0  ## Poisson's ratio

if len(tau)+1 != ni or len(E)!=ni:
    raise

nu=0  ## Poisson's ratio

Yc = 500
lc = .5

##plane strain assumption
## find the material constants lamb and mu for each spring in GKV

lamb = []; mu=[]
for i in range(ni):
    lamb.append(E[i]*nu/(1.+nu)/(1.-2*nu))
    mu.append(E[i]/2./(1+nu))

lamb = tuple(lamb)
mu = tuple(mu)

damage_calc = 1   ## set to true for performing damage calc.



unilateral = 1   ## bool for asymmetric tension/compression effects
split_choice = 1

if not unilateral:
    #### law 1: symmetric tension/compression  cons. law
    g1 = visclaw.GQuadratic()   ## degradation function g(d)
    #g2 = visclaw.G_Const_1()
    H = visclaw.HQuadratic()   ## softening function   h(d)
    law = visclaw.viscoElasticity2dPlaneStrain(ni,lamb,mu,tau,g1=g1,g2=g1,Yc = Yc,h=H,var_con=False)
else:
    ### law 2: aymmetric tension/compression cons. law    
    if split_choice is None: split_choice = 2      ## determines the type of split for free energy ( psi^+ and psi^-)
    g1 = visclaw.GQuadratic()    ## g1(d)
    g2 = visclaw.G_Const_1()     ## g2(d) = 1 
    H = visclaw.HQuadratic()    ## h(d)
    #H = visclaw.HQuadratic()
    #H = visclaw.HCubic(5)
    #H = visclaw.H_TLS_1()
    law = visclaw.viscoElasticity2dPlaneStrain_ASSIM(ni,lamb,mu,tau,split_choice=split_choice,Yc  = Yc,
                                                     g1 =g1,g2=g2,h=H)



##################################################
### Process the mesh and construct lip mesh  #####
##################################################


### start profiling
pr = cProfile.Profile()


logger.info('Starting Program 3 point bending')


mesh = simplexMesh.readGMSH(meshfilename)
if mesh_in_mm:  mesh.xy = mesh.xy/1000
nv = mesh.nvertices
nf =  mesh.ntriangles


logger.info('Mesh Loaded from file ' + meshfilename)
logger.info('Mesh size : nv: '+ str(nv) + ' nf: '+ str(nf))

### build lipMesh
lipmesh = dualLipMeshTriangle(mesh)
lipproj = lipprojector.damageProjector2D(lipmesh, verbose=True)
logger.info('LipMesh constructed from Mesh')
logger.info('LipMesh size : nv: '+ str(lipmesh.nvertices)+' nf: '+ str(lipmesh.ntriangles))

### build lipMesh
lipmesh = dualLipMeshTriangle(mesh)
lipproj = lipprojector.damageProjector2D(lipmesh, verbose=True)   ## projector for obtaining the bounds (upper  and lower bounds on damage)
logger.info('LipMesh constructed from Mesh')
logger.info('LipMesh size : nv: '+ str(lipmesh.nvertices)+' nf: '+ str(lipmesh.ntriangles))




##############################################################################
###################  INITIATE THE SOLVERS #####################################
#############################################################################

damage_solver = 'LF'     ### = 'LF' for Lip-Field solver
#damage_solver = 'PF'     ### = 'PF' for Phase-field AT2 solver


## create an instance of the mechanical class for viscoelasticity
mech = mech2d.Mechanics2D(mesh, law, lipproj,logger=logger,lc = lc)  


### DISPLACEMENT (AND INTERNAL STRAINS) SOLVER
solverdisp = None; solverdispoptions = None;
if not unilateral:
    ## load linear solver for displacements
     solverdisp = mech2d.Mechanics2D.solveDisplacementFixedDLinear
     solverdispoptions = {'linsolve':'cholmod'}
else:
    ## load non linear solver (works on Newton method with line search)
    solverdisp = mech2d.Mechanics2D.solve_u_eps_i_nonlinear
    solverdispoptions = {'linsolve':'cholmod','itmax':20, 'resmax':1.e-12,'res_energy_abs':1e-8}
    ## resmax - relative total force  on the free nodes;; 
    ## res_energy_abs- difference of incremental potential b/w 2 consective Newton iterations

### DAMAGE SOLVER
solverd = None; solverdoptions = None;
if damage_calc:
    if damage_solver == 'PF':
        solverd = mech2d.Mechanics2D.phase_field_AT2
        solverdoptions ={'linsolve':'cholmod'}
    else:
        solverd = mech2d.Mechanics2D.solveDLipBoundPatch 
        solverdoptions ={'mindeltad':1.e-3, 'fixpatchbound':False, 
                            'Patchsolver':'triangle', 'FMSolver':'triangle', 
                            'parallelpatch':False, 'snapthreshold':0.999,
                            #'kktsolveroptions': {'mode':'schur', 'linsolve':'cholmod'}
                            'kktsolveroptions': {'mode':'direct', 'linsolve':'umfpack'}
                            }


### Alternate Minimisation solver options (Staggered solver )
alternedsolver = mech.alternedSolver
alternedsolveroptions= {'abstole':1.e-8, 'reltole':1.e-6, 'deltadtol':1.e-6, 'max_iter': 30}


##trigger for loclaization

videl_d  =None

videl_d = []
a_low = 0     ;   a_left = L/2-(.5*.025)
a_uppe = h     ;   a_right = L/2 + (.5*.025)


for ele_num,ele in enumerate(mesh.triangles):
    for ver in ele:
        x_ver = mesh.xy[ver,0]
        y_ver = mesh.xy[ver,1]
        if x_ver> a_left and x_ver < a_right :
            videl_d.append(ele_num)
            break



#############################################################################
################### BCD'S ###################################################
#############################################################################

### idl0,idl1  classified edges created on mesh file for giving boundary conditions

idl0 = 8 #physical id of left support  
idl1 = 6 ## physical id of loading line



## get nodes lying on classified edges
vidl0 = mesh.getVerticesOnClassifiedEdges(idl0)
vidl1 = mesh.getVerticesOnClassifiedEdges(idl1)



if (len(vidl0)+len(vidl1)) ==0: 
    logger.warning('No Dirichlet nodes')   ##  stiffnes matrix doesn't have full rank
    raise

imposed_displacement= dict()

# Dirichlet Boundary conditions. degree of freedom are indexed by vertice Id  : ux(id) -> u(id*2), uy(id)   = u(8d*2+1)

## fix left support
for vid in vidl0:
    imposed_displacement[2*int(vid)] = 0.        ### x - direction
    imposed_displacement[2*int(vid)+1] = 0.      ### y-direction

## apply x-dispalcement in right support at time t=0
for vid in  vidl1:
    #imposed_displacement[2*int(vid)] =   0.     
    imposed_displacement[2*int(vid)+1] = 0.     ## >0 for t>0


##############################################################
## function to calculate the reaction force on the supports
##############################################################
def calc_reaction_force2(R):
    Fx_support = 0
    Fy_support = 0
    for vid2 in vidl0: 
        Fx_support += R[2*vid2]
        Fy_support += R[2*vid2+1]
    return (abs(Fx_support), abs(Fy_support))


##############################################################################
########################## analytical solution ###############################
############################################################################

## function to calculate anaytical solution of bar without damage
def bar_analytical_sol(U_dot,max_u_appl,dt=1e-1):
    u_appl= 0 
    ## list to store epsilon_i  
    ## epsilon_i[j][k] = epsilon_k for time step j
    epsilon_i = [[0]*ni]     ## internal strain  epsilon_0, epsilon_1,,...
    epsilon = [0]          ## total strain
    sigma =[0]             ## sigma (constant along the bar)
    t = [0]                ## time 
    while u_appl< max_u_appl:
        u_appl = U_dot*(t[-1]+dt)
        k1 = 0; k2=0
        for j in range(ni-1):
            k1 += tau[j] * epsilon_i[-1][j+1]/(tau[j] + dt)
            k2 += dt/((tau[j]+dt)*E[j+1])
        sigma_tdt =  ((U_dot*(t[-1]+dt)/L) -(k1))/(1/E[0] + k2)
        epsilon_tdt = U_dot*(t[-1]+dt)/L
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




## time iteration 

velocity = [1e-1, 8e-3, .1e-1, .5]           ## m/s
max_displacement_y = [.007, .001,.05,.5]      ## max_displacement applied for each speed
time_step = [1e-5,.001,.05,.1] 
tot_time = []     ## total time when the max_displacement  is attained


vi = int(0)
velocity = [velocity[vi]]
max_displacement_y = [max_displacement_y[vi]]
time_step = [time_step[vi]]



w2f = ['\nvelocity : '+str(velocity), '\nTime step : '+ str(time_step)]

## variables to store the results for force-displacement plot
## fem
uimp=[]
Fx  =[]
## anlaytic
uimp_ana = []
Fx_ana = []



## plotting options
plt.close('all')
onlineplot = False      ## for saving contour plots to pdf (takes too much space !!)
real_time_plot = False  ## for plotting  real time force-disp curves as the program runs
showmesh = False
real_time_plot_save = True  ## for saving as a png file the plot of force disp every time step in respath1


if not real_time_plot:
    matplotlib.use('Agg')  # Use the Agg backend (non-interactive)

if showmesh:
    mesh.plot()

if real_time_plot:
    plt.show()
    fig0, ax0 = plt.subplots(1,1)
    axes = plt.gca()
    axes.set_xlim(-1e-4, .007)
    max_data_plot = 3
    axes.set_ylim(-2e-2, max_data_plot)
    line, = axes.plot(0,0)
    fig0.suptitle('Stress Strain plot')
    ax0.set_xlabel(r'Strain')
    ax0.set_ylabel(r'Stress (Pa)')

## list to store time spent on disp and damage solver
timeu =[]
timed =[]

## list ot store number of iterations of alternate minimisation at each time step
alt_iter = []


## initital condition
u = np.zeros((nv,2))
d = np.zeros(nf)

## plot for damage
#figd, axd = plt.subplots() 

##time increment 

for ind, (speed, max_uy, DT) in enumerate(zip(velocity, max_displacement_y, time_step)):   
    logger.info('\n\n\n\n Simulation '+ str(ind))
    logger.info("\n Message " +str(ind) +": Solving for velocity = "+ str(speed)+ " m/s ")
    R,eps_i = mech.zeros()
    time = 0
    count = 0
    temp_disp, temp_force = [0], [0] 
    timeu.append(0); timed.append(0);
    alt_iter_tmp = []
    ## reinitialize dirichlet conditions to zero for different speed setting
    for vid in  vidl1:
            imposed_displacement[2*int(vid) ] = 0
    u_appl= 0           # initial applied displacement
    while u_appl < max_uy:
        u_appl= speed*(time+DT)
        #u_appl = np.sin((time+DT))
        for vid in  vidl1:
            imposed_displacement[2*int(vid) ] =   u_appl
        dmin = d.copy()
        dguess =d.copy()
        ##apply perturbations for localisation
        dguess[videl_d] += 5e-2
        pr.enable()
        res = alternedsolver( dmin =dmin, dguess =dguess, DT=DT, un= u.copy(), eps_i_n = eps_i.copy(), 
                    imposed_displacement = imposed_displacement, 
                    alternedsolveroptions = alternedsolveroptions,
                    solverdisp = solverdisp, solverdispoptions= solverdispoptions,
                    solverd= solverd, 
                    solverdoptions= solverdoptions,
                    damage_calc=damage_calc,adap_time_step=False)
        if not res['Converged'] :
            print('alterned didnt converge')
            #print ('alterned Solver DID NOT converge at step ui =', speed*time, ' message is ', res['info'])
            #logger.warning('alterned Solver DID NOT converge at step ui =', speed*time, ' message is ', res['info'])
            #u = res['u']
            #d = res['d']
            #mech.plots(u, d, eps1, u_appl, respath+'/'+expname+'LastBeforeFailure', showmesh = showmesh)
            #R = res['R']
            #break
        dtmp1 = res['d']
        dmind1 = np.linalg.norm(dtmp1-dmin)
        logger.info('\n End uimpminimize = '+ str(u_appl) + ' Conv :' + str(res['iter']) + ' iter, |dmin-d1|= ' + '%2.e'%dmind1+'\n\n')
        pr.disable()
        u = res['u']
        eps_i = res['eps_i']
        d = res['d']
        R = res['R']
        timeu[-1] += res['timeu']
        timed[-1] += res['timed']
        alt_iter_tmp.append(res['iter'])
        
        if onlineplot :     
            if np.max(d) > 0 :
                mech.plots(u, d,eps_i, u_appl, respath1+'/'+expname, showmesh = showmesh)
        


        ##updating local displacement curve
        temp_disp.append(u_appl/L)   ##strain_x
        temp_force.append(calc_reaction_force2(R)[0]*1/h) ##stress_x
        if real_time_plot:
            max_data_plot = max(max_data_plot, temp_force[-1])
            axes.set_ylim(-.2, 1.1*max_data_plot)
            axes.set_xlim(-1e-4, 1.1*temp_disp[-1])
            line.set_xdata(temp_disp)
            line.set_ydata(temp_force)       ## to scale the graph
            plt.draw()
            plt.pause(1e-10)   ## to update the plot
        if real_time_plot_save:
            plt.clf()
            plt.plot(temp_disp, temp_force)
            plt.xlabel('Strain')
            plt.ylabel('Stress (Pa)')
            plt.savefig(respath1+ '/force_disp.png')
        
        
        
        count+=1
        time+=DT
    
    
    tot_time.append(time)
    uimp.append(np.array(temp_disp))
    Fx.append(np.array(temp_force)*thickness)  
    alt_iter.append(alt_iter_tmp)
    
    ## analytic solution
    an_uimp, an_stress =    bar_analytical_sol(speed, max_uy,DT)
    uimp_ana.append(an_uimp)
    Fx_ana.append(an_stress)
            

def plot_f_d(ax):
    
    plot_stress_strain2(Fx,uimp,Fx_ana,uimp_ana,ax)
    #plot_force_displacement_ref_data()
    

def plot_stress_strain2(force,disp,force1,disp1,ax):
    
    for i,(f,u) in enumerate(zip(force,disp)):
        ax.plot(np.array(u),np.array(f).squeeze(), label = '2D FEM '+ str(velocity[i]*1000) + ' mm/s')
        
    for i,(f,u) in enumerate(zip(force1,disp1)):
        ax.plot(np.array(u),np.array(f).squeeze(),'*' ,label = 'Analytic 1D '+ str(velocity[i]*1000) + ' mm/s')
        
        


#plt.figure(2)
fig,ax = plt.subplots()
plot_f_d(ax)
plt.legend()
plt.xlabel('Strain ')
plt.ylabel('Stress (N/m2)')
fig.savefig(respath1+'/Force_displacement.pdf', format = 'pdf')
