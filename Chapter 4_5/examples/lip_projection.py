


# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 08:45:49 2022

@author: gopalsamy
"""



#######################################################################
#########  EXAMPLE : Lip-fiedl projection of a given damage field ###########################
## d-d_target in the sense of L2 norm
#######################################################################



import os
import sys
sys.path.append('../lib')
#sys.path.append('./.')

import viscoelasticity_law2D as visclaw
import viscoelasticymech2d   as mech2d
from mesh import simplexMesh, dualLipMeshTriangle
import gmshParser
import liplog
import logging

import numpy as np
import matplotlib.pylab as plt
import cProfile
import time as time_mod
import lipdamage as lipprojector




##############################################################################
#########################  DETAILs OF THE MESH ###############################
##############################################################################

## load the mesh
basemeshname = 'lip_proj_2'    ## name of the mesh
meshid = ''
meshfilename = '../msh/'+basemeshname + meshid +'.msh'


########### Dimensions of the geometry ######################

## Dimensions
thickness = 1e-3; b= thickness;  ## thickness in z-direction 
L = 1/1000 
h = 1/1000




##############################################################################
##########################  FOR WRITING RESULTS TO FILE ######################
##############################################################################

### name for outputfiles 

simulation_id = str(1)     ## give new id for performing simulation on same mesh with different setting (to save in a different folder)
expname = basemeshname+ simulation_id

### for outputfiles 
respath = '../tmp'
respath1 = os.path.join(respath,expname)
#os.mkdir(respath1)  ## Note : make sure a directory doesnt exist in the same name as expname


respath2 = os.path.join(respath1,'results')
#os.mkdir(respath2)


###  initiate logger to log the actions
log_file_name = expname+'_log'
mode = "a+"
# w: It is for write mode., r: It is for reading mode. a: It is for append mode.
# w+: Create the file if it does not exist and then open it in write mode.
# r+: Open the file in the read and write mode.
# a+: Create the file if it does not exist and open it in append mode
logger =  liplog.setLogger(respath = respath1, basefilename = log_file_name,mode= mode)






lc = .1  ## characterstic/regularization length (m)





##################################################
### Process the mesh and construct lip mesh  #####
##################################################


### start profiling
pr = cProfile.Profile()


logger.info('Starting Lip-projection')


gmsh_mesh = gmshParser.Mesh()
gmsh_mesh.read_msh(meshfilename)


mesh = simplexMesh.readGMSH(meshfilename)
nv = mesh.nvertices
nf =  mesh.ntriangles


logger.info('Mesh Loaded from file ' + meshfilename)
logger.info('Mesh size : nv: '+ str(nv) + ' nf: '+ str(nf))

### build lipMesh
lipmesh = dualLipMeshTriangle(mesh)
lipproj = lipprojector.damageProjector2D(lipmesh, verbose=True)
logger.info('LipMesh constructed from Mesh')
logger.info('LipMesh size : nv: '+ str(lipmesh.nvertices)+' nf: '+ str(lipmesh.ntriangles))







#############################################################################
################### d_target ###################################################
#############################################################################


## elements where d is enforced 1  (d_target)

videl_d1 = []
videl_d2 = []

"""
a_low = .2     ;   a_left = .2
a_uppe = 1-.2     ;   a_right = 1-.2

## damage is possible only when x, y in the  domain (of the geometry) satisfies following
##    a_low < y < a_upper   and  a_left < x < a_right    (else d = 0.)

for ele_num,ele in enumerate(mesh.triangles):
    for ver in ele:
        x_ver = mesh.xy[ver,0]
        y_ver = mesh.xy[ver,1]
        if y_ver<a_low or y_ver> a_uppe or x_ver< a_left or x_ver > a_right :
            videl_d1.append(ele_num)
            break

"""

phy_id1 = 100   ## elements physical tag where d is to be imposed 1
phys_id2 =101

for i in range(mesh.ntriangles):
    k = gmsh_mesh.getElementConnectivityAndPhysNumber(2,i)[1]
    if k == 100:
        videl_d1.append(i)
    if k == phys_id2:
        videl_d2.append(i)



d_target = np.zeros(lipmesh.nvertices)
d_target[videl_d1] = 1.

d_target[videl_d2] = .5








### Project d_target into Lip-space
d_lip = lipproj.lipProjClosestToTarget(dmin = np.zeros(mesh.ntriangles), dtarget = d_target, lc=lc, logger=logger)


"""
## Find upper and lower bounds
d_up = lipproj.lipProjFM(d_target, lc, lipmeasure = 'triangle', side ='up', logger = logger)

d_lo = lipproj.lipProjFM(d_target, lc, lipmeasure = 'triangle', side ='lo', logger = logger)
"""

"""
mesh.plotScalarField(d_target)
mesh.plotScalarField(d_lip)

mesh.plotScalarField(np.abs(d_up['d']-d_lo['d']))
"""


"""
#plt.figure(2)

fig,ax = plt.subplots()
plot_f_d(ax)
plt.xlabel('displacement (m)')
plt.ylabel('Force (N)')
fig.savefig(respath1+'/Force_displacement_.pdf', format = 'pdf') 
plt.legend()"""