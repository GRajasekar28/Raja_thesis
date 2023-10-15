# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:45:49 2022

@author: gopalsamy
"""


#### given the folder path to the results files in .npz format provides the post processing file for paraview


### create group vtk files

import os
import sys
sys.path.append('../lib')

import post_process_vtk_format as pp

import numpy as np
from mesh import simplexMesh

mesh_file_name = r'D:\VBox shared folder\visco v1.7\msh\old\beam_3_point_bending_r2.msh'  ## path to mesh file

 ## path where all the results file in .npz format are stored (each file corresponds to each time step)
folder_path = r'D:\VBox shared folder\visco v1.7\tmp\3_point_bending\280723\beam_3_point_bending_r21_speed_1mmps_gd2_dt_p05_stopping_criteria\results' 

  ## path to the force displacment data
fu_path = r'D:\VBox shared folder\visco v1.7\tmp\3_point_bending\280723\beam_3_point_bending_r21_speed_1mmps_gd2_dt_p05_stopping_criteria\force_disp.npz'

"""

mesh_file_name = r'D:\VBox shared folder\visco v1.7\msh\old\beam_offset_notch_3pb.msh'  ## path to mesh file

 ## path where all the results file in .npz format are stored (each file corresponds to each time step)
folder_path = r'D:\VBox shared folder\visco v1.7\tmp\3_point_offset_bending\beam_offset_notch_3pb1_speed_1mmps_gd2_dt_p05_split_3_new\results' 

  ## path to the force displacment data
fu_path = r'D:\VBox shared folder\visco v1.7\tmp\3_point_offset_bending\beam_offset_notch_3pb1_speed_1mmps_gd2_dt_p05_split_3_new\force_disp.npz'
"""



post_pr_file_name = 'output_post_proc'  ## final file name to open in paraview


model_type= 'LF' ##enter 'LF' for lipf-field



















##load mesh
mesh_file = simplexMesh.readGMSH(mesh_file_name)


## load data coressponding  to force displacement
ld1 = np.load(fu_path)
uimp = ld1['u']
Fy = ld1['Fy']
Fx = ld1['Fx']
cl = ld1['cl']  ## crack length


os.chdir(folder_path)


a = [i for i in os.listdir(folder_path) if os.path.isfile(i)]


##get applied dispalcements

u = []   ## list of applied displacements
for fle_name in a:
    u_appl = float(fle_name.split('_')[2])
    u.append(u_appl)

## do some sorting 
## (as the files names might not be in the order of u_appl)
srtd_indx = sorted(range(len(u)), key = lambda k: u[k])

srtd_files = [a[i] for i in srtd_indx]



## post_processes results stored in a new directory within given path
post_pr_dir = './post_proc'
os.mkdir(post_pr_dir)
os.chdir(post_pr_dir)



vtk_instance = pp.store_data_as_vtk(mesh_file)
t = [0] ## time
name = []
nodal_data= None

for it,fle_name in enumerate(srtd_files[:-1]):
    print('\rPreparing file '+ str(it+1) +'/'+ str(len(srtd_files)), end = '\r')
    ld = np.load('../'+fle_name, allow_pickle=1)
    
    #t.append(ld['time'])
    t.append(t[-1]+ld['DT'][0])
    u = ld['u']   ## nodal dispalcements
    ux = u[:,0] ; uy = u[:,1];
    d = ld['d']   ## face damage
    
    stress = ld['stress']
    sx = stress[:,0];  sy = stress[:,1]; sxy = stress[:,2];
    if model_type =='PF':
        d_nodes = ld['d_nodes']
        nodal_data = {'ux':ux, 'uy':uy,'d_vertex':d_nodes} ## For PF
    else:
        nodal_data = {'ux':ux, 'uy':uy}   ## For LF
    face_data = {'damage':d, 'sx':sx,'sy':sy,'sxy':sxy}  ## for PF face damage contains averaged values from vertices
    global_data = {'uF':np.array([uimp[it], Fy[it]])}
    
    name.append('sim'+str(it) )
    vtk_instance.save_vtk(path = name[it], point_data=nodal_data, cell_data= face_data, field_data=global_data)
    
    

vtk_instance.group_vtk(source =name, dest = post_pr_file_name, indices= t)
print('\nFind the file in '+ os.getcwd() + "\\" + post_pr_file_name + '.pvd')
    