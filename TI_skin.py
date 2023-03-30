# -*- coding: utf-8 -*-
"""
 example script that runs two simnibs tDCS simulations
 and calculates maximal amplitude of the TI envelope from the E-fields
 
 Created on Thu Jun 23 17:41:21 2022

@author: axthi
"""

from copy import deepcopy
import os
import numpy as np

from simnibs import sim_struct, run_simnibs, mesh_io
from simnibs.utils import TI_utils as TI

def _scalp_geo(m, fn_out, scalp_idx=1005):
        ''' write out scalp surface as geo file '''
        idx = (m.elm.tag1 == scalp_idx) & (m.elm.elm_type == 2)
        mesh_io.write_geo_triangles(m.elm[idx, :3]-1,
                                    m.nodes.node_coord, fn_out,
                                    name='scalp', mode='ba')

"""
     set up and run simulations for the two electrode pairs
"""

# specify general parameters
S = sim_struct.SESSION()
S.fnamehead = 'pp6.msh'
S.pathfem = 'TI'  # Directory for the simulation
S.fields = 'eEs'

tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (A)

electrode = tdcs.add_electrode()
electrode.channelnr = 1
electrode.centre = [153, 129, 293]  
electrode.shape = 'ellipse' 
electrode.dimensions = [5, 5] 
electrode.thickness = [0.5] 

electrode = tdcs.add_electrode()
electrode.channelnr = 2
electrode.centre = [149, 109, 290]
electrode.shape = 'ellipse'
electrode.dimensions = [5, 5]
electrode.thickness = [0.5]

# tdcs.cond[4].value = 1.89e7 # [S/m]
# tdcs.cond[4].name = 'metal_1'
# tdcs.cond[5].value = 1.89e7 # [S/m]
# tdcs.cond[5].name = 'metal_2'

# specify second electrode pair
tdcs = S.add_tdcslist(deepcopy(tdcs))
tdcs.electrode[0].centre = [244, 106, 281]
tdcs.electrode[1].centre = [240, 87, 278]

run_simnibs(S)


"""
    generate the TI field from the simulation results
"""
m1 = mesh_io.read_msh(os.path.join(S.pathfem, 'pp6_elec_skin_TDCS_1_scalar.msh'))
m2 = mesh_io.read_msh(os.path.join(S.pathfem, 'pp6_elec_skin_TDCS_2_scalar.msh'))

# remove all tetrahedra and triangles belonging to the electrodes so that
# the two meshes have same number of elements
tags_keep = np.hstack((np.arange(1,100), np.arange(1001,1100)))
m1=m1.crop_mesh(tags = tags_keep)
m2=m2.crop_mesh(tags = tags_keep)

# calculate the maximal amplitude of the TI envelope
ef1=m1.field['E']
ef2=m2.field['E']
print(ef1.value.shape)
TImax = TI.get_maxTI(ef1.value, ef2.value)
# make a new mesh for visualization of the field strengths
# and the amplitude of the TI envelope
mout = deepcopy(m1)
mout.elmdata = []
mout.add_element_field(ef1.norm(), 'magnE - pair 1')
mout.add_element_field(ef2.norm(), 'magnE - pair 2')                    
mout.add_element_field(TImax,'TImax')
# mout.add_element_field(TImax_v,'TImax_v')
mesh_io.write_msh(mout,os.path.join(S.pathfem, 'TI.msh'))
v = mout.view(
    visible_tags=[1002],
    visible_fields='TImax',    
    )
v.add_merge(os.path.join(S.pathfem, 'pp6_elec_skin_TDCS_1_el_currents.geo'))
v.add_view(ColormapNumber=10, ColormapAlpha=.8,
                   Visible=1)  # el_currents
v.add_merge(os.path.join(S.pathfem, 'pp6_elec_skin_TDCS_2_el_currents.geo'))
v.add_view(ColormapNumber=10, ColormapAlpha=.8,
                   Visible=1)  # el_currents
_scalp_geo(mout, os.path.join(S.pathfem, 'scalp.geo'))
v.add_merge(os.path.join(S.pathfem, 'scalp.geo'))
v.add_view(ColormapNumber=8, ColormapAlpha=.3,
                       Visible=0, ShowScale=0)
v.write_opt(os.path.join(S.pathfem, 'TI.msh'))
mesh_io.open_in_gmsh(os.path.join(S.pathfem, 'TI.msh'), True)
            
            

