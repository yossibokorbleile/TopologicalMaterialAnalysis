##
# @internal
# @file Circular_Max_Flow.py
# @brief Streamlit page for computing circular max flow.
# @version 0.1
# @date September 2025
# @author: Yossi Bokor Bleile
# @author: Matteo Pegoraro

import streamlit as st
import numpy
import matplotlib.pyplot as plt

from reeb_graph import Reeb_Graph
from reeb_aux import *
from utils import *
from diffusion_utils import *
from diffusion import *
from toma_functions import *

st.header("Circular Max Flow")

st.markdown("""
This page allows you to compute the circular max flow of a configuration, as implemented in the [Circular Max Flow](https://github.com/pego91C/ircular_Max_Flow) repository.
""")

st.markdown("""
It creates a static approximation of the backbone, by taking approximate locations of atoms in the backbone, and for each atom type a radius to thicken these atoms by. 
""")

# st.markdown("""
# Currently, it takes as input a single '.csv' file, that contains the coordinates of the atoms in the configuration in columns 'x', 'y', 'z', and then the radius to use for each atom in the column 'r'.
# """)

st.text_input("Input file", key="input_file", placeholder="path/to/input.csv")
st.text_input("File format (this should be plain text containing the format of the initial structure file):", key="file_format",placeholder="Auto") #specify format of the initial strutcure file

st.text_input("Backbone atoms", key="backbone_atoms_input", placeholder="Backbone atoms")

st.text_input("Flow atoms", key="flow_atoms_input", placeholder="Flow atoms")

st.text_input("Grid size in each dimension", key="grid_size_input", placeholder="Number of grid points in each dimension")

st.text_input("Fat radius", key="fat_input", placeholder="Fat radius")

st.checkbox("Periodic conditions", key="periodic")

st.text_input("Reeb stride", key="reeb_stride_input", placeholder="Reeb stride")

st.text_input("Stride", key="stride_input", placeholder="Stride")




def compute_circular_max_flow():
	st.session_state.atoms = []
	st.session_state.backbone_atoms = []
	st.session_state.flow_atoms = []

	for a in st.session_state.backbone_atoms_input.split(","):
		st.session_state.backbone_atoms.append(str(a).strip())
		st.session_state.atoms.append(str(a).strip())

	for a in st.session_state.flow_atoms_input.split(","):
		st.session_state.flow_atoms.append(str(a).strip())
		st.session_state.atoms.append(str(a).strip())
	# print("got atom types")
	st.session_state.n_backbone_types = len(st.session_state.backbone_atoms)
	st.session_state.n_flow_types = len(st.session_state.flow_atoms)
	st.session_state.grid_size = int(st.session_state.grid_size_input)
	st.session_state.fat = float(st.session_state.fat_input)
	st.session_state.reeb_stride = int(st.session_state.reeb_stride_input)
	st.session_state.stride = int(st.session_state.stride_input)

	st.session_state.atom_types, st.session_state.atom_coords, st.session_state.cell = sample_all_diffusion(st.session_state.input_file, st.session_state.file_format, st.session_state.reeb_stride)
	# print("loaded atom coords:", st.session_state.atom_coords)
	# st.session_state.cell = st.session_state.atom_coords[0].get_cell()[:]
	st.session_state.m = numpy.array([0,0,0])
	st.session_state.M = numpy.array([max(st.session_state.cell[:,0]), max(st.session_state.cell[:,1]), max(st.session_state.cell[:,2])])
	st.session_state.backbone_coords = []
	st.session_state.flow_coords = [] 
	st.session_state.backbone_idxs = []
	st.session_state.flow_idxs = []
	st.session_state.atom_types = numpy.array(st.session_state.atom_types)
	# print("atom types:")
	# print(st.session_state.atom_types)
	# print(type(st.session_state.atom_types))
	st.session_state.atom_idxs = []
	for i in range(len(st.session_state.atoms)):
		st.session_state.atom_idxs.append(numpy.where(numpy.isin(st.session_state.atom_types,[st.session_state.atoms[i]]))[0].tolist())
		# print(st.session_state.atom_idxs)
		if i < st.session_state.n_backbone_types:
			st.session_state.backbone_idxs.extend(st.session_state.atom_idxs[i])
			# print("backbone idxs:")
			# print(st.session_state.backbone_idxs)
		else:
			st.session_state.flow_idxs.extend(st.session_state.atom_idxs[i])
	# print("got backbone and flow idxs:", st.session_state.backbone_idxs, st.session_state.flow_idxs)
	st.session_state.backbone_atom_types = st.session_state.atom_types[st.session_state.backbone_idxs]
	st.session_state.flow_atom_types = st.session_state.atom_types[st.session_state.flow_idxs]
	# print("got backbone and flow coords")
	st.session_state.backbone_coords = st.session_state.atom_coords[:,st.session_state.backbone_idxs,:]
	st.session_state.flow_coords = st.session_state.atom_coords[:,st.session_state.flow_idxs,:]
	st.session_state.backbone_mean = point_cloud_frechet_mean_numba(st.session_state.backbone_coords, st.session_state.M, st.session_state.m, subsample=min([20,len(st.session_state.atom_coords)]), tol=0.001, maxiter = 20)   
	print("computed backbone means")
	st.session_state.M_flow = preprocess_PATHS(st.session_state.backbone_mean, st.session_state.backbone_coords, st.session_state.flow_coords, st.session_state.M, st.session_state.m) 
	# print("preprocessed flow coords")
	# st.session_state.D = flow_to_fmean_dist(st.session_state.M_flow, st.session_state.backbone_mean, st.session_state.M, st.session_state.m) 
	# print("computed flow to backbone distances")
	# print(st.session_state.D)

	st.session_state.radii = []
	for a in st.session_state.backbone_atoms:
		a_mean = st.session_state.backbone_mean[st.session_state.backbone_atom_types == a]
		st.session_state.radii.append(numpy.min(flow_to_fmean_dist(st.session_state.M_flow, a_mean, st.session_state.M, st.session_state.m) ))

	relax = [int(st.session_state.grid_size)//2,-(int(st.session_state.grid_size)//2+1)]
	covering = np.array([-1,1])
	st.session_state.backbone_radii = []
	for a in st.session_state.backbone_atom_types:
		st.session_state.backbone_radii.append(float(st.session_state.radii[st.session_state.backbone_atoms.index(a)]))
	
	
	reeb = Reeb_Graph(inputfile=None, backbone = st.session_state.backbone_mean, flow = st.session_state.flow_coords, radii = st.session_state.backbone_radii, M = st.session_state.M, m = st.session_state.m,
		grid_size = int(st.session_state.grid_size), 
		periodic = st.session_state.periodic,
		fat_radius = float(st.session_state.fat),
		covering = covering,
		reeb_stride = int(st.session_state.reeb_stride),
		transform_points = None,
		swap_res_grid_and_balls = True, #if this is set to false, you use the complementary of the thickening
		relax_z_axis = relax,
		verbose = False, save_RAM = True, stride=int(st.session_state.stride), MP=False)
	reeb.make_reeb_graph(plot=False)

	try:
		flow = circular_max_flow(reeb)*(reeb.unit_2d)
	except:
		flow=0

	st.session_state.max_flow = flow
	st.session_state.max_flow_computed = True

def test():
	st.session_state.input_file = "../examples/test_cmf.statsis2"
	st.session_state.backbone_atoms_input = "S,Si"
	st.session_state.flow_atoms_input = "Li"
	st.session_state.grid_size_input = "10"
	st.session_state.fat_input = "1"
	st.session_state.reeb_stride_input = "1"
	st.session_state.stride_input = "1"
	compute_circular_max_flow()

st.button("Test", on_click=test)

st.button("Compute circular max flow", on_click=compute_circular_max_flow)

if "max_flow_computed" in st.session_state and st.session_state.max_flow_computed:
	st.markdown(f"Computed max flow: {st.session_state.max_flow}.")