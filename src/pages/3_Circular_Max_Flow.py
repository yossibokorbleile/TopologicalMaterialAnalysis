##
# @internal
# @file Circular_Max_Flow.py
# @brief Streamlit page for computing circular max flow.
# @version 0.1
# @date July 2025
# @author: Yossi Bokor Bleile
# @author: Matteo Pegoraro

import streamlit as st
import numpy
import matplotlib.pyplot as plt

from reeb_graph import Reeb_Graph
from reeb_aux import *
from utils import *

st.header("Circular Max Flow")

st.markdown("""
This page allows you to compute the circular max flow of a configuration, as implemented in the [Circular Max Flow](https://github.com/pego91C/ircular_Max_Flow) repository.
""")

st.markdown("""
It creates a static approximation of the backbone, by taking approximate locations of atoms in the backbone, and for each atom type a radius to thicken these atoms by. 
""")

st.markdown("""
Currently, it takes as input a single '.csv' file, that contains the coordinates of the atoms in the configuration in columns 'x', 'y', 'z', and then the radius to use for each atom in the column 'r'.
""")

st.text_input("Input file", key="input_file", placeholder="path/to/input.csv")

st.text_input("Backbone atoms", key="backbone_atoms", placeholder="Backbone atoms")

st.text_input("Flow atom", key="flow_atom", placeholder="Flow atoms")

st.text_input("Grid size in each dimension", key="grid_size", placeholder="Number of grid points in each dimension")

st.text_input("Fat radius", key="fat", placeholder="Fat radius")

st.check_box("Periodic conditions", key="periodic", placeholder="Periodic")

st.text_input("Reeb stride", key="reeb_stride", placeholder="Reeb stride")

st.text_input("Stride", key="stride", placeholder="Stride")


def compute_circular_max_flow():
	relax = [int(st.session_state.grid_size)//2,-(int(st.session_state.grid_size)//2+1)]

	reeb = Reeb_Graph(X, M = M, m = m, radii = radii,
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

if "max_flow_computed" in st.session_state and st.session_state.max_flow_computed:
	st.markdown(f"Computed max flow: {st.session_state.max_flow}.")