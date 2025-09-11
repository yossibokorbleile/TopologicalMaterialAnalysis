##
# @internal
# @file Cobordisms.py
# @brief Streamlit page for analysing cobordisms.
# @version 0.1
# @date July 2025
# @author: Yossi Bokor Bleile

import streamlit as st
# import oineus
import numpy as np
import pandas
import os
import configparser
from ase import io, Atoms
import diode
import math
from colour import Color
from scipy.interpolate import interpn
from functools import cmp_to_key
from scipy.interpolate import interpn
import plotly.express as px
import plotly.graph_objects as go

from toma_functions import *
from cobordisms import *

check_params()

#Define the cobordism analysis function
def cobordism_analysis():
	"""! Compute the cobordism
	@brief Compute the cobordism
	
	This function computes the cobordism of the two structures.
	"""
	if st.session_state.test:
		print("test is true")
		st.session_state.file_path = "/home/ybleile/Seafile/tunnels/md_wrapped.statsis2"
		st.session_state.atoms=["S", "Si"]
		st.session_state.radii = [1.0, 1.0]
		st.session_state.sample_index = 0
		st.session_state.repeat_x = 1
		st.session_state.repeat_y = 1
		st.session_state.repeat_z = 1
		st.session_state.upper_threshold = 30.0
		st.session_state.lower_threshold = 10.0
		st.session_state.n_threads = 30
	else:
		print("test is false")
	sample = sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.atoms, st.session_state.radii)
	st.session_state.kicr_params.kernel = True
	st.session_state.kicr_params.cokernel = True
	st.session_state.kicr_params.image = True
	st.session_state.kicr_params.codomain = True
   
	upper = float(st.session_state.upper_threshold)
	lower = float(st.session_state.lower_threshold)

	simplices = weighted_alpha_diode(sample)
	
	A_verts = []
	B_verts = []
	for i in range(len(sample["z"])):
		if sample.loc[i, "z"] >= upper:
			A_verts.append(i)
		elif sample.loc[i, "z"] <= lower:
			B_verts.append(i)
	n_verts = len(sample["z"])
	perm_A, inv_A = generate_permutation(n_verts, A_verts)
	Xa, A = construct_pair(simplices, n_verts, A_verts, perm_A)
	kicr_A = oineus.compute_kernel_image_cokernel_reduction(Xa, A, st.session_state.kicr_params)
	perm_B, inv_B = generate_permutation(n_verts, B_verts)
	Xb, B = construct_pair(simplices, n_verts, B_verts, perm_B)
	kicr_B = oineus.compute_kernel_image_cokernel_reduction(Xb, B, st.session_state.kicr_params)
	perm_AB, inv_AB = generate_permutation(n_verts, A_verts+B_verts)
	Xab, AB = construct_pair(simplices, n_verts, A_verts+B_verts, perm_AB)
	kicr_AB = oineus.compute_kernel_image_cokernel_reduction(Xab, AB, st.session_state.kicr_params)

	sorted_ids_A = {}
	sorted_ids_B = {}
	sorted_ids_AB = {}
	for i in range(len(simplices)):
		# print("A  ",sorted(kicr_A.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_A.fil_K.simplex(i).vertices, inv_A))
		sorted_ids_A[tuple(sorted(convert_vertices(kicr_A.fil_K.simplex(i).vertices, inv_A)))] = i
		# print("B  ",sorted(kicr_B.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_B.fil_K.simplex(i).vertices, inv_B))
		sorted_ids_B[tuple(sorted(convert_vertices(kicr_B.fil_K.simplex(i).vertices, inv_B)))] = i
		# print("AB ",sorted(kicr_AB.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_AB.fil_K.simplex(i).vertices, inv_AB))
		sorted_ids_AB[tuple(sorted(convert_vertices(kicr_AB.fil_K.simplex(i).vertices, inv_AB)))] = i

	assert sorted_ids_A.keys() == sorted_ids_B.keys() == sorted_ids_AB.keys()

	d_phi, sorted_ids = d_phi_generator(kicr_AB, kicr_A, kicr_B, sorted_ids_AB, sorted_ids_A, sorted_ids_B, A_verts, B_verts, perm_AB, perm_A, perm_B, inv_AB, inv_A, inv_B)
	n_rows = max([max(c) for c in d_phi])+1
	dcmp_phi = oineus.Decomposition(d_phi, n_rows, False)
	st.session_state.params.verbose = True
	st.session_state.params.compute_v = True
	dcmp_phi.reduce(st.session_state.params)
	dgm_A = numpy.append(kicr_A.kernel_diagrams().in_dimension(0),kicr_A.kernel_diagrams().in_dimension(1),axis=0)
	numpy.append(dgm_A,kicr_A.kernel_diagrams().in_dimension(2),axis=0)

	idx_dgm_A = numpy.append(kicr_A.kernel_diagrams().index_diagram_in_dimension(0),kicr_A.kernel_diagrams().index_diagram_in_dimension(1),axis=0)
	numpy.append(idx_dgm_A,kicr_A.kernel_diagrams().index_diagram_in_dimension(2),axis=0)

	dgm_B = numpy.append(kicr_B.kernel_diagrams().in_dimension(0),kicr_B.kernel_diagrams().in_dimension(1),axis=0)
	numpy.append(dgm_B,kicr_B.kernel_diagrams().in_dimension(2),axis=0)

	idx_dgm_B = numpy.append(kicr_B.kernel_diagrams().index_diagram_in_dimension(0),kicr_B.kernel_diagrams().index_diagram_in_dimension(1),axis=0)
	numpy.append(idx_dgm_B,kicr_B.kernel_diagrams().index_diagram_in_dimension(2),axis=0)

	dgm_AB = numpy.append(kicr_AB.kernel_diagrams().in_dimension(0),kicr_AB.kernel_diagrams().in_dimension(1),axis=0)
	numpy.append(dgm_AB,kicr_AB.kernel_diagrams().in_dimension(2),axis=0)

	idx_dgm_AB = numpy.append(kicr_AB.kernel_diagrams().index_diagram_in_dimension(0),kicr_AB.kernel_diagrams().index_diagram_in_dimension(1),axis=0)
	numpy.append(idx_dgm_A,kicr_A.kernel_diagrams().index_diagram_in_dimension(2),axis=0)




	dgm_A_pd = pandas.DataFrame(numpy.concatenate((dgm_A, idx_dgm_A), axis=1), columns=["birth", "death", "birth idx", "death idx"])
	dgm_A_pd["birth idx"] = dgm_A_pd["birth idx"].astype(int)
	dgm_A_pd["death idx"] = dgm_A_pd["death idx"].astype(int)
	dgm_B_pd = pandas.DataFrame(numpy.concatenate((dgm_B, idx_dgm_B), axis=1), columns=["birth", "death", "birth idx", "death idx"])
	dgm_B_pd["birth idx"] = dgm_B_pd["birth idx"].astype(int)
	dgm_B_pd["death idx"] = dgm_B_pd["death idx"].astype(int)
	dgm_AB_pd = pandas.DataFrame(numpy.concatenate((dgm_AB, idx_dgm_AB), axis=1), columns=["birth", "death", "birth idx", "death idx"])
	dgm_AB_pd["birth idx"] = dgm_AB_pd["birth idx"].astype(int)
	dgm_AB_pd["death idx"] = dgm_AB_pd["death idx"].astype(int)
	print("kicr_A is ", kicr_A)
	print("kicr_B is ", kicr_B)
	print("kicr_AB is ", kicr_AB)

	birth_ids = find_births(dgm_A_pd, dgm_B_pd, dgm_AB_pd)
	# print("ids are:", birth_ids)
	index_dgm = match_points_index(sorted_ids, dcmp_phi.r_data, birth_ids)
	cobordisms_pd = convert_index_diagram(index_dgm, kicr_AB)
	print("Cobordism persistence diagram is: ", cobordisms_pd)
	st.session_state.cobordisms_pd = cobordisms_pd



st.title("Cobordism Analysis")
st.checkbox("test", key="test")
st.checkbox("Specify directory and file prefix for saving plots", key="custom_save")
if st.session_state.custom_save:
	st.session_state.save_directory = st.text_input("Save directory (this should be plain text containing the path to the directory to save the files):", key="save_directory_input", placeholder="path/to/save/directory")
	st.session_state.save_file_name = st.text_input("File name to use when saving the files (this should be plain text containing the prefix for the files):", key="save_file_name_input", placeholder="save_file_name")
comp_tab, plot_tab, vis_tab = st.tabs(["Computation", "Plots", "Visuatlisation"]) #create tabs for the various parts
st.session_state.mode = "multi"


file_path = comp_tab.text_input("Initial structure file (this should be plain text containing the path to the initial structure file):",key="file_path", placeholder="path/to/file.xyz") #specify initial structure file 
file_format = comp_tab.text_input("File format (this should be plain text containing the format of the initial structure file):", key="file_format",placeholder="Auto") #specify format of the initial strutcure file
if "processed_file" in st.session_state:
	if st.session_state.file_path != st.session_state["processed_file"]:
		st.session_state.processed = False

comp_tab.header("Configuration settings")
manual_config = comp_tab.checkbox("Manually specify configuration", key="manual_config")#manually set configuration
if not manual_config:
	st.session_state.config_file = comp_tab.text_input("Configuration file (this should be plain text containing the path to the configuration ini file):", key="configuration_file_input", placeholder="path/to/config.ini")
	st.session_state.config_name = comp_tab.text_input("Section name (this should be plain text containing the name of the section in the file):", key="configuration_name", placeholder="config")
else:
	try:
		st.session_state.atoms = [str(a).strip() for a in comp_tab.text_input("Atoms (this should be plain text of the atomic symbols in the structure separated by commas):", key="atoms_input", placeholder="H,C,N,O").split(",")]
	except:
		st.session_state.atoms = ["H","C","N","O"]
	try:
		st.session_state.radii = [float(r) for r in comp_tab.text_input("Radii (this should be plain text of the atomic radii to use for the computation, separated by commas, in the same order as the atoms):", key="radii_input", placeholder="0.3,0.7,1.0,1.2").split(",")]
	except:
		st.session_state.radii = [0.3,0.7,1.0,1.2]


comp_tab.markdown("Computation settings")
# manual_compute = comp_tab.checkbox("Manually specify settings for the computations (i.e number of threds, and if you want to compute kernel/image/cokernel)", key="maual_comp_config")
# if not manual_compute:
# 	not_same_config_file = comp_tab.checkbox("The computation settings are in a different configuration file.", key="not_same_config_file")
# 	if not_same_config_file:
# 		st.session_state.comp_file = comp_tab.text_input("Configuration file (this should be plain text containing the path to the configuration ini file):", key="comp_config_file", placeholder="path/to/config.ini")
# 	else:
# 		st.session_state.comp_file = st.session_state.config_file
# 	st.session_state.comp_name = comp_tab.text_input("Section name (this should be plain text containing the name of the section in the file):", key="comp_config_name", placeholder="comp")
# else:
try:
	st.session_state.sample_index = int(comp_tab.text_input("Sample index (this should be plain text containing the index of the samples to process):", key="sample_index_input", placeholder="0"))
except:
	st.session_state.sample_index = 0
try:
	st.session_state.repeat_x = int(comp_tab.text_input("Repitition in x-axis (this should be plain text containing the number of times to repeat the structure in the x-axis):", key="repeat_x_input", placeholder="1"))
except:
	st.session_state.repeat_x = 1
try:
	st.session_state.repeat_y = int(comp_tab.text_input("Repitition in y-axis (this should be plain text containing the number of times to repeat the structure in the y-axis):", key="repeat_y_input", placeholder="1"))
except:
	st.session_state.repeat_y = 1
try:
	st.session_state.repeat_z = int(comp_tab.text_input("Repitition in z-axis (this should be plain text containing the number of times to repeat the structure in the z-axis):", key="repeat_z_input", placeholder="1"))
except:
	st.session_state.repeat_z = 1
st.session_state.upper_threshold = comp_tab.text_input("Select upper threshold for the cobordism (this should be plain text containing the upper threshold for the cobordism):", key="upper_threshold_input", placeholder="")
st.session_state.lower_threshold = comp_tab.text_input("Select lower threshold for the cobordism (this should be plain text containing the lower threshold for the cobordism):", key="lower_threshold_input", placeholder="")
st.session_state.n_threads = comp_tab.text_input("Select number of threads to use (this should be plain text containing the number of threads to use):", key="n_threads_input", placeholder="")
if st.session_state.n_threads == "":
	st.session_state.params.n_threads = 4
else:
	st.session_state.params.n_threads = int(st.session_state.n_threads)
comp_tab.write("Number of threads is "+str(st.session_state.params.n_threads))

st.button("Compute cobordism", on_click=cobordism_analysis)



