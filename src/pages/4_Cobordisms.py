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

check_params()

st.title("Cobordism Analysis")

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
    st.session_state.sample_index = int(comp_tab.text_input("Sample start (this should be plain text containing the start index of the samples to process):", key="sample_index_input", placeholder="0"))
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



def generate_permutation(n_verts, sub_verts):
	permute_verts = [-1 for i in range(n_verts)]
	inverse_perm = [-1 for i in range(n_verts)]
	n_sub_verts = len(sub_verts)
	n_L = 0
	n_K = 0
	for i in range(n_verts):
		# print(i)
		if i in sub_verts:
			permute_verts[i] = n_L
			inverse_perm[n_L] = i
			n_L += 1
			# if i < len(sub_verts):
			# 	permute_verts[i] = i
			# 	inverse_perm[i] = i
			# else:
			# 	permute_verts[i] = i-len(sub_verts)
			# 	inverse_perm[i-len(sub_verts)] = i
				# permute_verts.append(i-n_verts+len(sub_verts))
		else:
			permute_verts[i] = n_sub_verts + n_K
			inverse_perm[n_sub_verts + n_K] = i
			n_K += 1
			# if i >= len(sub_verts):
			# 	permute_verts[i] = i
			# 	inverse_perm[i] = i
			# 	# permute_verts.append(i)
			# else:
			# 	permute_verts[i] = i+len(sub_verts)
			# 	inverse_perm[i+len(sub_verts)] = i
			# 	# permute_verts.append(i-len(sub_verts)+n_verts)
	return permute_verts, inverse_perm

def generate_permutation_AB(n_verts, sub_verts_A, sub_verts_B):
	permute_verts = []
	for i in range(n_verts):
		# print(i)
		if i in sub_verts_A:
			if i < len(sub_verts_A):
				permute_verts.append(i)
			else:
				permute_verts.append(i-n_verts+len(sub_verts_A))
		elif i in sub_verts_B:
			if i < len(sub_verts_A):
				permute.append(i+len(sub_verts_A))
			elif i < len(sub_verts_B)+len(sub_verts_A):
				permute_verts.append(i-n_verts+len(sub_verts_A)+len(sub_verts_B))
			else:
				permute_verts.append(i-len(sub_verts_A)-len(sub_verts_B)+n_verts)
		else:
			if i < len(sub_verts_A)+len(sub_verts_B):
				permute_verts.append(i-len(sub_verts_A)-len(sub_verts_B)+n_verts)
			else:
				permute_verts.append(i)
	return permute_verts


def dim_compare(x, y):
	"""! Comparison to compare list of simplicies to get them in the order for oineus
	@param x	simplex to compare
	@param y	simplex to compare

	@return -1 if x<y, 1 otherwise
	"""
	if len(x[0]) == len(y[0]):
		if x[0][0] <= y[0][0]:
			return -1
		else:
			return 1
	elif len(x[0]) < len(y[0]):
		return -1
	else:
		return 1

def translate_column(col, perm):
	for i in range(len(col)):
		col[i] = perm[col[i]]
	col.sort()
	return col

def construct_pair(simplices, n_verts, sub_verts, permute):
	simplices_not_sub = []
	simplices_sub = []
	# permute = generate_permutation(n_verts, sub_verts)
	simplices = sorted(simplices, key=cmp_to_key(dim_compare))
	sub_verts_set = set(sub_verts)
	for s in simplices:
		permuted_verts = [permute[i] for i in s[0]]
		if set(s[0]).issubset(sub_verts_set):
			simplices_sub.append([permuted_verts, s[1]])
		else:
			simplices_not_sub.append([permuted_verts, s[1]])
	set_sub = {tuple(s[0]) for s in simplices_sub}
	set_not_sub = {tuple(s[0]) for s in simplices_not_sub}
	assert set_sub.intersection(set_not_sub) == set()
	simplices_all = simplices_sub + simplices_not_sub
	simplices_sub = [[i,s[0],s[1]] for i, s in enumerate(simplices_sub)]
	simplices_all = [[i,s[0],s[1]] for i, s in enumerate(simplices_all)]
	return simplices_all, simplices_sub

def convert_vertices(vertices, permute):
	new = [permute[v] for v in vertices]
	new.sort()
	return new

def rref_d_phi(d_phi):
	columns = len(d_phi)
	rows = max([max(col) for col in d_phi])+1
	# print("rows: ", rows, " columns: ", columns)
	matrix = [[0 for j in range(columns)] for i in range(rows)]
	# print("matrix length: ", len(matrix))
	for i, col in enumerate(d_phi):
		for j in col:
			# print("i: ", i, " j: ", j)
			matrix[j][i] = 1
	matrix = sympy.Matrix(matrix)
	# print("matrix: ")
	# print(matrix)
	return matrix


def find_births(dgm_A, dgm_B, dgm_AB):
	"""! Find the points in Cok(d_phi)
	@param dgm_A	time and index diagram of A
	@param dgm_B	time and index diagram of B
	@param dgm_AB	time and index diagram of AB
	"""

	A_birth_ids = dgm_A['birth idx'].tolist()
	B_birth_ids = dgm_B['birth idx'].tolist()
	AB_birth_ids = dgm_AB['birth idx'].tolist()
	all_births = set(A_birth_ids) | set(B_birth_ids) | set(AB_birth_ids)
	# print(all_births)
	# print(A_birth_ids)
	# print(B_birth_ids)
	# print(AB_birth_ids)
	birth_ids = []
	birth_times = []
	for idx in all_births:
		# print("idx is ", idx)
		if idx in A_birth_ids and idx in B_birth_ids and idx in AB_birth_ids:
			# death_ids.append(idx)
			continue
		elif idx not in A_birth_ids and idx not in B_birth_ids and idx in AB_birth_ids:
			birth_ids.append(idx)
			# birth_times.append(dgm_AB['birth'][idx])
		elif idx in A_birth_ids and idx not in B_birth_ids and idx not in AB_birth_ids:
			# birth_ids.append(idx)
			# birth_times.append(dgm_A['birth'][idx])
			continue
		elif idx in A_birth_ids and idx in B_birth_ids and idx not in AB_birth_ids:
			# death_ids.append(idx)
			continue
		else:
			continue
	return birth_ids#, birth_times

def match_points_index(sorted_ids, cok_phi, birth_ids):
	"""! Match the points in Cok(d_phi)
	@param sorted_ids	indices of the simplices in the order they appear in the cokernel
	@param cok_phi		cokernel of d_phi
	@param birth_ids	indices of the births in the order they appear in the cokernel
	@param birth_times	times of the births in the order they appear in the cokernel
	"""	
	# print("birth ids are: ", birth_ids)
	duplicates = [item for item, count in Counter(sorted_ids).items() if count > 1]
	# print("Values that appear twice in sorted_ids:", duplicates)
	print(duplicates)

	# Get the indices for each duplicate value
	duplicate_indices = {}
	for i, value in enumerate(sorted_ids):
		if value in duplicates:
			if value in duplicate_indices:
				duplicate_indices[value].append(i)
			else:
				duplicate_indices[value] = [i]

	# print("Indices for duplicate values:", duplicate_indices)
	
	matched_births = {}
	for i in duplicates:
		# print("i is ", i,  " and ", duplicate_indices[i])
		# print(cok_phi[duplicate_indices[i][1]])
		if cok_phi[duplicate_indices[i][1]][-1] in birth_ids:
			matched_births[cok_phi[duplicate_indices[i][1]][-1]] = i
	# print("matched births are: ", matched_births)
	inf_births = []
	for i in birth_ids:
		if i not in matched_births.keys():
			inf_births.append(i)
	# print("inf births are: ", inf_births)
	index_dgm = []
	for i in matched_births.keys():
		index_dgm.append([i, matched_births[i]])
	for i in inf_births:
		index_dgm.append([i, -1])

	return index_dgm

def convert_index_diagram(index_dgm, kicr):
	"""! Convert the index diagram to the original diagram
	@param index_dgm	index diagram
	@param kicr		KICR object
	"""
	
	persistence_diagram = []
	for pt in index_dgm:
		print("looking at ", pt)
		print(kicr.fil_K.simplex(int(pt[0])).value)
		if pt[1] != -1:
			persistence_diagram.append([kicr.fil_K.simplex(int(pt[0])).value, kicr.fil_K.simplex(int(pt[1])).value, pt[0], pt[1]])
		else:
			persistence_diagram.append([kicr.fil_K.simplex(int(pt[0])).value, math.inf, pt[0], -1])
	return persistence_diagram
	
def d_phi_generator(kicr_AB, sorted_ids_A, sorted_ids_B, A_verts, B_verts, inv_AB):
	d_phi = []
	sorted_ids = []
	simplices = kicr_AB.fil_K.simplices()
	for i in range(len(simplices)):
		simplex = simplices[i]
		og_verts = sorted(convert_vertices(simplex.vertices, inv_AB))
		sorted_id_A = sorted_ids_A[tuple(og_verts)]
		sorted_id_B = sorted_ids_B[tuple(og_verts)]
		print(i, " looking at simplex ", simplex, " with sorted_A ", sorted_id_A, " and sorted_B ", sorted_id_B)
		in_A = True
		in_B = True
		new_col_A = []
		new_col_B = []
		for v in og_verts:
			if v not in A_verts:
				in_A = False
				break
			if v not in B_verts:
				in_B = False
				break
		zero_col_A = False
		zero_col_B = False
		lowest_one_in_A = False
		lowest_one_in_B = False
		if kicr_A.decomposition_im.r_data[sorted_id_A] == []:
			zero_col_A = True
			# print("zero column in R_im^A")
			cur_column = kicr_A.decomposition_im.v_data[sorted_id_A]
			# print(cur_column)
			new_col_A = []
			for s in cur_column:
				new_col_A.append(sorted_ids_AB[tuple(convert_vertices(kicr_A.fil_K.simplices()[kicr_A.new_order_to_old()[s]].vertices, inv_A))])
			new_col_A.sort()
			# print(new_col_A)
		if not zero_col_A:
			# print("non-zero column in R_im^A")
			lowest_one = kicr_A.new_order_to_old()[kicr_A.decomposition_im.r_data[sorted_id_A][-1]]
			# print("lowest one in R_im^A is ", lowest_one)
			lowest_one_vertices = convert_vertices(kicr_A.fil_K.simplices()[lowest_one].vertices, inv_A)
			lowest_one_in_A = True
			for v in lowest_one_vertices:
				if v not in A_verts:
					lowest_one_in_A = False
					break
			if lowest_one_in_A:
				# print("lowest one in R_im^A is in A")
				cur_column = kicr_A.decomposition_im.v_data[sorted_id_A]
				# print(cur_column)
				new_col_A = []
				for s in cur_column:
					new_col_A.append(sorted_ids_AB[tuple(convert_vertices(kicr_A.fil_K.simplices()[kicr_A.new_order_to_old()[s]].vertices, inv_A))])
				new_col_A.sort()
			else: 
				continue
				# print("lowest one in R_im^A is not in A")	
		# print("sorted_id_B is ", sorted_id_B)
		if kicr_B.decomposition_im.r_data[sorted_id_B] != []:
			# print("non-zero column in R_im^B")
			lowest_one = kicr_A.new_order_to_old()[kicr_B.decomposition_im.r_data[sorted_id_B][-1]]
			# print("lowest one in R_im^B is ", lowest_one)
			lowest_one_vertices = convert_vertices(kicr_B.fil_K.simplices()[lowest_one].vertices, inv_B)
			lowest_one_in_B = True
			for v in lowest_one_vertices:
				if v not in B_verts:
					lowest_one_in_B = False
					break
			if lowest_one_in_B:
				# print("lowest one in R_im^B is in B")
				cur_column = kicr_B.decomposition_im.v_data[sorted_id_B]
				# print(cur_column)
				new_col_B = []
				for s in cur_column:
					new_col_B.append(sorted_ids_AB[tuple(convert_vertices(kicr_B.fil_K.simplices()[kicr_B.new_order_to_old()[s]].vertices, inv_B))])
				new_col_B.sort()
			else:
				# print("lowest one in R_im^B is not in B")	
				continue
		# print("zero_col_A is ", zero_col_A)
		# print("zero_col_B is ", zero_col_B)
		if new_col_A != [] and new_col_B != []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			if in_A and not in_B:
				# print("in_A and not in_B")
				d_phi.append(new_col_A)
				sorted_ids.append(i)
				d_phi.append(new_col_B)
				sorted_ids.append(i)
			elif not in_A and in_B:
				# print("not in_A and in_B")
				d_phi.append(new_col_B)
				sorted_ids.append(i)
				d_phi.append(new_col_A)
				sorted_ids.append(i)
		elif new_col_A != [] and new_col_B == []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			d_phi.append(new_col_A)
			sorted_ids.append(i)
		elif new_col_A == [] and new_col_B != []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			d_phi.append(new_col_B)
			sorted_ids.append(i)
	return d_phi, sorted_ids

def cobordism_analysis():
    sample = sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.atoms, st.session_state.radii)
   
    upper = st.session_state.upper_threshold
    lower = st.session_state.lower_threshold

    simplices = weighted_alpha_diode(sample)
    
    A_verts = []
    B_verts = []
    for i in range(len(sample["z"])):
        if sample["z"].iloc[i] >= upper:
            A_verts.append(i)
        elif sample["z"].iloc[i] <= lower:
            B_verts.append(i)
    n_verts = len(sample["z"])
    perm_A, inv_A = generate_permutation(n_verts, A_verts)
    Xa, A = construct_pair(alpha, n_verts, A_verts, perm_A)
    kicr_A = oineus.compute_kernel_image_cokernel_reduction(Xa, A, KICR_params)
    perm_B, inv_B = generate_permutation(n_verts, B_verts)
    Xb, B = construct_pair(alpha, n_verts, B_verts, perm_B)
    kicr_B = oineus.compute_kernel_image_cokernel_reduction(Xb, B, KICR_params)
    perm_AB, inv_AB = generate_permutation(n_verts, A_verts+B_verts)
    Xab, AB = construct_pair(alpha, n_verts, A_verts+B_verts, perm_AB)
    kicr_AB = oineus.compute_kernel_image_cokernel_reduction(Xab, AB, KICR_params)

    sorted_ids_A = {}
    sorted_ids_B = {}
    sorted_ids_AB = {}
    for i in range(len(alpha)):
        # print("A  ",sorted(kicr_A.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_A.fil_K.simplex(i).vertices, inv_A))
        sorted_ids_A[tuple(sorted(convert_vertices(kicr_A.fil_K.simplex(i).vertices, inv_A)))] = i
        # print("B  ",sorted(kicr_B.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_B.fil_K.simplex(i).vertices, inv_B))
        sorted_ids_B[tuple(sorted(convert_vertices(kicr_B.fil_K.simplex(i).vertices, inv_B)))] = i
        # print("AB ",sorted(kicr_AB.fil_K.simplex(i).vertices), " becomes ", convert_vertices(kicr_AB.fil_K.simplex(i).vertices, inv_AB))
        sorted_ids_AB[tuple(sorted(convert_vertices(kicr_AB.fil_K.simplex(i).vertices, inv_AB)))] = i

    assert sorted_ids_A.keys() == sorted_ids_B.keys() == sorted_ids_AB.keys()

    d_phi, sorted_ids = d_phi_generator(kicr_AB, sorted_ids_A, sorted_ids_B,A_verts, B_verts, inv_AB)
    
    dgm_A_1 = kicr_A.kernel_diagrams().in_dimension(1)
    idx_dgm_A_1 = kicr_A.kernel_diagrams().index_diagram_in_dimension(1)
    dgm_B_1 = kicr_B.kernel_diagrams().in_dimension(1)
    idx_dgm_B_1 = kicr_B.kernel_diagrams().index_diagram_in_dimension(1)
    dgm_AB_1 = kicr_AB.kernel_diagrams().in_dimension(1)
    idx_dgm_AB_1 = kicr_AB.kernel_diagrams().index_diagram_in_dimension(1)

    dgm_A_pd = pandas.DataFrame(numpy.concatenate((dgm_A_1, idx_dgm_A_1), axis=1), columns=["birth", "death", "birth idx", "death idx"])
    dgm_B_pd = pandas.DataFrame(numpy.concatenate((dgm_B_1, idx_dgm_B_1), axis=1), columns=["birth", "death", "birth idx", "death idx"])
    dgm_AB_pd = pandas.DataFrame(numpy.concatenate((dgm_AB_1, idx_dgm_AB_1), axis=1), columns=["birth", "death", "birth idx", "death idx"])

    birth_ids = find_births(dgm_A_pd, dgm_B_pd, dgm_AB_pd)
    # print("ids are:", birth_ids)
    index_dgm = match_points_index(sorted_ids, dcmp_phi.r_data, birth_ids)
    cobordisms_pd = convert_index_diagram(index_dgm, kicr_AB)
    print("Cobordism persistence diagram is: ", cobordisms_pd)
    st.session_state.cobordisms_pd = cobordisms_pd

