##
# @internal
# @file process.py
# @brief functions to process the data using oineus and diode.
# @version 0.5
# @date December 2024

import streamlit as st
import oineus
import numpy as np
import pandas as pd
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

from toma_io import *


def read_sample(structure_file : str, configuration : str):
	"""! import a specified sample range from a configuration file
	
	@param file_path    path to the file to use for configurations
	@param sample_range    name of the structure to use
	
	@result sample_start, sample_end, sample_step  first sample, last sample, step between each one
	"""
	config = configparser.ConfigParser() #config parser to read the file
	structures = config.read(configuration_file) #load the file which contains the information
	sample_start = int(config.get(configuration, "START")) #read sample start
	sample_end = int(config.get(configuration, "END")) #read sample end
	sample_step = int(config.get(configuration, "STEP")) #read time step
	print("Have read the following settings for sampling:")
	print("sample_start:", sample_start) #print sample start
	print("sample_end: ", sample_end) #print sample end
	print("sample_step: ", sample_step) #print sample step
	return sample_start, sample_end, sample_step

def read_oineus_settings(structure_file : str, setting_name : str):
	"""! import settings for kernel/image/cokernel 
	@param file_path	path to the ini file containing the settings
	@param settings_name	name of the settings to use
	
	@result params 	oineus.ReductionParams with the settings to use
	"""
	config = configparser.ConfigParser() #config parser to read the file
	settings = config.read(structure_file) #load file containing the information
	params = oineus.ReductionParams()
	params.n_threads = int(config.get(setting_name, "N_THREADS"))
	if (str(config.get(setting_name, "KERNEL")) == "TRUE") or (str(config.get(setting_name, "KERNEL")) == "True") or (str(config.get(setting_name, "KERNEL")) == "T") or (str(config.get(setting_name, "KERNEL")) == "true") or (str(config.get(setting_name, "KERNEL")) == "t"): #if kernel should be calculated
		params.kernel = True
	if (str(config.get(setting_name, "IMAGE")) == "TRUE") or (str(config.get(setting_name, "IMAGE")) == "True") or (str(config.get(setting_name, "IMAGE")) == "T") or (str(config.get(setting_name, "IMAGE")) == "true") or (str(config.get(setting_name, "IMAGE")) == "t"): #if image should be calculated
		params.image = True
	if (str(config.get(setting_name, "COKERNEL")) == "TRUE") or (str(config.get(setting_name, "COKERNEL")) == "True") or (str(config.get(setting_name, "COKERNEL")) == "T") or (str(config.get(setting_name, "COKERNEL")) == "true") or (str(config.get(setting_name, "COKERNEL")) == "t"): #if image should be calculated 
		params.cokernel = True
	if (str(config.get(setting_name, "VERBOSE")) == "TRUE") or (str(config.get(setting_name, "VERBOSE")) == "True") or (str(config.get(setting_name, "VERBOSE")) == "T") or (str(config.get(setting_name, "VERBOSE")) == "true") or (str(config.get(setting_name, "VERBOSE")) == "t"): #verbose settings
		params.verbose = True
	return params


# def load_atom_file(file_path : str, format : str, index):
# 	"""! load the file containing the initial configureation
# 	@param file_path file to load
# 	@return atoms     what we obtain from ase.io.read
# 	"""
# 	atoms = io.read(file_path, format=format, index = index)#, format = format) #just read the atoms file
# 	return atoms

def sample_at(file_path : str, format : str, sample_index, repeat_x : int, repeat_y : int, repeat_z : int, atom_list, radius_list):
	"""! Sample a structure at a particular time, with cell repetitions.
	@param atoms		initial configuration
	@param sample_index 	time to sample at
	@param repeat_x 	repetition in x dir
	@param repeat_y 	repetition in y dir
	@param repeat_z 	repetition in z dir
	@parm atom_list 	list of atoms in the config
	@param radius_list	list of radii to use for the atoms

	@return points		data frame of the points including the repetitions with columns 'Atom', 'x', 'y', 'z', 'w'
	"""
	print("Repeating as follows: ",repeat_x, repeat_y, repeat_z)
	if format == "Auto":
		sample = io.read(file_path, index=sample_index).repeat((repeat_x, repeat_y, repeat_z))
	else:
		sample = io.read(file_path, format=format, index=sample_index).repeat((repeat_x, repeat_y, repeat_z)) #get the sample of the atomic structure at sample_index and repeat it as apprpriate
	#sample = atoms[sample_index].repeat((repeat_x, repeat_y, repeat_z)) #get the sample of the atomic structure at sample_index and repeat it as apprpriate
	coord = sample.get_positions() #get the coordinates of the atoms once repeated
	cell = sample.get_cell() #get the cell size
	dfpoints = pd.DataFrame(np.column_stack([sample.get_chemical_symbols(), coord]), columns=["Atom", "x", "y", "z"]) #combine the atomic symbols with their location into a pd.DataFrame
	atoms_found = dfpoints["Atom"].unique() #get a list of all the atoms found in the structure
	remove = [] #any atoms in the sample which were not in the atoms list are removed 
	print("Found atoms of the following types: ", atoms_found, " and atom list is ", atom_list) #print all of the ones we found
	for a in atoms_found:
		if a not in atom_list:
			remove.append(a)
	if len(remove) != 0:
		print("We are going to remove atoms of following types ", remove, " as they were not specified in the configuration.")
		dfpoints = dfpoints[dfpoints["Atom"].isin(atom_list)] #remove the corresponding atoms
	conditions = [(dfpoints["Atom"]==atom_list[i]) for i in range(len(atom_list))] #set conditions to select the radius  by atom type
	choice_weights = [radius_list[i]**2 for i in range(len(radius_list))] #the weights are the radius squared 
	print("Conditions are size ", len(conditions), " and choice weights are ", choice_weights)
	dfpoints["w"] = np.select(conditions, choice_weights) #set the weights in the dataframe
	dfpoints["x"] = pd.to_numeric(dfpoints["x"]) #ensure that everything is numeric and not string
	dfpoints["y"] = pd.to_numeric(dfpoints["y"]) #ensure that everything is numeric and not string
	dfpoints["z"] = pd.to_numeric(dfpoints["z"]) #ensure that everything is numeric and not string
	return dfpoints #return the dataframe of points


def weighted_alpha_diode(points):
	"""! Use diode to fill the weighted alpha shapes
	@param points 	pd.DataFrame with columns 'x', 'y', 'z' for coordinates, and column 'w' with the weights.

	@return weighted alpha shape from diode.
	"""
	return diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())

def convert_simps_to_oineus(simplices : list): 
	"""! Diode is set to create simplices for dionysus, so we need to convert them to the correct type for oineus.
	@param simplices 	a list of simplices from diode
	
	@return oin_simps	a list of oineus simplices
	"""
	oin_simps = [oineus.Simplex_double(s[0], s[1]) for s in simplices]
	return oin_simps

def oineus_compare(x, y):
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

def sub_complex(points : pd.DataFrame, z_upper : float, z_lower : float):
	"""! Given the points, and the upper and lower thresholds in the 'z'-component. 

	@param points		pd.DataFrame containing of the points.
	@param z_upper		float giving the upper threshold, any point above this is in the subcomplex
	@param z_lower		float giving the lower threshold, any point below this is in the subcomplex

	@return sub_comp	list containing the indices of the points on which we build the subcomplex
	"""
	print("The upper threshold is {} and the lower threshold is {}".format(z_upper, z_lower))
	sub_comp = []
	for i in points.index.values:
		if (points["z"][i] >= z_upper) or (points["z"][i]    <= z_lower):
			sub_comp.append(True)
		else:
			sub_comp.append(False)
	return sub_comp     

def oineus_filtration(points : pd.DataFrame, params : oineus.ReductionParams):
	"""! Given a set of points, compute the oineus.filtration of the alpha complex
 	
	@param points		pd.DataFrame containing points and their weights
	@param params		oineus.ReductionParams which contains the settings for oineus
 
	@return K			oineus.filtration
  	"""
	simplices = diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())
	for i in range(len(simplices)):
		simplices[i] = [sorted(simplices[i][0]), simplices[i][1]]
	simplices = sorted(simplices, key=cmp_to_key(oineus_compare))
	K = [[i,s[0],s[1]] for i, s in enumerate(simplices)]
	K = oineus.list_to_filtration(K)
	return K
		
def oineus_pair(points : pd.DataFrame, sub : list):
	"""! Given a set of points, and the points that are in the subset L, construct the complexes and map between them. The subcomplex L will consists of all simplices whose vertex are in the subset.

	@param points		pd.DataFrame containing the points and their weights
	@param sub			a list containing the indices of the points on which we construct the subcomplex

	@return K			list of simplices for the entire complex, as needed by oineus
	@return L			list of simplices for the subcomplex, as needed by oineus
	@return L_to_K		list which tells you how to map the simplices in L to the simplices in K
	"""
	points["sub"]=sub
	points.sort_values(by="sub", ascending=False)
	simplices = diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())
	for i in range(len(simplices)):
		simplices[i] = [sorted(simplices[i][0]), simplices[i][1]]
	simplices = sorted(simplices, key=cmp_to_key(oineus_compare))
	L = []
	not_L = []
	for i,s in enumerate(simplices):
		if len(s[0])==1:
			if sub[s[0][0]]==True:
				L.append([s[0], s[1]])
			else:
				not_L.append([s[0],s[1]])
		else:
			sub_complex = True
			for v in s[0]:
				if sub[v] == False:
					sub_complex=False
					break
			if sub_complex == True:
				L.append([s[0], s[1]])
			else:
				not_L.append([s[0],s[1]])
	K = []
	for s in L:
		K.append(s)
	for s in not_L:
		K.append(s)
	L = [[i,s[0],s[1]] for i, s in enumerate(L)]
	K = [[i,s[0],s[1]] for i, s in enumerate(K)]
	return K, L#, L_to_K
 
def oineus_process(points : pd.DataFrame, params : oineus.ReductionParams):
	"""! Given some points with weights, and the number of threads to use, obtain the persistent homology of the weighted alpha complex of these points, using oineus.
	
	@param points		pd.DataFrame of the points, with colums 'x','y','z','w'
	@param params		oineus.ReudctionParams
	
	@return dgm_1 		pd.DataFrame of the indexed dimension 1 diagram
	@return dgm_2		pd.DataFrame of the indexed dimension 2 diagram
	"""
	filt = oineus_filtration(points, params) #get the filtration for oineus
	dcmp =  oineus.Decomposition(filt, False) #initialise the decomposition without cohomology
	dcmp.reduce(params) #reduce the matrix
	dgms = dcmp.diagram(filt) #initialise the diagrams
	dgm_1 = pd.DataFrame(np.hstack([dgms.in_dimension(1), dgms.index_diagram_in_dimension(1)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	dgm_1["birth simplex"]=dgm_1["birth simplex"].astype(int) #convert indices to int
	dgm_1["death simplex"]=dgm_1["death simplex"].astype(int) #convert indices to int
	dgm_2 = pd.DataFrame(np.hstack([dgms.in_dimension(2), dgms.index_diagram_in_dimension(2)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 2 diagram
	dgm_2["birth simplex"]=dgm_2["birth simplex"].astype(int) #convert indices to int
	dgm_2["death simplex"]=dgm_2["death simplex"].astype(int) #convert indices to int
	return dcmp, filt, dgm_1, dgm_2
	
def oineus_kernel_image_cokernel(points : pd.DataFrame, params : oineus.ReductionParams, upper_threshold : float, lower_threshold : float):
	"""! Given points, and parameters for oineus, calculate the kernel/image/cokernel persistence as desired.

	@param points			pd.DataFrame of the points, with columns 'x','y','z','w' corresponding to the coordinates and weights respectively
	@param kernel 			boolean parameter to set if kernel persistence is calculated
	@param image 			boolean parameter to set if image persistence is calculated
	@param cokernel 		boolean parameter to set if cokernel persistence is calculated
	@param n_threads		number of threads to use in oineus
	@param upper_threshold	float, z-coordinate above which points are in the subcomplex 
	@param lower_threshold	float z-coordinate below which points are in the subcomplex


	@return kicr			oineus object which contains the kernel, image, cokernel persistence diagrams as required, can also calculate ones that weren't initially specificed
	"""
	print("started oineus_kernel_image_cokernel")
	sub = sub_complex(points, upper_threshold, lower_threshold)
	K, L = oineus_pair(points, sub)
	L = oineus.list_to_filtration_double(L)
	K = oineus.list_to_filtration_double(K)
	print("about to reduce")
	kicr = oineus.KerImCokReduced_double(K,L,params,False)
	print("reduced")
	dgm_0 = pd.DataFrame(np.hstack([kicr.codomain_diagrams().in_dimension(0), kicr.codomain_diagrams().index_diagram_in_dimension(0)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 0 diagram
	print("got dgm_0")
	dgm_0["birth simplex"]=dgm_0["birth simplex"].astype(int) #convert indices to int
	dgm_0["death simplex"]=dgm_0["death simplex"].astype(int) #convert indices to int
	dgm_1 = pd.DataFrame(np.hstack([kicr.codomain_diagrams().in_dimension(1), kicr.codomain_diagrams().index_diagram_in_dimension(1)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	print("got dgm_1")
	dgm_1["birth simplex"]=dgm_1["birth simplex"].astype(int) #convert indices to int
	dgm_1["death simplex"]=dgm_1["death simplex"].astype(int) #convert indices to int
	dgm_2 = pd.DataFrame(np.hstack([kicr.codomain_diagrams().in_dimension(2), kicr.codomain_diagrams().index_diagram_in_dimension(2)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	print("got dgm_2")
	dgm_2["birth simplex"]=dgm_2["birth simplex"].astype(int) #convert indices to int
	dgm_2["death simplex"]=dgm_2["death simplex"].astype(int) #convert indices to int
	print("finished oineus_kernel_image_cokernel")
	return kicr, dgm_0, dgm_1, dgm_2

def calculate_APF(dgm): 
	"""! Calcualte the APF from a diagram 
	@param dgm 		the diargam you want to calculate the APF for
	@return APF		the APF as a list of coordiantes
	"""
	lifetime = abs(dgm["death"] - dgm["birth"]) #get the lifetime of each point
	mean_age = (dgm["death"] + dgm["birth"])/2 #get the mean age
	APF = np.transpose(np.vstack([mean_age, lifetime])) #np array of the mean age and lifetime
	APF = APF[APF[:,0].argsort()] #sort the np array by ascending mean age
	for i in range(1, np.shape(APF)[0], 1):
			APF[i,1] = APF[i,1] + APF[i-1,1] #TODO: remove duplicates and only keep the last value of each one
	return pd.DataFrame(APF, columns = ["mean age", "lifetime"])

# Function to compute the persistent homology 
def compute():
	st.session_state.params = oineus.ReductionParams()
	if not st.session_state["manual_config"]:
		load_configuration_settings()#streamlit_functions.load_configuration_settings()
	if not st.session_state["maual_comp_config"]:
		load_computation_settings()#streamlit_functions.load_computation_settings()
	st.session_state.sample_indices = []
	st.session_state.dgms_0 = []
	st.session_state.dgms_1 = []
	st.session_state.dgms_2 = []
	st.session_state.APFs_0 = []
	st.session_state.APFs_1 = []
	st.session_state.APFs_2 = []
	st.session_state.atom_locations_list = []
	st.session_state.dcmps = [] 
	st.session_state.filts = []
	if "kernel" not in st.session_state:
		st.session_state.kernel = False
		st.session_state.params.kernel = False
	else:
		print("KERNEL is ", st.session_state.kernel)
		st.session_state.params.kernel = st.session_state.kernel
	if "image" not in st.session_state:
		st.session_state.image = False
		st.session_state.params.image = False
	else:
		st.session_state.params.image= st.session_state.image
	if "cokernel" not in st.session_state:
		st.session_state.cokernel = False
		st.session_state.params.cokernel = False
	else:
		st.session_state.params.cokernel = st.session_state.cokernel
	if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
		st.session_state.kicrs = []
	if st.session_state["kernel"]:
		st.session_state.kernel_dgms_0 = []
		st.session_state.kernel_dgms_1 = []
		st.session_state.kernel_dgms_2 = []
		st.session_state.kernel_APFs_0 = []
		st.session_state.kernel_APFs_1 = []
		st.session_state.kernel_APFs_2 = []
	if st.session_state["image"]:
		st.session_state.image_dgms_0 = []
		st.session_state.image_dgms_1 = []
		st.session_state.image_dgms_2 = []
		st.session_state.image_APFs_0 = []
		st.session_state.image_APFs_1 = []
		st.session_state.image_APFs_2 = []
	if st.session_state["cokernel"]:
		st.session_state.cokernel_dgms_0 = []
		st.session_state.cokernel_dgms_1 = []
		st.session_state.cokernel_dgms_2 = []
		st.session_state.cokernel_APFs_0 = []
		st.session_state.cokernel_APFs_1 = []
		st.session_state.cokernel_APFs_2 = []
	if st.session_state.sample_end == "Auto":
		st.session_state.sample_end = 1
	for s in range(st.session_state.sample_start, st.session_state.sample_end, st.session_state.sample_step):
		st.session_state.sample_index = s
		st.session_state.sample_indices.append(s)
		st.session_state.atom_locations = sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.atoms, st.session_state.radii)
		print("atom_locations is:")
		print(st.session_state.atom_locations)
		st.session_state.atom_locations_list.append(st.session_state.atom_locations)
		if st.session_state.params.kernel or st.session_state.params.image or st.session_state.params.cokernel:
			top_pt = max(st.session_state.atom_locations["z"])
			bot_pt = min(st.session_state.atom_locations["z"])
			height = abs(top_pt - bot_pt)
			if st.session_state.thickness == "Auto":
				ut= top_pt - 0.1*height
				lt = bot_pt + 0.1*height
			else:
				ut = top_pt - st.session_state.thickness*height
				lt = bot_pt + st.session_state.thickness*height
			st.session_state.kicr, st.session_state.dgm_0, st.session_state.dgm_1, st.session_state.dgm_2 = oineus_kernel_image_cokernel(st.session_state.atom_locations, st.session_state.params, ut, lt)
			st.session_state.kicrs.append(st.session_state.kicr)
			st.session_state.dgms_0.append(st.session_state.dgm_0)
			st.session_state.dgms_1.append(st.session_state.dgm_1)
			st.session_state.dgms_2.append(st.session_state.dgm_2)
			st.session_state.dcmps.append(st.session_state.kicr.decomposition_f)
			st.session_state.filts.append(st.session_state.kicr.fil_K)
			st.session_state.APFs_0.append(calculate_APF(st.session_state.dgm_0))#process.calculate_APF(st.session_state.dgm_0))
			st.session_state.APFs_1.append(calculate_APF(st.session_state.dgm_1))#process.calculate_APF(st.session_state.dgm_1))
			st.session_state.APFs_2.append(calculate_APF(st.session_state.dgm_2))#process.calculate_APF(st.session_state.dgm_2))
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				st.session_state.kicrs.append(st.session_state.kicr)
			if st.session_state["kernel"]:
				kernel_dgm_0 = pd.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(0), columns=["birth", "death"])
				kernel_dgm_1 = pd.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(1), columns=["birth", "death"])
				kernel_dgm_2 = pd.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(2), columns=["birth", "death"])
				print("kernel_dgm_2 is:")
				print(kernel_dgm_2)
				st.session_state.kernel_dgms_0.append(kernel_dgm_0)
				st.session_state.kernel_dgms_1.append(kernel_dgm_1)
				st.session_state.kernel_dgms_2.append(kernel_dgm_2)
				st.session_state.kernel_APFs_0.append(calculate_APF(kernel_dgm_0))#process.calculate_APF(kernel_dgm_0))
				st.session_state.kernel_APFs_1.append(calculate_APF(kernel_dgm_1))#process.calculate_APF(kernel_dgm_1))
				st.session_state.kernel_APFs_2.append(calculate_APF(kernel_dgm_2))#process.calculate_APF(kernel_dgm_2))
			if st.session_state["image"]:
				image_dgm_0 = pd.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(0), columns=["birth", "death"])
				image_dgm_1 = pd.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(1), columns=["birth", "death"])
				image_dgm_2 = pd.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(2), columns=["birth", "death"])
				st.session_state.image_dgms_0.append(image_dgm_0)
				st.session_state.image_dgms_1.append(image_dgm_1)
				st.session_state.image_dgms_2.append(image_dgm_2)
				st.session_state.image_APFs_0.append(calculate_APF(image_dgm_0))#	process.calculate_APF(image_dgm_0))
				st.session_state.image_APFs_1.append(calculate_APF(image_dgm_1))#process.calculate_APF(image_dgm_1))
				st.session_state.image_APFs_2.append(calculate_APF(image_dgm_2))#process.calculate_APF(image_dgm_2))
			if st.session_state["cokernel"]:
				cokernel_dgm_0 = pd.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				cokernel_dgm_1 = pd.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				cokernel_dgm_2 = pd.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				st.session_state.cokernel_dgms_0.append(cokernel_dgm_0)
				st.session_state.cokernel_dgms_1.append(cokernel_dgm_1)
				st.session_state.cokernel_dgms_2.append(cokernel_dgm_2)
				st.session_state.cokernel_APFs_0.append(calculate_APF(cokernel_dgm_0))#process.calculate_APF(cokernel_dgm_0))
				st.session_state.cokernel_APFs_1.append(calculate_APF(cokernel_dgm_1))#process.calculate_APF(cokernel_dgm_1))
				st.session_state.cokernel_APFs_2.append(calculate_APF(cokernel_dgm_2))#process.calculate_APF(cokernel_dgm_2))
		else:
			st.session_state.dcmp, st.session_state.filt, st.session_state.dgm_0, st.session_state.dgm_1, st.session_state.dgm_2 = 	oineus_process(st.session_state.atom_locations, st.session_state.params)
			st.session_state.dgms_0.append(st.session_state.dgm_0)
			st.session_state.dgms_1.append(st.session_state.dgm_1)
			st.session_state.dgms_2.append(st.session_state.dgm_2)
			st.session_state.APFs_0.append(calculate_APF(st.session_state.dgm_0))#process.calculate_APF(st.session_state.dgm_0))
			st.session_state.APFs_1.append(calculate_APF(st.session_state.dgm_1))#process.calculate_APF(st.session_state.dgm_1))
			st.session_state.APFs_2.append(calculate_APF(st.session_state.dgm_2))#process.calculate_APF(st.session_state.dgm_2))
			st.session_state.dcmps.append(st.session_state.dcmp) 
			st.session_state.filts.append(st.session_state.filt)

	if len(st.session_state.sample_indices) != len(st.session_state.dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.APFs_2):
		st.markdown("*WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
	if "kernel_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_2):
			st.markdown("KERNEL *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
			st.write(len(st.session_state.kernel_dgms_0))
			st.write(len(st.session_state.sample_indices))
	if "image_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.image_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.image_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.image_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_2):
			st.markdown("IMAGE *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.") 
	if "cokernel_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_2):
			st.markdown("COKERNEL *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
	st.session_state["processed_file"] =  st.session_state.file_path
	st.session_state.processed = True
	# print("Finished!")


###define various functions needed for later
def test():
	print(os.getcwd())
	st.session_state["manual_comp_config"] = True
	st.session_state.processed=True
	st.session_state.config_file = "../examples/structure-types.ini"
	st.session_state.file_path = "../examples/ZIF_test.xyz"
	st.session_state.config_name = "ZIF-TEST"
	st.session_state.comp_name = "ZIF-TEST"
	st.session_state.sample_start = 0
	st.session_state.sample_end = 2
	st.session_state.sample_step = 1
	st.session_state.repeat_x = 1
	st.session_state.repeat_y = 1
	st.session_state.repeat_z = 1
	st.session_state.thickness=0.1
	st.session_state.kernel = True
	st.session_state.image = True
	st.session_state.cokernel = True
	compute()