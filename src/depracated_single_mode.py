# import sys
# sys.path.insert(0, '..')
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

st.session_state.params = oineus.ReductionParams()

def read_configuration(configuration_file : str, configuration : str):
	"""! import a specified structure from a configuration file
	
	@param configuration_file    path to the file to use for configurations
	@param configuration    name of the structure to use
	
	@result atoms, radii, repeat_x, repeat_y, repeat_z  the atom symbols, radii to use for the atoms, repetition in x-axis, repetition in y-axis, repetition in z-axis
	"""
	config = configparser.ConfigParser() #to be able to read the configuration file we need a parse
	config.read(configuration_file) #load the file containing the structures
	atoms = [str(a).strip() for a in config.get(configuration, "ATOMS").split(",")] #read the atoms in the selected structure
	radii = [float (r) for r in config.get(configuration, "RADII").split(",")] #get the radii to use for the atoms
	print("Proceeding with the following:\nAtoms:") #print what has been loaded into the program, should probably display this in the GUI so that it can be confirmed
	for i, a in enumerate(atoms):
		print("{} with radius {}".format(a, radii[i])) #print atom and the radius
	repeat_x = int(config.get(configuration, "REPEAT_X")) #read repitition in x-axis
	print("repeating in x-axis: {}".format(repeat_x)) #print repitition in x-axis
	repeat_y = int(config.get(configuration, "REPEAT_Y")) #read repitition in y-axis
	print("repeating in y-axis: {}".format(repeat_x)) #print repitition in y-axis
	repeat_z = int(config.get(configuration, "REPEAT_Z")) #read repitition in z-axis
	print("repeating in z-axis: {}".format(repeat_z)) #print repitition in z-axis
	sample_index = int(config.get(configuration, "SAMPLE_INDEX")) #read repitition in z-axis
	print("sampling at index: {}".format(sample_index)) #print repitition in z-axis
	return atoms, radii, repeat_x, repeat_y, repeat_z, sample_index

def read_computation_settings(settings_file : str, settings_name):
	config = configparser.ConfigParser() #to be able to read the configuration file we need a parse
	config.read(settings_file) #load the file containing the structures
	try:
		n_threads = int(config.get(settings_name,"N_THREADS"))
	except:
		n_threads = 4
	try:
		save_plots = bool(config.get(settings_name,"SAVE_PLOTS"))
	except:
		save_plots = False
	try:
		if config.get(settings_name,"KERNEL") == "TRUE" or config.get(settings_name,"KERNEL") == "True" or config.get(settings_name,"KERNEL") == "true" or config.get(settings_name,"KERNEL") == "t" or config.get(settings_name,"KERNEL") == "T":
			kernel = True
		else:
			kernel = False
	except:
		kernel = False
	try:
		if config.get(settings_name,"IMAGE") == "TRUE" or config.get(settings_name,"IMAGE") == "True" or config.get(settings_name,"IMAGE") == "true" or config.get(settings_name,"IMAGE") == "t" or config.get(settings_name,"IMAGE") == "T":
			image = True
		else:
			image = False
	except:
		image = False
	try:
		if config.get(settings_name,"COKERNEL") == "TRUE" or config.get(settings_name,"COKERNEL") == "True" or config.get(settings_name,"COKERNEL") == "true" or config.get(settings_name,"COKERNEL") == "t" or config.get(settings_name,"COKERNEL") == "T":
			cokernel = True
		else:
			cokernel = False
	except:
		cokernel = False
	try:
		thickness = float(config.get(settings_name, "THICKNESS"))
	except:
		thickness = "Auto"
	try:
		lower_threshold = float(config.get(settings_name, "LOWER_THRESHOLD"))
	except:
		lower_threshold = 1 # DONT KNOW A GOOD DEFAULT
	
	return n_threads, save_plots, kernel, image, cokernel, thickness, lower_threshold


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

#def persistent_homology_filt_dionysus(simplices : list): #no longer needed
#	"""! Get the filtration and persistence module from a list of simplices (using dionysus), and remove any simplicies in dimensions above 3.
#	@param simplices 	list of simplices from dionysus
#
#	@return filt, m 	the dionysus filtration, dionysus persistent homology
#	"""
#	restricted_simps = []
#	for s in simplices:
#		if len(s[0]) <= 4:
#			restricted_simps.append(s)
#	filt = dionysus.Filtration(restricted_simps)
#	m = dionysus.homology_persistence(filt, progress = True)
#	return filt, m

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
	#L_to_K = []
	#id_L = 0
	#K_to_L = [-1 for i in range(len(simplices))]
	for i,s in enumerate(simplices):
		if len(s[0])==1:
			if sub[s[0][0]]==True:
				#K_to_L[i]= id_L            
				L.append([s[0], s[1]])
				#L_to_K.append(i)
				#id_L +=1
			else:
				not_L.append([s[0],s[1]])
		else:
			sub_complex = True
			for v in s[0]:
				if sub[v] == False:
					sub_complex=False
					break
			if sub_complex == True:
				#verts = [K_to_L[v] for v in s[0]]
				L.append([s[0], s[1]])
				#L_to_K.append(i)
			  	#K_to_L[i] = id_L
				#id_L +=1
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
	L = oineus.list_to_filtration_float(L)
	K = oineus.list_to_filtration_float(K)
	print("about to reduce")
	kicr = oineus.KerImCokReduced_float(K,L,params,False)
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
	return kicr, dgm_0,dgm_1, dgm_2

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
# import streamlit_functions
def load_configuration_settings():
	st.session_state.atoms, st.session_state.radii, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.sample_index = read_configuration(st.session_state["config_file"], st.session_state["config_name"])

def load_computation_settings():
	st.session_state.n_threads, st.session_state.save_plots, st.session_state.kernel, image, st.session_state.cokernel, st.session_state.thickness, st.session_state.lower_threshold = read_computation_settings(st.session_state["comp_file"], st.session_state["comp_name"])



st.session_state.loaded = False
st.session_state.processed = False
st.session_state.plotted = False
st.session_state.params = oineus.ReductionParams()

st.header("ToMA Single Mode")
file_path = st.text_input("Intial structure file:",key="file_path") #specify initial structure file 
file_format = st.text_input("File format:", key="file_format",placeholder="Auto") #specify format of the initial strutcure file


st.header("Configuration settings")
manual_config = st.checkbox("Manually specify configuration", key="manual_config")#manually set configuration
if not manual_config:
	st.session_state.config_file = st.text_input("Configuration file:", key="configuration_file")
	st.session_state.config_name = st.text_input("Configuration name:", key="configuration_name")
else:
	st.session_state.atoms = st.text_input("Atoms:", key="atoms_input")
	st.session_state.radii = st.text_input("Radii:", key="radii_input")
	st.session_state.sample_index = st.text_input("Sample index", key="sample_index_input")
	st.session_state.repeat_x = st.text_input("Repitition in x-axis:", key="repeat_x_input")
	st.session_state.repeat_y = st.text_input("Repitition in y-axis:", key="repeat_y_input")
	st.session_state.repeat_z = st.text_input("Repitition in z-axis:", key="repeat_z_input")


st.markdown("Computation settings")
manual_compute = st.checkbox("Manually specify settings for the computations (i.e number of threds, and if you want to compute kernel/image/cokernel)", key="maual_comp_config")
if not manual_compute:
	same_config_file = st.checkbox("The computation settings are in the same configuration file.", key="same_config_file")
	if not same_config_file:
		st.session_state.comp_file = st.text_input("Configuration file:", key="comp_config_file")
	else:
		st.session_state.comp_file = st.session_state.config_file
	st.session_state.comp_name = st.text_input("Configuration name:", key="comp_config_name")
else:
	st.session_state.params.kernel = st.checkbox("Compute kernel persistence", key="kernel_check")
	st.session_state.params.image = st.checkbox("Compute image persistence", key="image_check")
	st.session_state.params.cokernel = st.checkbox("Compute cokernel persistence", key="cokernel_check")
	st.text_input("Select thickness of top and bottom layer:", key="thickness_input", placeholder="Automatic detection")
	st.session_state.n_threads = st.text_input("Select number of threads to use:", key="n_threads_input", placeholder="4")
	if st.session_state.n_threads == "":
		st.session_state.params.n_threads = 4
	else:
		st.session_state.params.n_threads = int(st.session_state.n_threads)
	st.write("Number of threads is", st.session_state.params.n_threads)

	if st.session_state["thickness_input"] == "":
		st.session_state.thickness = "Auto"
	else:
		st.session_state.thickness= float(st.session_state["thickness_input"])

st.header("Plot generation")
st.markdown("Please select which of the following plots you would like to generate.")
pd_checks = st.columns(3)
with pd_checks[0]:
	st.checkbox("Dimension 0 Persistence Diagram", key="pd0")
with pd_checks[1]:
	st.checkbox("Dimension 1 Persistence Diagram", key="pd1")
with pd_checks[2]:
	st.checkbox("Dimension 2 Persistence Diagram", key="pd2")
pd0 = st.session_state["pd0"]
pd1 = st.session_state["pd1"]
pd2 = st.session_state["pd2"]

apf_checks = st.columns(3)
with apf_checks[0]:
	st.checkbox("Dimension 0 Accumulated Persistence Function", key="apf0")
with apf_checks[1]:
	st.checkbox("Dimension 1 Accumulated Persistence Function", key="apf1")
with apf_checks[2]:
	st.checkbox("Dimension 2 Accumulated Persistence Function", key="apf2")
apf0 = st.session_state["apf0"]
apf1 = st.session_state["apf1"]
apf2 = st.session_state["apf2"]

st.write("The file selected to analyse is ", file_path)
if file_format == "":
	st.write("File format will automatically detected.")
else:
	st.write("File format is", file_format)

	
def compute():
	if not st.session_state["manual_config"]:
		load_configuration_settings()
	if not st.session_state["maual_comp_config"]:
		load_computation_settings()
	st.session_state.atom_locations  = sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.atoms, st.session_state.radii)
	st.session_state.processed=True
	top_pt = max(st.session_state.atom_locations["z"])
	bot_pt = min(st.session_state.atom_locations["z"])
	height = abs(top_pt - bot_pt)
	if st.session_state.thickness == "Auto":
		ut= top_pt - 0.1*height
		lt = bot_pt + 0.1*height
	else:
		ut = top_pt - st.session_state.thickness*height
		lt = bot_pt + st.session_state.thickness*height
	st.session_state.kicr, st.session_state.dgm_1, st.session_state.dgm_2 =  oineus_kernel_image_cokernel(st.session_state.atom_locations, st.session_state.params, ut, lt)
	st.session_state.processed = True

def save():
	print("saving")
	if st.session_state.processed:
		# Standard persistence diagrams
		if pd0:
			fig_pd_0 = plot_PD(st.session_state.dgm_0, st.session_state.file_path+" PD0")
			fig_pd_0.write_image("PD0.png")
		if pd1:
			fig_pd_1 = plot_PD(st.session_state.dgm_1, st.session_state.file_path+" PD1") 
			fig_pd_1.write_image("PD1.png")
			# Kernel/Image/Cokernel diagrams if they were computed
			if st.session_state.kernel or st.session_state.image or st.session_state.cokernel:
				fig_kic_pd_1 = plot_kernel_image_cokernel_PD(
					st.session_state.kicr, 
					1, 
					True,
					st.session_state.kernel,
					st.session_state.image, 
					st.session_state.cokernel,
					st.session_state.file_path+" codomain/kernel/image/cokernel dimension 1"
				)
				fig_kic_pd_1.write_image("KIC_PD1.png")
		if pd2:
			fig_pd_2 = plot_PD(st.session_state.dgm_2, st.session_state.file_path+" PD2")
			fig_pd_2.write_image("PD2.png")
			# Kernel/Image/Cokernel diagrams if they were computed
			if st.session_state.kernel or st.session_state.image or st.session_state.cokernel:
				fig_kic_pd_2 = plot_kernel_image_cokernel_PD(
					st.session_state.kicr,
					2,
					True, 
					st.session_state.kernel,
					st.session_state.image,
					st.session_state.cokernel,
					st.session_state.file_path+" codomain/kernel/image/cokernel dimension 2"
				)
				fig_kic_pd_2.write_image("KIC_PD2.png")
		# Accumulated persistence functions
		if apf0:
			apf_0 = calculate_APF(st.session_state.dgm_0)
			fig_apf_0 = plot_APF(apf_0, st.session_state.file_path+" APF0")
			fig_apf_0.write_image("APF0.png")
		if apf1:
			apf_1 = calculate_APF(st.session_state.dgm_1)
			fig_apf_1 = plot_APF(apf_1, st.session_state.file_path+" APF1")
			fig_apf_1.write_image("APF1.png")
		if apf2:
			apf_2 = calculate_APF(st.session_state.dgm_2)
			fig_apf_2 = plot_APF(apf_2, st.session_state.file_path+" APF2")
			fig_apf_2.write_image("APF2.png")

st.button("Process", key="process", on_click=compute)
st.button("Save", key="save", on_click=save)

if st.session_state.processed:
	st.write("DONE")



# 			params.n_threads = int(values_main["n_threads"])
# 			if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 				if values_main["upper_threshold"] == "Automatic":
# 					upper_threshold = math.floor(max(atom_locations["z"])) 
# 				else:
# 					upper_threshold = float(values_main["upper_threshold"])
# 				if values_main["lower_threshold"] == "Automatic":
# 						lower_threshold = math.ceil(min(atom_locations["z"]))
# 				else:
# 					lower_threshold = float(values_main["lower_threshold"])
# 				params.kernel = values_main["kernel"]
# 				params.image = values_main["image"]
# 				params.cokernel = values_main["cokernel"]
# 				kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(atom_locations, params, upper_threshold, lower_threshold)
# 				if values_main["APF1"]:
# 					APF_1 = calculate_APF(dgm_1)
# 				if values_main["APF2"]:
# 					APF_2 = calculate_APF(dgm_2)	
# 			else:
# 				dcmp, filt, dgm_1, dgm_2 = oineus_process(atom_locations, params)
# 				if values_main["APF1"]:
# 					APF_1 = calculate_APF(dgm_1)
# 				if values_main["APF2"]:
# 					APF_2 = calculate_APF(dgm_2)
# 			processed = True
   
# 		if event_main == "Plot":
# 			if processed == False:
# 				sg.popup_error("File not processed, please Process it.")
# 			else:
# 				if values_main['PD1'] == True:
# 					if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 						try:
# 							fig_kic_pd_1 = plot_kernel_image_cokernel_PD(kicr, 1, True, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" codmain/kernel/image/cokernel dimension 1 sample "+str(s))
# 							fig_kic_pd_1.show()
# 						except:
# 							print("Kernel/image/cokernel persistence has not been calculated.")
# 					else:
# 						fig_pd_1 = plot_PD(dgm_1, file_path+" PD1 sample "+str(s))
# 						fig_pd_1.show()
# 				if values_main['PD2'] == True:
# 					if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 						try:
# 							fig_kic_pd_2 = plot_kernel_image_cokernel_PD(kicr, 2, True, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" codmain/kernel/image/cokernel dimension 2 sample "+str(s))
# 							fig_kic_pd_2.show()
# 						except:
# 							print("Kernel/image/cokernel persistence has not been calculated.")
# 					else:
# 						fig_pd_2 = plot_PD(dgm_2, file_path+" PD2 sample "+str(s))
# 						fig_pd_2.show()
# 				if values_main['APF1'] == True:
# 					try:
# 						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 sample "+str(s))
# 						fig_apf_1.show()
# 					except:
# 						APF_1 = calculate_APF(dgm_1)
# 						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 sample "+str(s))
# 						fig_apf_1.show()
# 				if values_main['APF2'] == True:
# 					try:
# 						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 sample "+str(s))
# 						fig_apf_2.show()	
# 					except:
# 						APF_2 = calculate_APF(dgm_2)
# 						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 sample "+str(s))
# 						fig_apf_2.show()
# 				plotted = True
				
		
# 		if event_main == "Save":
# 			if processed == False:
# 				sg.popup_error("File not processed, please process it.")
# 			else:
# 				dir = os.path.dirname(file_path)
# 				file_name = os.path.splitext(os.path.split(file_path)[1])[0]
# 				try:
# 					pd.DataFrame(np.column_stack([births[1], deaths[1]]), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_PD_1.csv", header = None)
# 					pd.DataFrame(APF_1, columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_APF_1.csv", header=  None)
# 				except:
# 					print("Persistence diagram in dimension 1 is empty.")
# 				try:
# 					pd.DataFrame(np.column_stack([births[2], deaths[2]]), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_PD_2.csv", header = None)
# 					pd.DataFrame(APF_2, columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_APF_2.csv", header = None)
# 				except:
# 					print("Persistence diagram in dimension 2 is empty.")
# 				if values_main["kernel"]:
# 					pd.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_1.csv", header = None)
# 					pd.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_2.csv", header = None)
# 				if values_main["image"]:
# 					pd.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_image_PD_1.csv", header = None)
# 					pd.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_image_PD_2.csv", header = None)
# 				if values_main["cokernel"]:
# 					pd.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_1.csv", header = None)
# 					pd.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_2.csv", header = None)
# 				if plotted == False:
# 					if values_main['APF1'] == True:
# 						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 "+str(s))
# 					if values_main['APF2'] == True:
# 						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 "+str(s))
# 					if values_main['PD1'] == True:
# 						fig_pd_1 = plot_PD(births[1], deaths[1], file_path+" PD1 "+str(s))
# 						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 							fig_kic_pd_1 = plot_kernel_image_cokernel_PD(kicr, 1, True, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" codmain/kernel/image/cokernel dimension 1"+str(s))
# 					if values_main['PD2'] == True:
# 						fig_pd_2 = plot_PD(births[2], deaths[2], file_path+" PD2 "+str(s))
# 						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 							fig_kic_pd_2 = plot_kernel_image_cokernel_PD(kicr, 2, True, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" codmain/kernel/image/cokernel dimension 2"+str(s))
# 				else:
# 					if values_main["PD1"]:
# 						try:
# 							fig_pd_1.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_PD_1.png")
# 						except:
# 							print("saving failed")
# 						try: #values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 							fig_kic_pd_1.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_kic_PD_1.png")
# 						except:
# 							print("saving failed")
# 					if values_main["PD2"]:
# 						try:
# 							fig_pd_2.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_PD_2.png")
# 						except:
# 							print("saving failed")
# 						try:#if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
# 							fig_kic_pd_2.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_kic_PD_2.png")
# 						except:
# 							print("saving failed")
# 					if values_main["APF1"]:
# 						try:
# 							fig_apf_1.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_APF_1.png")
# 						except:
# 							print("saving failed")
# 					if values_main["APF2"]:
# 						try:
# 							fig_apf_2.write_image(dir+"/"+file_name+"_sample_"+str(s)+"_APF_2.png")
# 						except:
# 							print("saving failed")
		
# 		if event_main == "Visualisation":
# 			if processed == False:
# 				sg.popup_error("File not processed, please process it.")
# 			else:
# 				dfVis = generate_visulisation_df(dgm_1, dcmp.r_data, filt, atom_locations, atoms)
# 				visualisation_table_layout = [[sg.Table(values=dfVis.values.tolist(), headings=dfVis.columns.values.tolist(), auto_size_columns=True, num_rows = 50, display_row_numbers=True, selected_row_colors="red on yellow", enable_events=True)], [sg.Button("Display selected", key="Display"), sg.Checkbox("Plot neighbours", key="neighbours", font=ButtonFont)]]
# 				visualisation_table_window = sg.Window("AMA: 1-cycle representatives table", visualisation_table_layout, resizable=True)
# 				while True:
# 					event_visualisation, value_visualisation = visualisation_table_window.read()
# 					if event_visualisation == "Display":
# 						vis = generate_display(atom_locations, dfVis, value_visualisation[0][0], filt, value_visualisation["neighbours"])
# 						vis.show()
# 					if event_visualisation == "Exit" or event_visualisation == sg.WIN_CLOSED:
# 						break
# 		if event_main == "Exit" or event_main == sg.WIN_CLOSED:
# 			window_main.close()
# 			entry_window()
# 			break
# 		if event_main == "Multi":
# 			window_main.close()
# 			multi_mode()
# 			break
# 		if event_main == "Batch":
# 			window_main.close()
# 			batch_mode()
# 			break
# 		if event_main == "Quit":
# 			window_main.close()
# 			break