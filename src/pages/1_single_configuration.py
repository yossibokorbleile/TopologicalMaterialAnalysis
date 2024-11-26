##
# @internal
# @file multi_mode.py
# @brief Streamlit page for analysing a multiple samples.
# @version 0.1
# @date July 2024

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
from scipy.interpolate import interpn
import plotly.express as px
import plotly.graph_objects as go




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
	return atoms, radii, repeat_x, repeat_y, repeat_z

def read_computation_settings(settings_file : str, settings_name):
	"""! import a specified settings from a file
	
	@param settings_file    path to the file to use for configurations
	@param settings_name    name of the structure to use
	
	@result n_threads, save_plots, kernel, image, cokernel, thickness  the number of threads for Oineus to use, wether plots should be saved, compute kernel persistence, compute image persistence, compute cokernel persistence, thickness settings
	"""
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
		if mode_config.get(settings_name,"IMAGE") == "TRUE" or mode_config.get(settings_name,"IMAGE") == "True" or mode_config.get(settings_name,"IMAGE") == "true" or mode_config.get(settings_name,"IMAGE") == "t" or mode_config.get(settings_name,"IMAGE") == "T":
			image = True
		else:
			image = False
	except:
		image = False
	try:
		if mode_config.get(settings_name,"COKERNEL") == "TRUE" or mode_config.get(settings_name,"COKERNEL") == "True" or mode_config.get(settings_name,"COKERNEL") == "true" or mode_config.get(settings_name,"COKERNEL") == "t" or mode_config.get(settings_name,"COKERNEL") == "T":
			cokernel = True
		else:
			cokernel = False
	except:
		cokernel = False
	try:
		thickness = float(mode_config.get(settings_name, "THICKNESS"))
	except:
		thickness = "Auto"
	return n_threads, save_plots, kernel, image, cokernel, thickness


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
# import streamlit_functions
def load_configuration_settings():
	st.session_state.atoms, st.session_state.radii, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z = read_configuration(st.session_state["config_file"], st.session_state["config_name"])

def load_computation_settings():
	st.session_state.n_threads, st.session_state.save_plots, st.session_state.kernel, image, st.session_state.cokernel, st.session_state.thickness = load_computation_settings(st.session_state["comp_file"], st.session_state["comp_name"])



def plot_APF(APF : np.array, name : str):
	"""! Plot an accumulated persistence function
	
	@param APF - np.array with 2 columns of coordinates which define the APF
	@param name - title for the plot
	
	@result a plotly.express figure
	"""
	
	fig = px.line(x=APF["mean age"], y=APF["lifetime"], labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_APFs(APFs : list, APF_names : list, fig_name : str):#, APF_colour, APF_label):
	"""! Plot a set accumulated persistence function, with automatic colour differentiation.
	
	@param APFs - accumlated persistence functions to plot
	
	@result a matplotlib figure
	"""
	assert len(APFs) == len(APF_names)
	fig = go.Figure(labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=fig_name)
	last_pt = math.ceil(max([APFs[i][-1,0] for i in range(len(APFs))])*1.1)
	for i in range(len(APFs)):
		APFs[i] = np.vstack([APFs[i], [last_pt, APFs[i][-1,1]]])
	for i in range(len(APFs)):
		fig.add_trace(go.Scatter(x=APFs[i][:,0], y=APFs[i][:,1], mode="lines", name=APF_names[i]))
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_PD(dgm, name : str):
	"""! Plot a persistence diagram, with a specific colour
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.
	
	@param dgm 	- 	pandas.DataFrame of the diagram
	@param name    - name to use as title of the plot
	
	@result a plotly.express figure
	"""
	try:
		max_val = max(dgm["death"][dgm["death"] != math.inf])
	except:
		max_val = max(dgm["birth"])
	birth = []
	death = []
	inf_fin = []
	fig = go.Figure()
	for i in range(dgm.shape[0]):
		if dgm["death"].iloc[i] == math.inf:
			birth.append(dgm["birth"].iloc[i])
			death.append(max_val*1.1)
			inf_fin.append("inf")
		else:
			birth.append(dgm["birth"].iloc[i])
			death.append(dgm["death"].iloc[i])
			inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="inf_fin", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig


def plot_PDs(dgms, name : str):
	"""! Plot several persistence diagrams, with  automatic colour choices
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.

	@param dgms - list of diagrams
	@param name - title to use for the plot

	@results a plotly.express figure
	"""
	birth = []
	death = []
	samp = []
	inf_fin = []
	vals = []
	for dgm in dgms:
		dgm_vals = []
		for d in dgm["death"]:
			if d != math.inf:
				dgm_vals.append(d)
		if len(dgm_vals) !=0:
			vals.append(max(dgm_vals))
	if len(vals) != 0:
		max_val = max(vals)
	else:
		max_val = max([max(b) for b in dgm["birth"] for dgm in dmgs])
	fig = go.Figure()
	for i in range(len(dgms)):
		for j in range(len(dgms[i]["death"])):
			if dgms[i]["death"].iloc[j] == math.inf:
				birth.append(dgs[i]["birth"].iloc[j])
				death.append(max_val*1.1)
				samp.append(str(i))
				inf_fin.append("inf")
			else:
				birth.append(dgms[i]["birth"].iloc[j])
				death.append(dgms[i]["death"].iloc[j])
				samp.append(str(i))
				inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "sample":samp, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", color="sample", symbol="inf_fin", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig


def plot_kernel_image_cokernel_PD(kicr, d : int, codomain : bool, kernel : bool, image : bool, cokernel : bool, name : str):
	"""! Plot kernel, image, cokernel on same figure
	@param kicr 	oineus::KerImCokReduced 
	@param d	 	the dimension to extract (either 1 or 2)
	@param kernel	bool to plot kernel
	@param image	bool to plot image
	@param cokernel	bool to plot cokernel
	@return figu	figure with the chosen PD diagrams
	"""
	fig = go.Figure()
	max_val = -math.inf
	if codomain:
		codomain_pd = kicr.codomain_diagrams().in_dimension(d)
		print("codomain diagram has {} points".format(codomain_pd.shape[0]))
		if math.inf in codomain_pd[:,1] and max_val < max([d for d in codomain_pd[:,1] if d !=math.inf]):
			max_val = max(codomain_pd[:,1])
	if kernel:
		kernel_pd = kicr.kernel_diagrams().in_dimension(d)
		print("kernel diagram has {} points".format(kernel_pd.shape[0]))
		if math.inf in kernel_pd[:,1] and max_val < max([d for d in kernel_pd[:,1] if d !=math.inf]):
			max_val = max(kernel_pd[:,1])
	if image:
		image_pd = kicr.image_diagrams().in_dimension(d)
		print("image diagram has {} points".format(image_pd.shape[0]))
		if math.inf in image_pd[:,1]  and max_val < max([d for d in image_pd[:,1] if d !=math.inf]):
			max_val = max(image_pd[:,1])
	if cokernel:
		print("cokernel diagram has {} points".format(cokernel_pd.shape[0]))
		cokernel_pd = kicr.cokernel_diagrams().in_dimension(d)
		if math.inf in cokernel_pd[:,1] and max_val < max([d for d in cokernel_pd[:,1] if d !=math.inf]):
			max_val =  max(cokernel_pd[:,1])
	birth = []
	death = []
	pt_type = []
	inf_fin = []
	if codomain:
		for i in range(codomain_pd.shape[0]):
			if codomain_pd[i,1] == math.inf:
				birth.append(codomain_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("codomain")
				inf_fin.append("inf")
			else:
				birth.append(codomain_pd[i,0])
				death.append(codomain_pd[i,1])
				pt_type.append("codomain")
				inf_fin.append("fin")
	if kernel:
		for i in range(kernel_pd.shape[0]):
			if kernel_pd[i,1] == math.inf:
				birth.append(kernel_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("kernel")
				inf_fin.append("inf")
			else:
				birth.append(kernel_pd[i,0])
				death.append(kernel_pd[i,1])
				pt_type.append("kernel")
				inf_fin.append("fin")
	if image:
		for i in range(image_pd.shape[0]):
			if image_pd[i,1] == math.inf:
				birth.append(image_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("image")
				inf_fin.append("inf")
			else:
				birth.append(image_pd[i,0])
				death.append(image_pd[i,1])
				pt_type.append("image")
				inf_fin.append("fin")
	if cokernel:
		for i in range(cokernel_pd.shape[0]):
			if cokernel_pd[i,1] == math.inf:
				birth.append(cokernel_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("cokernel")
				inf_fin.append("inf")
			else:
				birth.append(cokernel_pd[i,0])
				death.append(cokernel_pd[i,1])
				pt_type.append("cokernel")
				inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "pt_type":pt_type, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="inf_fin", color="pt_type", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def write_files(dgm, file_path, save_plots=False, plot_name = ""):
	"""! Save a digram as a csv file, with the option to save a plot as well. 
	
	@param dgm 	the diagram to be saved, and plotted
	@param dir 	directory in which to save the files
	@param name	name of the structure
	@param save_plot save the plots as well 
	"""
	
	if len(dgm["birth"]) == 0:
		print("Diagram {} is empty.".format(plot_name))
		dgm.to_csv(file_path+".csv")
		pd.DataFrame([[0, 0]], columns=["mean age", "lifetime"]).to_csv(file_path.replace("PD","APF")+".csv")
	else:
		dgm.to_csv(file_path+".csv")
		APF = calculate_APF(dgm)
		pd.DataFrame(APF, columns=["mean age", "lifetime"]).to_csv(file_path.replace("PD","APF")+".csv")
		if save_plots:
			fig = plot_PD(dgm, plot_name)
			fig.write_html(file_path+".html")
			fig = plot_APF(APF, plot_name)
			fig.write_html(file_path.replace("PD", "APF")+".html")
	return True


st.header("Single Configuration Mode")
comp_tab, plot_tab, vis_tab = st.tabs(["Computation", "Plots", "Visuatlisation"]) #create tabs for the various parts
st.session_state.mode = "multi"


###define various functions needed for later
def test():
	print(os.getcwd())
	st.session_state["maual_comp_config"] = True
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

# function to generate plots
def generate_plots():
	st.session_state.fig_pds_0 = []
	st.session_state.fig_pds_1 = []
	st.session_state.fig_pds_2 = []
	st.session_state.fig_apfs_0 = []
	st.session_state.fig_apfs_1 = []
	st.session_state.fig_apfs_2 = []
	if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
		st.session_state.fig_kic_pds_0 = []
		st.session_state.fig_kic_pds_1 = []
		st.session_state.fig_kic_pds_2 = []
		st.session_state.fig_pds_1 = []
		st.session_state.fig_pds_2 = []
		st.session_state.fig_kernel_apfs_0 = []
		st.session_state.fig_kernel_apfs_1 = []
		st.session_state.fig_kernel_apfs_2 = []
		st.session_state.fig_image_apfs_0 = []
		st.session_state.fig_image_apfs_1 = []
		st.session_state.fig_image_apfs_2 = []
		st.session_state.fig_cokernel_apfs_0 = []
		st.session_state.fig_cokernel_apfs_1 = []
		st.session_state.fig_cokernel_apfs_2 = []
		st.session_state.fig_apfs_0 = []
		st.session_state.fig_apfs_1 = []
		st.session_state.fig_apfs_2 = []
	if "kernel" not in st.session_state:
		st.session_state["kernel"] = False
	if "image" not in st.session_state:
		st.session_state["image"] = False
	if "cokernel" not in st.session_state:
		st.session_state["cokernel"] = False
	for i,s in enumerate(st.session_state.sample_indices):
		if st.session_state["pd0"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_0.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 0, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 0 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 0.")
			st.session_state.fig_pds_0.append(plot_PD(st.session_state.dgms_0[i], st.session_state.file_path+" PD0 sample "+str(s)))
		if st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_1.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 1, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 1 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 1.")
			st.session_state.fig_pds_1.append(plot_PD(st.session_state.dgms_1[i], st.session_state.file_path+" PD1 sample "+str(s)))
		if st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_2.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 2, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 2 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 2.")
			st.session_state.fig_pds_2.append(plot_PD(st.session_state.dgms_2[i], st.session_state.file_path+" PD2 sample "+str(s)))
		if st.session_state["apf0"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_0.append(plot_APF(st.session_state.kernel_APFs_0[i], file_path+" kernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 0.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_0.append(plot_APF(st.session_state.image_APFs_0[i], file_path+" image APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 0.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_0.append(plot_APF(st.session_state.cokernel_APFs_0[i], file_path+" cokernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 0.")
			st.session_state.fig_apfs_0.append(plot_APF(st.session_state.APFs_0[i], file_path+" APF0 sample "+str(s)))
		if st.session_state["apf1"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_1.append(plot_APF(st.session_state.kernel_APFs_1[i], file_path+" kernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 1.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_1.append(plot_APF(st.session_state.image_APFs_1[i], file_path+" image APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 1.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_1.append(plot_APF(st.session_state.cokernel_APFs_1[i], file_path+" cokernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 1.")
			st.session_state.fig_apfs_1.append(plot_APF(st.session_state.APFs_1[i], file_path+" APF1 sample "+str(s)))
		if st.session_state["apf2"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_2.append(plot_APF(st.session_state.kernel_APFs_2[i], file_path+" kernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 2.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_2.append(plot_APF(st.session_state.image_APFs_2[i], file_path+" image APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 2.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_2.append(plot_APF(st.session_state.cokernel_APFs_2[i], file_path+" cokernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 2.")
			st.session_state.fig_apfs_2.append(plot_APF(st.session_state.APFs_2[i], file_path+" APF2 sample "+str(s)))
	st.session_state.plots_generated = True

#function to display plots
def display_plots():
	if "plots_generated" not in st.session_state or not st.session_state["plots_generated"]:
		generate_plots()
	plot_tab.write(st.session_state.sample_end)
	plot_tab.write(st.session_state.sample_indices)
	for i,s in enumerate(st.session_state.sample_indices):
		if st.session_state["pd0"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_0[i])
			# else: 
			plot_tab.plotly_chart(st.session_state.fig_pds_0[i])
		if st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_1[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_1[i])
		if st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_2[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_2[i])
		if st.session_state["apf0"]:
			if st.session_state["kernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_kernel_apfs_0[i])
				except:
					print("No kernel APF in dimension 0 to plot for sample index ", s)
			if st.session_state["image"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_image_apfs_0[i])
				except:
					print("No image APF in dimension 0 to plot for sample index ", s)
			if st.session_state["cokernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_cokernel_apfs_0[i])
				except:
					print("No cokernel APF in dimension 0 to plot for sample index ", s)
			plot_tab.plotly_chart(st.session_state.fig_apfs_0[i])
		if st.session_state["apf1"]:
			if st.session_state["kernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_kernel_apfs_1[i])
				except:
					print("No kernel APF in dimension 1 to plot for sample index ", s)
			if st.session_state["image"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_image_apfs_1[i])
				except:
					print("No image APF in dimension 1 to plot for sample index ", s)
			if st.session_state["cokernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_cokernel_apfs_1[i])
				except:
					print("No cokernel APF in dimension 1 to plot for sample index ", s)
			plot_tab.plotly_chart(st.session_state.fig_apfs_1[i])
		if st.session_state["apf2"]:
			if st.session_state["kernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_kernel_apfs_2[i])
				except:
					print("No kernel APF in dimension 2 to plot for sample index ", s)
			if st.session_state["image"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_image_apfs_2[i])
				except:
					print("No image APF in dimension 2 to plot for sample index ", s)
			if st.session_state["cokernel"]:
				try:
					plot_tab.plotly_chart(st.session_state.fig_cokernel_apfs_2[i])
				except:
					print("No cokernel APF in dimension 2 to plot for sample index ", s)
			plot_tab.plotly_chart(st.session_state.fig_apfs_2[i])


# function to save plots
def save_plots():
	dir_name = os.path.dirname(file_path)
	file_name = os.path.splitext(os.path.split(file_path)[1])[0]
	print(dir_name)
	print(file_name)
	if "plots_generated" not in st.session_state or not st.session_state["plots_generated"]:
		generate_plots()
	for i, s in enumerate(st.session_state.sample_indices):
		st.session_state.sample_index = s
		if st.session_state["pd0"] == True:
			try:
				st.session_state.fig_pds_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_PD_0.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_PD_0.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_pds_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_PD_0.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_pds_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_PD_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_PD_0.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_pds_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_PD_0.png")
		if st.session_state["pd1"] == True:
			try:
				st.session_state.fig_pds_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_PD_1.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_PD_1.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_pds_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_PD_1.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_pds_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_PD_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_PD_1.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_pds_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_PD_1.png")
		if st.session_state["pd2"] == True:
			try:
				st.session_state.fig_pds_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_PD_2.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_PD_2.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_pds_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_2.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_PD_2.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_pds_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_PD_2.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_PD_2.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_pds_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_PD_0.png")
		if st.session_state["apf0"] == True:
			try:
				st.session_state.fig_apfs_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_APF_0.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_APF_0.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_apfs_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_APF_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_APF_0.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_apfs_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_APF_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_APF_0.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_apfs_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_APF_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_APF_0.png")
		if st.session_state["apf1"] == True:
			try:
				st.session_state.fig_apfs_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_APF_1.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_APF_1.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_apfs_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_APF_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_APF_1.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_apfs_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_APF_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_APF_1.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_apfs_1[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_APF_1.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_APF_1.png")
		if st.session_state["apf2"] == True:
			try:
				st.session_state.fig_apfs_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_APF_2.png")
			except:
				print("Error saving "+file_name+"_sample_"+str(s)+"_APF_2.png")
			if st.session_state.params.kernel:
				try:
					st.session_state.fig_kernel_apfs_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_kernel_APF_2.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_kernel_APF_2.png")
			if st.session_state.params.image:
				try:
					st.session_state.fig_image_apfs_2[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_image_APF_2.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_image_APF_2.png")
			if st.session_state.params.cokernel:
				try:
					st.session_state.fig_cokernel_apfs_0[i].write_image(dir_name+"/"+file_name+"_sample_"+str(s)+"_cokernel_APF_0.png")
				except:
					print("Error saving "+file_name+"_sample_"+str(s)+"_cokernel_APF_0.png")


### lets set up the computation tab
file_path = comp_tab.text_input("Initial structure file:",key="file_path") #specify initial structure file 
file_format = comp_tab.text_input("File format:", key="file_format",placeholder="Auto") #specify format of the initial strutcure file
if "processed_file" in st.session_state:
	if file_path != st.session_state["processed_file"]:
		st.session_state.processed = False

comp_tab.header("Configuration settings")
manual_config = comp_tab.checkbox("Manually specify configuration", key="manual_config")#manually set configuration
if not manual_config:
	st.session_state.config_file = comp_tab.text_input("Configuration file:", key="configuration_file")
	st.session_state.config_name = comp_tab.text_input("Configuration name:", key="configuration_name")
else:
	st.session_state.atoms = comp_tab.text_input("Atoms:", key="atoms_input")
	st.session_state.radii = comp_tab.text_input("Radii:", key="radii_input")


comp_tab.markdown("Computation settings")
manual_compute = comp_tab.checkbox("Manually specify settings for the computations (i.e number of threds, and if you want to compute kernel/image/cokernel)", key="maual_comp_config")
if not manual_compute:
	same_config_file = comp_tab.checkbox("The computation settings are in the same configuration file.", key="same_config_file")
	if not same_config_file:
		st.session_state.comp_file = comp_tab.text_input("Configuration file:", key="comp_config_file")
	else:
		st.session_state.comp_file = st.session_state.config_file
	st.session_state.comp_name = comp_tab.text_input("Configuration name:", key="comp_config_name")
else:
	st.session_state.sample_start = comp_tab.text_input("Sample start", key="sample_start_input")
	st.session_state.sample_end = comp_tab.text_input("Sample end", key="sample_end_input")
	st.session_state.sample_step = comp_tab.text_input("Sample step", key="sample_step_input")
	st.session_state.repeat_x = comp_tab.text_input("Repitition in x-axis:", key="repeat_x_input")
	st.session_state.repeat_y = comp_tab.text_input("Repitition in y-axis:", key="repeat_y_input")
	st.session_state.repeat_z = comp_tab.text_input("Repitition in z-axis:", key="repeat_z_input")
	st.session_state.kernel = comp_tab.checkbox("Compute kernel persistence", key="kernel_check")
	st.session_state.image = comp_tab.checkbox("Compute image persistence", key="image_check")
	st.session_state.cokernel = comp_tab.checkbox("Compute cokernel persistence", key="cokernel_check")
	st.session_state.thickness = comp_tab.text_input("Select thickness of top and bottom layer:", key="thickness_input", placeholder="Automatic detection")
	st.session_state.n_threads = comp_tab.text_input("Select number of threads to use:", key="n_threads_input", placeholder="4")
	if st.session_state.n_threads == "":
		st.session_state.params.n_threads = 4
	else:
		st.session_state.params.n_threads = int(st.session_state.n_threads)
	comp_tab.markdown(f"Number of threads is".format(st.session_state.params.n_threads))

	if st.session_state["thickness_input"] == "":
		st.session_state.thickness = "Auto"
	else:
		st.session_state.thickness= float(st.session_state["thickness_input"])

comp_tab.markdown(f"The file selected to analyse is "+file_path)	
if file_format == "":
	comp_tab.markdown("File format will automatically detected.")
else:
	comp_tab.markdown("File format is", file_format)


comp_tab.button("Process", key="process", on_click=compute)


### Set up the plot tab
plot_tab.header("Plot generation")
if "processed" in st.session_state and st.session_state["processed"]:
	plot_tab.markdown("Please select which of the following plots you would like to generate.")
	pd_checks = plot_tab.columns(3)
	with pd_checks[0]:
		st.checkbox("Dimension 0 Persistence Diagram", key="pd0")
	with pd_checks[1]:
		st.checkbox("Dimension 1 Persistence Diagram", key="pd1")
	with pd_checks[2]:
		st.checkbox("Dimension 2 Persistence Diagram", key="pd2")

	apf_checks = plot_tab.columns(3)
	with apf_checks[0]:
		st.checkbox("Dimension 0 Accumulated Persistence Function", key="apf0")
	with apf_checks[1]:
		st.checkbox("Dimension 1 Accumulated Persistence Function", key="apf1")
	with apf_checks[2]:
		st.checkbox("Dimension 2 Accumulated Persistence Function", key="apf2")
else:
	plot_tab.markdown("Persistent homology has not been computed, so the plots can not be generated yet. Please proces the file, and then return to this tab.")

plot_buttons = plot_tab.columns(3)
with plot_buttons[0]:
	st.button("Generate plots", key="generate_plots", on_click=generate_plots)
with plot_buttons[1]:
	st.button("Display plots", key="display_plots", on_click=display_plots)
with plot_buttons[2]:
	st.button("Save plots", key="save_Plots", on_click=save_plots)

###Set up visualisation tab
vis_tab.markdown("In this tab, you can select *representatives* of homology classes to visualise.")	
vis_tab.checkbox("Visualisation", key="visualisation")
if 'selected_row' not in st.session_state:
	st.session_state.selected_row = None
if "processed" in st.session_state and st.session_state["processed"] and st.session_state["visualisation"] and not st.session_state.params.kernel and not st.session_state.params.image and not st.session_state.params.cokernel:
	vis_tab.write(st.session_state.params)
	vis_tab.write(st.session_state.dcmps)
	selected_sample = vis_tab.radio("Selection which sample you want to explore", st.session_state.sample_indices)
	st.session_state.selected_sample_index = st.session_state.sample_indices.index(selected_sample)
	vis_tab.write(st.session_state.selected_sample_index)
	# st.session_state.dimension = vis_tab.radio("What dimension cycles do you want to visualise:", [1, 2])
	st.session_state.dimension = 1
	vis_tab.checkbox("Display neighbouring atoms.", key="neighbours")
	if st.session_state.dimension == 1:
		vis_tab.markdown("Visulisation of representative 1-cycles.")
		st.session_state.vis_dgm = st.session_state.dgms_1[st.session_state.selected_sample_index]
	elif st.session_state.dimension == 2:
		vis_tab.markdown("Visualising representatives of 2-cycles is still underdevelopment, reverting to visualisation 1-cycles.")
		st.session_state.dimension = 1
		# vis_tab.markdown("Visulisation of representative 2-cycles.")
		# st.session_state.vis_dgm = st.session_state.dgm_2
	st.session_state.dfVis = visualisation.generate_visulisation_df(st.session_state.vis_dgm, st.session_state.dcmps[st.session_state.selected_sample_index].r_data, st.session_state.filts[st.session_state.selected_sample_index], st.session_state.atom_locations_list[st.session_state.selected_sample_index], st.session_state.atoms)
	to_display = ["birth", "death", "lifetime"]
	for a in st.session_state.atoms:
		to_display.append(a)
	viz = vis_tab.dataframe(st.session_state.dfVis[to_display], on_select="rerun")
	if viz.selection.rows == []:
		st.session_state.selected_row = None
	else:
		st.session_state.selected_row = viz.selection.rows
	# vis_tab.write("Selected Row:")
	# vis_tab.write(st.session_state.selected_row)
	if st.session_state.selected_row != None:
		for cycle_id in st.session_state.selected_row:
			vis_tab.plotly_chart(visualisation.generate_display(st.session_state.atom_locations, st.session_state.dgm_1, cycle_id, st.session_state.filt, neighbours = st.session_state["neighbours"]))
elif st.session_state.params.kernel or st.session_state.params.image or st.session_state.params.cokernel:
	vis_tab.markdown("Visulation of kernel/image/cokernel persistent homology is not yet available.")
else:
	vis_tab.markdown("Persistent homology has not been computed, so representative cycles are unknow. Please proces the file, and then return to this tab.")

st.button("test", key="test", on_click=test)

if "params" in st.session_state:
	print(st.session_state.params)

