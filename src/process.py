##
# @internal
# @file process.py
# @brief Functions for processing the structures.
# Various wrappers and functions for obtaining filtrations, persistent homology objects and diagrams from point configurations.

from ase import io, Atoms
import numpy 
import pandas
import diode
#import dionysus
import oineus
import math
from colour import Color
from scipy.interpolate import interpn
from functools import cmp_to_key
 

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
	repeat_y = int(config.get(sconfiguratione, "REPEAT_Y")) #read repitition in y-axis
	print("repeating in y-axis: {}".format(repeat_x)) #print repitition in y-axis
	repeat_z = int(config.get(configuration, "REPEAT_Z")) #read repitition in z-axis
	print("repeating in z-axis: {}".format(repeat_z)) #print repitition in z-axis
	return atoms, radii, repeat_x, repeat_y, repeat_z

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
	dfpoints = pandas.DataFrame(numpy.column_stack([sample.get_chemical_symbols(), coord]), columns=["Atom", "x", "y", "z"]) #combine the atomic symbols with their location into a pandas.DataFrame
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
	dfpoints["w"] = numpy.select(conditions, choice_weights) #set the weights in the dataframe
	dfpoints["x"] = pandas.to_numeric(dfpoints["x"]) #ensure that everything is numeric and not string
	dfpoints["y"] = pandas.to_numeric(dfpoints["y"]) #ensure that everything is numeric and not string
	dfpoints["z"] = pandas.to_numeric(dfpoints["z"]) #ensure that everything is numeric and not string
	return dfpoints #return the dataframe of points


def weighted_alpha_diode(points):
	"""! Use diode to fill the weighted alpha shapes
	@param points 	pandas.DataFrame with columns 'x', 'y', 'z' for coordinates, and column 'w' with the weights.

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

def sub_complex(points : pandas.DataFrame, z_upper : float, z_lower : float):
	"""! Given the points, and the upper and lower thresholds in the 'z'-component. 

	@param points		pandas.DataFrame containing of the points.
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

def oineus_filtration(points : pandas.DataFrame, params : oineus.ReductionParams):
	"""! Given a set of points, compute the oineus.filtration of the alpha complex
 	
	@param points		pandas.DataFrame containing points and their weights
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
		
def oineus_pair(points : pandas.DataFrame, sub : list):
	"""! Given a set of points, and the points that are in the subset L, construct the complexes and map between them. The subcomplex L will consists of all simplices whose vertex are in the subset.

	@param points		pandas.DataFrame containing the points and their weights
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
 
def oineus_process(points : pandas.DataFrame, params : oineus.ReductionParams):
	"""! Given some points with weights, and the number of threads to use, obtain the persistent homology of the weighted alpha complex of these points, using oineus.
	
	@param points		pandas.DataFrame of the points, with colums 'x','y','z','w'
	@param params		oineus.ReudctionParams
	
	@return dgm_1 		pandas.DataFrame of the indexed dimension 1 diagram
	@return dgm_2		pandas.DataFrame of the indexed dimension 2 diagram
	"""
	filt = oineus_filtration(points, params) #get the filtration for oineus
	dcmp =  oineus.Decomposition(filt, False) #initialise the decomposition without cohomology
	dcmp.reduce(params) #reduce the matrix
	dgms = dcmp.diagram(filt) #initialise the diagrams
	dgm_1 = pandas.DataFrame(numpy.hstack([dgms.in_dimension(1), dgms.index_diagram_in_dimension(1)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	dgm_1["birth simplex"]=dgm_1["birth simplex"].astype(int) #convert indices to int
	dgm_1["death simplex"]=dgm_1["death simplex"].astype(int) #convert indices to int
	dgm_2 = pandas.DataFrame(numpy.hstack([dgms.in_dimension(2), dgms.index_diagram_in_dimension(2)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 2 diagram
	dgm_2["birth simplex"]=dgm_2["birth simplex"].astype(int) #convert indices to int
	dgm_2["death simplex"]=dgm_2["death simplex"].astype(int) #convert indices to int
	return dcmp, filt, dgm_1, dgm_2
	
def oineus_kernel_image_cokernel(points : pandas.DataFrame, params : oineus.ReductionParams, upper_threshold : float, lower_threshold : float):
	"""! Given points, and parameters for oineus, calculate the kernel/image/cokernel persistence as desired.

	@param points			pandas.DataFrame of the points, with columns 'x','y','z','w' corresponding to the coordinates and weights respectively
	@param kernel 			boolean parameter to set if kernel persistence is calculated
	@param image 			boolean parameter to set if image persistence is calculated
	@param cokernel 		boolean parameter to set if cokernel persistence is calculated
	@param n_threads		number of threads to use in oineus
	@param upper_threshold	float, z-coordinate above which points are in the subcomplex 
	@param lower_threshold	float z-coordinate below which points are in the subcomplex


	@return kicr			oineus object which contains the kernel, image, cokernel persistence diagrams as required, can also calculate ones that weren't initially specificed
	"""
	sub = sub_complex(points, upper_threshold, lower_threshold)
	K, L = oineus_pair(points, sub)
	L = oineus.list_to_filtration_float(L)
	K = oineus.list_to_filtration_float(K)
	kicr = oineus.KerImCokReduced_float(K,L,params,False)
	dgm_1 = pandas.DataFrame(numpy.hstack([kicr.codomain_diagrams().in_dimension(1), kicr.codomain_diagrams().index_diagram_in_dimension(1)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	dgm_1["birth simplex"]=dgm_1["birth simplex"].astype(int) #convert indices to int
	dgm_1["death simplex"]=dgm_1["death simplex"].astype(int) #convert indices to int
	dgm_2 = pandas.DataFrame(numpy.hstack([kicr.codomain_diagrams().in_dimension(2), kicr.codomain_diagrams().index_diagram_in_dimension(2)]), columns = ["birth", "death", "birth simplex", "death simplex"]) #get the indexed dimension 1 diagram
	dgm_2["birth simplex"]=dgm_2["birth simplex"].astype(int) #convert indices to int
	dgm_2["death simplex"]=dgm_2["death simplex"].astype(int) #convert indices to int
	return kicr, dgm_1, dgm_2

def calculate_APF(dgm): 
	"""! Calcualte the APF from a diagram 
	@param dgm 		the diargam you want to calculate the APF for
	@return APF		the APF as a list of coordiantes
	"""
	lifetime = abs(dgm["death"] - dgm["birth"]) #get the lifetime of each point
	mean_age = (dgm["death"] + dgm["birth"])/2 #get the mean age
	APF = numpy.transpose(numpy.vstack([mean_age, lifetime])) #numpy array of the mean age and lifetime
	APF = APF[APF[:,0].argsort()] #sort the numpy array by ascending mean age
	for i in range(1, numpy.shape(APF)[0], 1):
			APF[i,1] = APF[i,1] + APF[i-1,1] #TODO: remove duplicates and only keep the last value of each one
	return pandas.DataFrame(APF, columns = ["mean age", "lifetime"])


	