##
# @internal
# @file process.py
# @brief Functions for processing the structures.
# Various wrappers and functions for obtaining filtrations, persistent homology objects and diagrams from point configurations.

from ase import io, Atoms
import numpy 
import pandas
import diode
import dionysus
import oineus
import math
from colour import Color
from scipy.interpolate import interpn
from functools import cmp_to_key
import configparser

def read_configuration(structure_file : str, structure : str):
	"""! import a specified structure from a configuration file
	
	@param file_path    path to the file to use for configurations
	@param structure    name of the structure to use
	
	@result atoms, radii, repeat_x, repeat_y, repeat_z  the atom symbols, radii to use for the atoms, repetition in x-axis, repetition in y-axis, repetition in z-axis
	"""
	config = configparser.ConfigParser() 
	config.read(structure_file)
	atoms = [str(a).strip() for a in config.get(structure, "ATOMS").split(",")]
	radii = [float (r) for r in config.get(structure, "RADII").split(",")]
	print("Proceeding with the following:")
	print("Atoms:")
	for i, a in enumerate(atoms):
		print("{} with radius {}".format(a, radii[i]))
	repeat_x = int(config.get(structure, "REPEAT_X"))
	print("repeating in x-axis: {}".format(repeat_x))
	repeat_y = int(config.get(structure, "REPEAT_Y"))
	print("repeating in y-axis: {}".format(repeat_x))
	repeat_z = int(config.get(structure, "REPEAT_Z"))
	print("repeating in z-axis: {}".format(repeat_z))
	return atoms, radii, repeat_x, repeat_y, repeat_z

def read_sample(structure_file : str, sample_range : str):
	"""! import a specified sample range from a configuration file
	
	@param file_path    path to the file to use for configurations
	@param sample_range    name of the structure to use
	
	@result sample_start, sample_end, sample_step  first sample, last sample, step between each one
	"""
	config = configparser.ConfigParser() 
	structures = config.read(structure_file)
	sample_start = int(config.get(sample_range, "START"))
	sample_end = int(config.get(sample_range, "END"))
	sample_step = int(config.get(sample_range, "STEP"))
	return sample_start, sample_end, sample_step

def read_kernel_image_cokernel(structure_file : str, setting_name : str):
	"""! import settings for kernel/image/cokernel 
	@param file_path	path to the ini file containing the settings
	@param settings_name	name of the settings to use
	
	@result
	"""
	config = configparser.ConfigParser() 
	settings = config.read(structure_file)
	kernel = False
	image = False
	cokernel = False
	verbose = False
	n_threads = int(config.get(setting_name, "N_THREADS"))
	if (str(config.get(setting_name, "KERNEL")) == "TRUE") or (str(config.get(setting_name, "KERNEL")) == "True") or (str(config.get(setting_name, "KERNEL")) == "T"):
		kernel = True
	if (str(config.get(setting_name, "IMAGE")) == "TRUE") or (str(config.get(setting_name, "IMAGE")) == "True") or (str(config.get(setting_name, "IMAGE")) == "T"):
		image = True
	if (str(config.get(setting_name, "COKERNEL")) == "TRUE") or (str(config.get(setting_name, "COKERNEL")) == "True") or (str(config.get(setting_name, "COKERNEL")) == "T"):
		cokernel = True
	if (str(config.get(setting_name, "VERBOSE")) == "TRUE") or (str(config.get(setting_name, "VERBOSE")) == "True") or (str(config.get(setting_name, "VERBOSE")) == "T"):
		verbose = True
	return n_threads, kernel, image, cokernel, verbose


def load_atom_file(file_path : str, format : str, index = ":"):
	"""! load the file containing the initial configureation
	@param file_path file to load
	
	@return atoms     what we obtain from ase.io.read
	"""
	atoms = io.read(file_path, index = index, format = format)
	return atoms

def sample_at(atoms, sample_index, repeat_x : int, repeat_y : int, repeat_z : int, atom_list, radius_list):
	"""! Sample a structure at a particular time, with cell repetitions.
	@param atoms		initial configuration
	@param sample_index 	time to sample at
	@param repeat_x 	repetition in x dir
	@param repeat_y 	repetition in y dir
	@param repeat_z 	repetition in z dir
	@parm atom_list 	list of atoms in the config
	@param radius_list	list of radii to use for the atoms

	@return points		data frame of the points including the repetitions
	"""
	sample  = atoms[sample_index].repeat((repeat_x, repeat_y, repeat_z))
	coord = sample.get_positions()
	cell = sample.get_cell()
	data = numpy.column_stack([sample.get_chemical_symbols(), coord])
	dfpoints = pandas.DataFrame(data, columns=["Atom", "x", "y", "z"])
	atoms_found = dfpoints["Atom"].unique()
	remove = []
	print(atom_list)
	for a in atoms_found:
		if a not in atom_list:
			remove.append([a])
			print("We are going to remove the following atom types:", a)
	if len(remove) != 0:
		dfpoints = dfpoints[dfpoints["Atom"].isin(atom_list)]
	conditions = [(dfpoints["Atom"]==atom_list[i]) for i in range(len(atom_list))]
	choice_weights = [radius_list[i]**2 for i in range(len(radius_list))]
	dfpoints["w"] = numpy.select(conditions, choice_weights)
	dfpoints["x"] = pandas.to_numeric(dfpoints["x"])
	dfpoints["y"] = pandas.to_numeric(dfpoints["y"])
	dfpoints["z"] = pandas.to_numeric(dfpoints["z"])
	points = dfpoints[["x", "y", "z", "w", "Atom"]]
	return points


def weighted_alpha_diode(points):
	"""! Use diode to fill the weighted alpha shapes
	@param points 	pandas.DataFrame with columns 'x', 'y', 'z' for coordinates, and column 'w' with the weights.

	@return weighted alpha shape from diode.
	"""
	return diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())

def persistent_homology_filt_dionysus(simplices : list):
	"""! Get the filtration and persistence module from a list of simplices (using dionysus), and remove any simplicies in dimensions above 3.
	@param simplices 	list of simplices from dionysus

	@return filt, m 	the dionysus filtration, dionysus persistent homology
	"""
	restricted_simps = []
	for s in simplices:
		if len(s[0]) <= 4:
			restricted_simps.append(s)
	filt = dionysus.Filtration(restricted_simps)
	m = dionysus.homology_persistence(filt, progress = True)
	return filt, m

def extract_diagrams_from_dionysus(filt, m):
	"""! Given a dionysus fitlration and persistent homology, extract the persistence diagrams
	@param filt		dionysus filtration
	@param m		dionysus persistent homology

	@return dgm		nx2 numpy.array of diagram
	"""
	dionysus_diagrams= dionysus.init_diagrams(m, filt)
	dgms = []
	for i in range(len(dionysus_diagrams)):
		dgm_i = numpy.empty((0,2), float)
		for p in dionysus_diagrams[i]:
			dgm_i = numpy.append(dgm_i, [p.birth, p.death], axis=0)
		dgms.append(dgm_i)
	return dgms
	
	
	
def persistent_homology_diagrams_from_points_dionysus(points):#extract the persistent homology of the frame 
	"""! Given a set of points with weights, return the persistence diagrams from dionysus.
	@param points	pandas.DataFrame of points, columns 'x','y','z' coordinates, and column 'w' the weights.

	@return dgms 	dionysus persistence diagrams
	"""
	# Calculate the alpha shape using Diode
	simplices = weighted_alpha_diode(points)
	# Generate the filtration using Dionysus
	filt, m = persistent_homology_filt_dionysus(simplices)
	# Get persistent homologypersistent_homology_dionysus(filtration)
	d_dgms = extract_diagrams_from_dionysus(filt, m)
	dgms = extract_diagrams_from_dionysus
	return dgms

def aggregate_diagrams(dgms): #aggregate the diagrams into one big diagram
	"""! Given a list of diagrams, combine these into a single diagram.

	@param dgms		list of diagrams as numpy.arrays.

	@return dgm		numpy.array of the combined diagram.
	"""
	dgm = numpy.concatenate(dgms, axis=0)
	return dgm
	#for i in range(len(births)):
	#	for j in range(len(births[i])):
	#		birth.append(births[i][j])  
	#		death.append(deaths[i][j]) 
	#birth = numpy.array(birth)
	#death = numpy.array(death)
	#return birth, death

def calculate_APF(dgm): #TODO: decide what to do with points at infinity
	"""! Calcualte the APF from a diagram 
	@param dgm 		the dionysus diargam you want to calculate the APF for

	@return APF		the APF as a list of coordiantes
	"""
	lifetime = dgm[:,1] - dgm[:,0]
	mean_age = (dgm[:,1] + dgm[:,0])/2
	APF = numpy.transpose(numpy.vstack([mean_age, lifetime]))
	APF = APF[APF[:,0].argsort()]
	for i in range(1, numpy.shape(APF)[0], 1):
			APF[i,1] = APF[i,1] + APF[i-1,1]
	return APF

#def oineus_compare_long(x, y):
#	"""! 
#	"""
#	if len(x[0]) == len(y[0]):
#		for i in range(len(x[0])):
#			if x[0][i] < y[0][i]:
#				return -1
#			return 1
#	elif len(x[0]) < len(y[0]):
#		return -1
#	else:
#		return 1


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
		
def oineus_pair(points : pandas.DataFrame, sub : list):
	"""! Given a set of points, and the points that are in the subset L, construct the complexes and map between them. The subcomplex L will consists of all simplices whose vertex are in the subset.

	@param points		pandas.DataFrame containing the points and their weights
	@param sub			a list containing the indices of the points on which we construct the subcomplex

	@return K			list of simplices for the entire complex, as needed by oineus
	@return L			list of simplices for the subcomplex, as needed by oineus
	@return L_to_K		list which tells you how to map the simplices in L to the simplices in K
	"""
	simplices = diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())
	for i in range(len(simplices)):
		simplices[i] = [sorted(simplices[i][0]), simplices[i][1]]
	simplices = sorted(simplices, key=cmp_to_key(oineus_compare))
	K = []
	L = []
	L_to_K = []
	id_L = 0
	K_to_L = [-1 for i in range(len(simplices))]
	for i,s in enumerate(simplices):
		if len(s[0])==1:
			K.append([i,s[0],s[1]])
			if sub[s[0][0]]==True:
				K_to_L[i]= id_L            
				L.append([id_L, [K_to_L[s[0][0]]], s[1]])
				L_to_K.append(i)
				id_L +=1
		else:
			K.append([i, s[0],s[1]])
			sub_complex = True
			for v in s[0]:
				if sub[v] == False:
					sub_complex=False
					break
			if sub_complex == True:
				verts = [K_to_L[v] for v in s[0]]
				L.append([id_L, verts, s[1]])
				L_to_K.append(i)
				K_to_L[i] = id_L
				id_L +=1
	return K, L, L_to_K


def sub_complex(points : pandas.DataFrame, z_upper : float, z_lower : float):
	"""! Given the points, and the upper and lower thresholds in the 'z'-component. 

	@param points		pandas.DataFrame containing of the points.
	@param z_upper		float giving the upper threshold, any point above this is in the subcomplex
	@param z_lower		float giving the lower threshold, any point below this is in the subcomplex

	@return sub_comp	list containing the indices of the points on which we build the subcomplex
	"""
	sub_comp = []
	for i in range(points.shape[0]):
		if (points["z"][i] >= z_upper) or (points["z"][i]    <= z_lower):
			sub_comp.append(True)
		else:
			sub_comp.append(False)
	return sub_comp     

def kernel_image_cokernel(points : pandas.DataFrame, kernel : bool, image : bool, cokernel : bool, n_threads : int, upper_threshold : float, lower_threshold : float):
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
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	sub = sub_complex(points, upper_threshold, lower_threshold)
	K, L, L_to_K = oineus_pair(points, sub)
	kicr = oineus.compute_kernel_image_cokernel_reduction(K, L, L_to_K, params)
	return kicr