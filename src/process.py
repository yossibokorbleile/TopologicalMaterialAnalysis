##
# @internal
# @file process.py
# @brief Functions for processing the structures.
# Various wrappers and functions for obtaining filtrations, persistent homology objects and diagrams from point configurations.

from ase import io
import numpy 
import pandas
import diode
import dionysus
import oineus
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from colour import Color
from matplotlib import cm
from matplotlib.colors import LogNorm
from scipy.interpolate import interpn
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets  import RectangleSelector
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


def load_atom_file(file_path : str, format : str):
	"""! load the file containing the initial configureation
	@param file_path file to load
	
	@return xyz     what we obtain from ase.io.read
	"""
	atoms = io.read(file_path, index = ":", format = format)
	return atoms

def sample_at(xyz, sample_time, repeat_x : int, repeat_y : int, repeat_z : int, atom_list, radius_list):
	"""! Sample a structure at a particular time, with cell repetitions.
	@param xyz		initial configuration
	@param sample_time 	time to sample at
	@param repeat_x 	repetition in x dir
	@param repeat_y 	repetition in y dir
	@param repeat_z 	repetition in z dir
	@parm atom_list 	list of atoms in the config
	@param radius_list	list of radii to use for the atoms

	@return points		data frame of the points including the repetitions
	"""
	xyz[sample_time] = xyz[sample_time].repeat((repeat_x, repeat_y, repeat_z))
	coord = xyz[sample_time].get_positions()
	cell = xyz[sample_time].get_cell()
	data = numpy.column_stack([xyz[sample_time].get_chemical_symbols(), coord])
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
	return diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())

def persistent_homology_filt_dionysus(simplices, dim):
	restricted_simps = []
	for s in simplices:
		if len(s[0]) <= 4:
			restricted_simps.append(s)
	filt = dionysus.Filtration(restricted_simps)
	m = dionysus.homology_persistence(filt, progress = True)
	return filt, m

def extract_diagrams_dionysus(filt, m):
	return dionysus.init_diagrams(m, filt)
	
def get_birth_death(dionysus_diagrams):    
	births = []
	deaths = []
	for i in range(len(dionysus_diagrams)):
		births_i = []
		deaths_i = []
		for p in dionysus_diagrams[i]:
			births_i.append(p.birth)
			deaths_i.append(p.death)    
		births.append(numpy.array(births_i))
		deaths.append(numpy.array(deaths_i))
	return births, deaths
	
def persistent_homology_diagrams_from_points_dionysus(points):#extract the persistent homology of the frame 
	# Calculate the alpha shape using Diode
	simplices = weighted_alpha_diode(points)
	# Generate the filtration using Dionysus
	filt, m = persistent_homology_filt_dionysus(simplices)
	# Get persistent homologypersistent_homology_dionysus(filtration)
	dgms = extract_diagrams_dionysus(filt, m)
	return dgms

def aggregate_diagrams(births, deaths): #aggregate the diagrams into one big diagram
	birth = []
	death = []
	for i in range(len(births)):
		for j in range(len(births[i])):
			birth.append(births[i][j])  
			death.append(deaths[i][j]) 
	birth = numpy.array(birth)
	death = numpy.array(death)
	return birth, death

def calculate_APF(births, deaths): #TODO: decide what to do with points at infinity
	lifetime = deaths - births
	mean_age = (deaths - births)/2
	APF = numpy.transpose(numpy.vstack([mean_age, lifetime]))
	APF = APF[APF[:,0].argsort()]
	for i in range(1, numpy.shape(APF)[0], 1):
			APF[i,1] = APF[i,1] + APF[i-1,1]
	return APF

def oineus_compare_long(x, y):
	if len(x[0]) == len(y[0]):
		for i in range(len(x[0])):
			if x[0][i] < y[0][i]:
				return -1
			return 1
	elif len(x[0]) < len(y[0]):
		return -1
	else:
		return 1


def oineus_compare(x, y):
	if len(x[0]) == len(y[0]):
		if x[0][0] <= y[0][0]:
			return -1
		else:
			return 1
	elif len(x[0]) < len(y[0]):
		return -1
	else:
		return 1
		
def oineus_pair(points, sub):
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


def sub_complex(points, z_upper, z_lower):
	sub_comp = []
	for i in range(points.shape[0]):
		if (points["z"][i] >= z_upper) or (points["z"][i]    <= z_lower):
			sub_comp.append(True)
		else:
			sub_comp.append(False)
	return sub_comp     

def kernel_image_cokernel(points, kernel, image, cokernel, n_threads, upper_threshold, lower_threshold):
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	sub = sub_complex(points, upper_threshold, lower_threshold)
	K, L, L_to_K = oineus_pair(points, sub)
	kicr = oineus.compute_kernel_image_cokernel_reduction(K, L, L_to_K, params)
	return kicr
#	if params.kernel:
#		pd_1 = kicr.kernel_diagrams().in_dimension(1)
#		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_kernel_PD_1.csv")
#		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_kernel_PD_1.png")
#		matplotlib.pyplot.close()
#		pd_2 = kicr.kernel_diagrams().in_dimension(2)
#		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_kernel_PD_2.csv")
#		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_kernel_PD_2.png")
#		matplotlib.pyplot.close()
#		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
#		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
#		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_kernel_APF_1.csv")
#		fig = plot_APF(APF_1, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_kernel_APF_1.png")
#		matplotlib.pyplot.close()
#		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_kernel_APF_2.csv")
#		fig = plot_APF(APF_2, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_kernel_APF_2.png")
#		matplotlib.pyplot.close()
#	if params.image:
#		pd_1 = kicr.image_diagrams().in_dimension(1)
#		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_image_PD_1.csv")
#		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_image_PD_1.png")
#		matplotlib.pyplot.close()
#		pd_2 = kicr.kernel_diagrams().in_dimension(2)
#		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_image_PD_2.csv")
#		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_image_PD_2.png")
#		matplotlib.pyplot.close()
#		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
#		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
#		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_image_APF_1.csv")
#		fig = plot_APF(APF_1, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_image_APF_1.png")
#		matplotlib.pyplot.close()
#		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_image_APF_2.csv")
#		fig = plot_APF(APF_2, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_image_APF_2.png")
#		matplotlib.pyplot.close()
#	if params.cokernel:
#		pd_1 = kicr.kernel_diagrams().in_dimension(1)
#		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_cokernel_PD_1.csv")
#		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_cokernel_PD_1.png")
#		matplotlib.pyplot.close()
#		pd_2 = kicr.kernel_diagrams().in_dimension(2)
#		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_cokernel_PD_2.csv")
#		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
#		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_cokernel_PD_2.png")
#		matplotlib.pyplot.close()
#		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
#		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
#		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_cokernel_APF_1.csv")
#		fig = plot_APF(APF_1, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_cokernel_APF_1.png")
#		matplotlib.pyplot.close()			
#		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_cokernel_APF_2.csv")
#		fig = plot_APF(APF_2, 'blue')
#		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_cokernel_APF_2.png")
#		matplotlib.pyplot.close()