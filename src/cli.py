##
# @internal
# @file cli.py
# contains the functions for the command line interface
from process import *
from plots import *
import os
import matplotlib
#from mpi4py import MPI 
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#nprocs = comm.Get_size()


def batch_mode(parent_dir : str, file_ext : str, format : str, structure_file : str, structure : str, sample_range : str, save_all : bool = False):
	"""! Batch mode for CLI usage. 
	
	@param parent_dir 	base directory from which to initialise the search
	@param file_ext		extension used for the files containing the initial configurations, files with any other extension will be ignored
	@param structure_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param structure 	specific structure to load from the configuration file
	@param sample_range		specific sample settings to use, should be listed in the structure file	
	"""
	
	atom_list, radii, repeat_x, repeat_y, repeat_z = read_configuration(structure_file, structure)
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []
	sample_settings = read_sample(structure_file, sample_range)
	sample_every = list(range(sample_settings[0], sample_settings[1], sample_settings[2]))
	config_files = []
	for root, dirs, files in os.walk(parent_dir):
		for f in files:
			if f.endswith(file_ext):
				config_files.append(os.path.join(root, f))
	print("Have found the following configuration files:")
	for f in config_files:
		print(f)
	#if rank == 0:
	#	distances = np.empty([n,n], dtype= float)
	#	data = np.arange(0,len(config_files),1)
	#	print(data)
	#	# determine the size of each sub-task
	#	ave, res = divmod(len(config_files), nprocs)
	#	counts = [ave + 1 if p < res else ave for p in range(nprocs)]
	#	# determine the starting and ending indices of each sub-task
	#	starts = [sum(counts[:p]) for p in range(nprocs)]
	#	ends = [sum(counts[:p+1]) for p in range(nprocs)]
	#	# converts data into a list of arrays 
	#	data = [data[starts[p]:ends[p]] for p in range(nprocs)]
	#else:
	#	data = None
	#data = comm.scatter(data, root=0)
	#print('Process {} has data:'.format(rank), data)
	#for config in data:
	for config in config_files:
		print("Looking at {}".format(config))
		atoms = load_atom_file(config, format)
		dir = os.path.dirname(config)
		if not os.path.exists(os.path.join(dir, "PD1")):
			os.mkdir(os.path.join(dir, "PD1"))
		if not os.path.exists(os.path.join(dir, "PD2")):
			os.mkdir(os.path.join(dir, "PD2"))
		if not os.path.exists(os.path.join(dir, "APF1")):
			os.mkdir(os.path.join(dir, "APF1"))
		if not os.path.exists(os.path.join(dir, "APF2")):
			os.mkdir(os.path.join(dir, "APF2"))
		config_name = os.path.splitext(os.path.split(config)[1])[0]
		births_1 = []
		deaths_1 = []
		births_2 = []
		deaths_2 = []
		APFs_1 = []
		APFs_2 = []
		for s in sample_every:
			print("looking at sample {}".format(s))
			points = sample_at(atoms, s, repeat_x, repeat_y, repeat_z, atom_list, radii)
			dgms = persistent_homology_diagrams_from_points_dionysus(points)
			births_s, deaths_s = get_birth_death(dgms)
			APF_1 = calculate_APF(births_s[1], deaths_s[1])
			APF_2 = calculate_APF(births_s[2], deaths_s[2])
			if save_all:
				pandas.DataFrame(numpy.column_stack([births_s[1], deaths_s[1]])).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.csv")
				fig = plot_PD(births_s[1], deaths_s[1], 'blue')
				matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(numpy.column_stack([births_s[2], deaths_s[2]])).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.csv")
				fig = plot_PD(births_s[2], deaths_s[2], 'blue')
				matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.csv")
				fig = plot_APF(APF_1, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.csv")
				fig = plot_APF(APF_2, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.png")
				matplotlib.pyplot.close()
			births_1.append(births_s[1])
			deaths_1.append(deaths_s[1])
			births_2.append(births_s[2])
			deaths_2.append(deaths_s[2])
			APFs_1.append(APF_1)
			APFs_2.append(APF_2)
		fig = plot_PDs(births_1, deaths_1)
		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_all_samples_PD_1.png")
		matplotlib.pyplot.close()
		fig = plot_PDs(births_2, deaths_2)
		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_all_samples_PD_2.png")
		matplotlib.pyplot.close()
		fig = plot_APFs(APFs_1)
		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_all_samples_AFP_1.png")
		matplotlib.pyplot.close()
		fig = plot_APFs(APFs_2)
		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_all_samples_AFP_2.png")
		matplotlib.pyplot.close()

def batch_mode_kernel(parent_dir : str, file_ext : str, format : str, structure_file : str, structure : str, sample_range : str, kernel_settings : str, save_all : bool = False):
	"""! Batch mode kernel for CLI usage. 
	
	@param parent_dir 	base directory from which to initialise the search
	@param file_ext		extension used for the files containing the initial configurations, files with any other extension will be ignored
	@param structure_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param structure 	specific structure to load from the configuration file
	@param sample_range		specific sample settings to use, should be listed in the structure file	
	"""
	import oineus
	atom_list, radii, repeat_x, repeat_y, repeat_z = read_configuration(structure_file, structure)
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []
	sample_settings = read_sample(structure_file, sample_range)
	sample_every = list(range(sample_settings[0], sample_settings[1], sample_settings[2]))
	config_files = []
	for root, dirs, files in os.walk(parent_dir):
		for f in files:
			if f.endswith(file_ext):
				config_files.append(os.path.join(root, f))
	print("Have found the following configuration files:")
	for f in config_files:
		print(f)
	for config in config_files:
		print("looking at {}".format(config))
		atoms = load_atom_file(config, format)
		dir = os.path.dirname(config)
		if not os.path.exists(os.path.join(dir, "PD1")):
			os.mkdir(os.path.join(dir, "PD1"))
		if not os.path.exists(os.path.join(dir, "PD2")):
			os.mkdir(os.path.join(dir, "PD2"))
		if not os.path.exists(os.path.join(dir, "APF1")):
			os.mkdir(os.path.join(dir, "APF1"))
		if not os.path.exists(os.path.join(dir, "APF2")):
			os.mkdir(os.path.join(dir, "APF2"))
		config_name = os.path.splitext(os.path.split(config)[1])[0]
		births_1 = []
		deaths_1 = []
		births_2 = []
		deaths_2 = []
		APFs_1 = []
		APFs_2 = []
		for s in sample_every:
			print("looking at sample {}".format(s))
			points = sample_at(atoms, s, repeat_x, repeat_y, repeat_z, atom_list, radii)
			dgms = persistent_homology_diagrams_from_points_dionysus(points)
			births_s, deaths_s = get_birth_death(dgms)
			APF_1 = calculate_APF(births_s[1], deaths_s[1])
			APF_2 = calculate_APF(births_s[2], deaths_s[2])
			if save_all:
				pandas.DataFrame(numpy.column_stack([births_s[1], deaths_s[1]])).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.csv")
				fig = plot_PD(births_s[1], deaths_s[1], 'blue')
				matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(numpy.column_stack([births_s[2], deaths_s[2]])).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.csv")
				fig = plot_PD(births_s[2], deaths_s[2], 'blue')
				matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.csv")
				fig = plot_APF(APF_1, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.csv")
				fig = plot_APF(APF_2, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.png")
				matplotlib.pyplot.close()
			births_1.append(births_s[1])
			deaths_1.append(deaths_s[1])
			births_2.append(births_s[2])
			deaths_2.append(deaths_s[2])
			APFs_1.append(APF_1)
			APFs_2.append(APF_2)
			params = oineus.ReductionParams()
			params.n_threads, params.kernel, params.image, params.cokernel, params.verbose = read_kernel_image_cokernel(structure_file, kernel_settings)
			spread = math.floor(max(points["z"]))-math.ceil(min(points["z"]))
			sub = sub_complex(points, math.floor(max(points["z"]))-0.05*spread, math.ceil(min(points["z"]))+0.05*spread)
			K, L, L_to_K = oineus_pair(points, sub)
			kicr = oineus.compute_kernel_image_cokernel_reduction(K, L, L_to_K, params)
			if params.kernel:
				pd_1 = kicr.kernel_diagrams().in_dimension(1)
				try:
					pandas.DataFrame(pd_1).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_kernel_PD_1.csv")
					fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_kernel_PD_1.png")
					matplotlib.pyplot.close()
				except:
					print("Kernel diagram in dimension 1 is empty.")
				pd_2 = kicr.kernel_diagrams().in_dimension(2)
				try:
					pandas.DataFrame(pd_2).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_kernel_PD_2.csv")
					fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_kernel_PD_2.png")
					matplotlib.pyplot.close()
				except:
					print("Kernel diagram in dimension 2 is empty.")
				APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
				APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
				pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_kernel_APF_1.csv")
				fig = plot_APF(APF_1, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_kernel_APF_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_kernel_APF_2.csv")
				fig = plot_APF(APF_2, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_kernel_APF_2.png")
				matplotlib.pyplot.close()
			if params.image:
				pd_1 = kicr.image_diagrams().in_dimension(1)
				try:
					pandas.DataFrame(pd_1).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_image_PD_1.csv")
					fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_image_PD_1.png")
					matplotlib.pyplot.close()
				except:
					print("Image diagram in dimension 1 is empty.")
				pd_2 = kicr.kernel_diagrams().in_dimension(2)
				try:
					pandas.DataFrame(pd_2).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_image_PD_2.csv")
					fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_image_PD_2.png")
					matplotlib.pyplot.close()
				except:
					print("Image diagram in dimension 2 is empty.")
				APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
				APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
				pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_image_APF_1.csv")
				fig = plot_APF(APF_1, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_image_APF_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_image_APF_2.csv")
				fig = plot_APF(APF_2, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_image_APF_2.png")
				matplotlib.pyplot.close()
			if params.cokernel:
				pd_1 = kicr.kernel_diagrams().in_dimension(1)
				try:
					pandas.DataFrame(pd_1).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_cokernel_PD_1.csv")
					fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_cokernel_PD_1.png")
					matplotlib.pyplot.close()
				except:
					print("Cokernel diagram in dimension 1 is empty.")
				pd_2 = kicr.kernel_diagrams().in_dimension(2)
				try:
					pandas.DataFrame(pd_2).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_cokernel_PD_2.csv")
					fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
					matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_cokernel_PD_2.png")
					matplotlib.pyplot.close()
				except:
					print("Cokernel diagram in dimension 2 is empty.")
				APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
				APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
				pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_cokernel_APF_1.csv")
				fig = plot_APF(APF_1, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_cokernel_APF_1.png")
				matplotlib.pyplot.close()
				pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_cokernel_APF_2.csv")
				fig = plot_APF(APF_2, 'blue')
				matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_cokernel_APF_2.png")
				matplotlib.pyplot.close()
		
		
def multi_mode(config_file : str, structure_file : str, format : str, structure : str, sample_range : str, save_all : bool = False):
	"""! Multi mode for CLI usage. 
	
	@param config_file	intial configuration file
	@param structure_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param structure 	specific structure to load from the configuration file
	@param sample_range		specific sample settings to use, should be listed in the structure file
	@param save_all		boolean to save images for each sample
 	"""
	atom_list, radii, repeat_x, repeat_y, repeat_z = read_configuration(structure_file, structure)
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []
	sample_settings = read_sample(structure_file, sample_range)
	sample_every = list(range(sample_settings[0], sample_settings[1], sample_settings[2]))
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []	
	atoms = load_atom_file(config_file, format)
	dir = os.path.dirname(config_file)
	if not os.path.exists(os.path.join(dir, "PD1")):
		os.mkdir(os.path.join(dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(dir, "APF2"))
	config_name = os.path.splitext(os.path.split(config_file)[1])[0]
	for s in sample_every:
		print("looking at sample {}".format(s))
		points = sample_at(atoms, s, repeat_x, repeat_y, repeat_z, atom_list, radii)
		dgms = persistent_homology_diagrams_from_points_dionysus(points)
		births_s, deaths_s = get_birth_death(dgms)
		APF_1 = calculate_APF(births_s[1], deaths_s[1])
		APF_2 = calculate_APF(births_s[2], deaths_s[2])
		if save_all:
			pandas.DataFrame(numpy.column_stack([births_s[1], deaths_s[1]])).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.csv")
			fig = plot_PD(births_s[1], deaths_s[1], 'blue')
			matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.png")
			matplotlib.pyplot.close()
			pandas.DataFrame(numpy.column_stack([births_s[2], deaths_s[2]])).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.csv")
			fig = plot_PD(births_s[2], deaths_s[2], 'blue')
			matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(s)+"_PD_2.png")
			matplotlib.pyplot.close()
			pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.csv")
			fig = plot_APF(APF_1, 'blue')
			matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.png")
			matplotlib.pyplot.close()
			pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.csv")
			fig = plot_APF(APF_2, 'blue')
			matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(s)+"_APF_2.png")
			matplotlib.pyplot.close()
		births_1.append(births_s[1])
		deaths_1.append(deaths_s[1])		
		APFs_1.append(APF_1)
		APFs_2.append(APF_2)
	fig = plot_PDs(births_1, deaths_1)
	matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_samples_PD_1.png")
	matplotlib.pyplot.close()
	births_2.append(births_s[2])
	deaths_2.append(deaths_s[2])
	fig = plot_PDs(births_2, deaths_2)
	matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_samples_PD_2.png")
	matplotlib.pyplot.close()
	fig = plot_APFs(APFs_1)
	matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_samples_AFP_1.png")
	matplotlib.pyplot.close()
	fig = plot_APFs(APFs_2)
	matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_samples_AFP_2.png")
	matplotlib.pyplot.close()
	
 
 
def single_mode(config_file : str, format : str, structure_file : str, structure : str, sample_time : int):
	"""! single mode for CLI usage. 
	
	@param config_file	intial configuration file
	@param structure_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param structure 	specific structure to load from the configuration file
	@param sample_time	time at which to sample	
	"""
	print("Hi", structure_file)
	atom_list, radii, repeat_x, repeat_y, repeat_z = read_configuration(structure_file, structure)
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []	
	atoms = load_atom_file(config_file, format)
	dir = os.path.dirname(config_file)
	if not os.path.exists(os.path.join(dir, "PD1")):
		os.mkdir(os.path.join(dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(dir, "APF2"))
	config_name = os.path.splitext(os.path.split(config_file)[1])[0]
	print("looking at sample {}".format(sample_time))
	points = sample_at(atoms, sample_time, repeat_x, repeat_y, repeat_z, atom_list, radii)
	dgms = persistent_homology_diagrams_from_points_dionysus(points)
	births_s, deaths_s = get_birth_death(dgms)
	APF_1 = calculate_APF(births_s[1], deaths_s[1])
	APF_2 = calculate_APF(births_s[2], deaths_s[2])
	pandas.DataFrame(numpy.column_stack([births_s[1], deaths_s[1]])).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_PD_1.csv")
	fig = plot_PD(births_s[1], deaths_s[1], 'blue')
	matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_PD_1.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(numpy.column_stack([births_s[2], deaths_s[2]])).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_PD_2.csv")
	fig = plot_PD(births_s[2], deaths_s[2], 'blue')
	matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_PD_2.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_APF_1.csv")
	fig = plot_APF(APF_1, 'blue')
	matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_APF_1.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_APF_2.csv")
	fig = plot_APF(APF_2, 'blue')
	matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_APF_2.png")
	matplotlib.pyplot.close()
  
def single_mode_kernel(config_file : str, format : str, structure_file : str, structure : str, sample_time : int, kernel_settings : str):
	"""! single mode for CLI usage. 
	
	@param config_file	intial configuration file
	@param structure_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param structure 	specific structure to load from the configuration file
	@param sample_time	time at which to sample	
	@param kernel 		set to true if you want to do any kernel/image/cokernel persistence, (Default = False)
	"""
	import oineus
	atom_list, radii, repeat_x, repeat_y, repeat_z = read_configuration(structure_file, structure)
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []
	births_1 = []
	deaths_1 = []
	births_2 = []
	deaths_2 = []
	APFs_1 = []
	APFs_2 = []	
	atoms = load_atom_file(config_file, format)
	dir = os.path.dirname(config_file)
	if not os.path.exists(os.path.join(dir, "PD1")):
		os.mkdir(os.path.join(dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(dir, "APF2"))
	config_name = os.path.splitext(os.path.split(config_file)[1])[0]
	print("looking at sample {}".format(sample_time))
	points = sample_at(atoms, sample_time, repeat_x, repeat_y, repeat_z, atom_list, radii)
	dgms = persistent_homology_diagrams_from_points_dionysus(points)
	births_s, deaths_s = get_birth_death(dgms)
	APF_1 = calculate_APF(births_s[1], deaths_s[1])
	APF_2 = calculate_APF(births_s[2], deaths_s[2])
	pandas.DataFrame(numpy.column_stack([births_s[1], deaths_s[1]])).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_PD_1.csv")
	fig = plot_PD(births_s[1], deaths_s[1], 'blue')
	matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_PD_1.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(numpy.column_stack([births_s[2], deaths_s[2]])).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_PD_2.csv")
	fig = plot_PD(births_s[2], deaths_s[2], 'blue')
	matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_PD_2.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_APF_1.csv")
	fig = plot_APF(APF_1, 'blue')
	matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_APF_1.png")
	matplotlib.pyplot.close()
	pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_APF_2.csv")
	fig = plot_APF(APF_2, 'blue')
	matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_APF_2.png")
	matplotlib.pyplot.close()
	params = oineus.ReductionParams()
	params.n_threads, params.kernel, params.image, params.cokernel, params.verbose = read_kernel_image_cokernel(structure_file, kernel_settings)
	sub = sub_complex(points, math.floor(max(points["z"])), math.ceil(min(points["z"])))
	K, L, L_to_K = oineus_pair(points, sub)
	kicr = oineus.compute_kernel_image_cokernel_reduction(K, L, L_to_K, params)
	if params.kernel:
		pd_1 = kicr.kernel_diagrams().in_dimension(1)
		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_kernel_PD_1.csv")
		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_kernel_PD_1.png")
		matplotlib.pyplot.close()
		pd_2 = kicr.kernel_diagrams().in_dimension(2)
		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_kernel_PD_2.csv")
		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_kernel_PD_2.png")
		matplotlib.pyplot.close()
		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_kernel_APF_1.csv")
		fig = plot_APF(APF_1, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_kernel_APF_1.png")
		matplotlib.pyplot.close()
		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_kernel_APF_2.csv")
		fig = plot_APF(APF_2, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_kernel_APF_2.png")
		matplotlib.pyplot.close()
	if params.image:
		pd_1 = kicr.image_diagrams().in_dimension(1)
		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_image_PD_1.csv")
		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_image_PD_1.png")
		matplotlib.pyplot.close()
		pd_2 = kicr.kernel_diagrams().in_dimension(2)
		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_image_PD_2.csv")
		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_image_PD_2.png")
		matplotlib.pyplot.close()
		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_image_APF_1.csv")
		fig = plot_APF(APF_1, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_image_APF_1.png")
		matplotlib.pyplot.close()
		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_image_APF_2.csv")
		fig = plot_APF(APF_2, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_image_APF_2.png")
		matplotlib.pyplot.close()
	if params.cokernel:
		pd_1 = kicr.kernel_diagrams().in_dimension(1)
		pandas.DataFrame(numpy.column_stack(pd_1)).to_csv(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.csv")
		fig = plot_PD(pd_1[:,0], pd_1[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD1/"+config_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.png")
		matplotlib.pyplot.close()
		pd_2 = kicr.kernel_diagrams().in_dimension(2)
		pandas.DataFrame(numpy.column_stack(pd_2)).to_csv(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.csv")
		fig = plot_PD(pd_2[:,0], pd_2[:,1], 'blue')
		matplotlib.pyplot.savefig(dir+"/PD2/"+config_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.png")
		matplotlib.pyplot.close()
		APF_1 = calculate_APF(pd_1[:,0], pd_1[:,1])
		APF_2 = calculate_APF(pd_2[:,0], pd_2[:,1])
		pandas.DataFrame(APF_1).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.csv")
		fig = plot_APF(APF_1, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF1/"+config_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.png")
		matplotlib.pyplot.close()
		pandas.DataFrame(APF_2).to_csv(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.csv")
		fig = plot_APF(APF_2, 'blue')
		matplotlib.pyplot.savefig(dir+"/APF2/"+config_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.png")
		matplotlib.pyplot.close()



  
  
