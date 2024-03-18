##
# @internal
# @file cli.py
# contains the functions for the command line interface
from process import *
from plots import *
import os
import matplotlib


#def process_sample(atom_locations, sample_time : int,  repeat_x, repeat_y, repeat_z, params, upper_threshold = "Auto", lower_threshold = "Auto"):


def single_mode(structure_file : str, file_format : str, configuration_file : str, configuration : str, sample_time : int,  n_threads : int, save_plots : bool, kernel : bool, image : bool, cokernel : bool, upper_threshold, lower_threshold):
	"""! Single mode for CLI usage. 
	
	@param structure_file	intial structure file
	@param file_format 		format the structure file is in
	@param configuration_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param configuration 	specific structure to load from the configuration file
	@param sample_time		time at which to sample	
	@param n_threads		number of threads for Oineus
	@param kernel			True to calculate kernel persistence
	@param image			True to calculate image persistence
	@param cokernel			True to calculate cokernel persistence
	@param upper_threshold	the height above which points are in the subcomplex, set to 'Auto' for automatic selection
	@param lower_threshold	the height below which points are in the subcomplex, set to 'Auto' for automatic selection
	"""
	dir = os.path.abspath(structure_file).split("/")[0:-1]
	print(dir)
	if not os.path.exists(os.path.join(dir, "PD1")):
		os.mkdir(os.path.join(*dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(*dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(*dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(*dir, "APF2"))
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(configuration_file, configuration)
	atom_locations = load_atom_file(file_path, file_format)
	points = sample_at(atom_locations, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
	if params.kernel or params.image or params.cokernel:
		if upper_threshold == "Auto":
			upper_threshold = math.floor(max(points["z"]))
		if lower_threshold == "Auto":
			lower_threshold =  math.ceil(min(points["z"]))
		kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
		if params.kernel:
			pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"])
			pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.csv")
			fig = plot_PD(pd_1, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.png")
			matplotlib.pyplot.close()
			pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"])
			pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.csv")
			fig = plot_PD(pd, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.png")
			matplotlib.pyplot.close()
		if params.image:
			pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"])
			pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.csv")
			fig = plot_PD(pd_1, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.png")
			matplotlib.pyplot.close()
			pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"])
			pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.csv")
			fig = plot_PD(pd, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.png")
			matplotlib.pyplot.close()
		if params.cokernel:
			pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"])
			pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.csv")
			fig = plot_PD(pd_1, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.png")
			matplotlib.pyplot.close()
			pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth","death"])
			pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.csv")
			fig = plot_PD(pd, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.png")
			matplotlib.pyplot.close()
			APF = calculate_APF(pd)
			APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.csv")
			fig = plot_APF(APF, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.png")
			matplotlib.pyplot.close()
	else:
		dcmp, filt, dgm_1, dgm_2 = oineus_process(points, params)
		dgm_1.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
		dgm_2.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
		fig = plot_PD(dgm_1, structure_name)
		matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.png")
		matplotlib.pyplot.close()
		fig = plot_PD(dgm_2, structure_name)
		matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_2.png")
		matplotlib.pyplot.close()
		APF_1 = calculate_APF(dgm_1)
		APF_1.to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.csv")
		fig = plot_APF(APF_1, structure_name)
		matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.png")
		matplotlib.pyplot.close()
		APF_2 = calculate_APF(dgm_2)
		APF_2.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.csv")
		fig = plot_APF(APF_2, structure_name)
		matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.png")
		matplotlib.pyplot.close()


def multi_mode(structure_file : str, file_format : str, configuration_file : str, configuration : str, sample_start, sample_end, sample_step,  n_threads : int, kernel : bool, image : bool, cokernel : bool,upper_threshold = "Auto", lower_threshold = "Auto"):
	"""! Multi mode for CLI usage. 
	
	@param structure_file	intial structure file
	@param file_format 		format the structure file is in
	@param configuration_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param configuration 	specific structure to load from the configuration file
	@param sample_start		time at which to start samples
	@param sample_end		time at which to end samples
	@param sample_step		step between samples
	@param n_threads		number of threads for Oineus
	@param kernel			True to calculate kernel persistence
	@param image			True to calculate image persistence
	@param cokernel			True to calculate cokernel persistence
	@param upper_threshold	the height above which points are in the subcomplex, set to 'Auto' for automatic selection
	@param lower_threshold	the height below which points are in the subcomplex, set to 'Auto' for automatic selection
	"""

	
	dir = os.path.dirname(structure_file)
	if not os.path.exists(os.path.join(dir, "PD1")):
		os.mkdir(os.path.join(dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(dir, "APF2"))
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(configuration_file, configuration)
	atom_locations = load_atom_file(file_path, file_format)
	structure_name = os.path.splitext(os.path.split(structure_file)[1])[0]
	for sample_time in range(sample_start, sample_end, sample_step):
		print("looking at sample {}".format(s))
		points = sample_at(atom_locations, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
		if params.kernel or params.image or params.cokernel:
			if upper_threshold == "Auto":
				upper_threshold = math.floor(max(points["z"]))
			if lower_threshold == "Auto":
				lower_threshold =  math.ceil(min(points["z"]))
			kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
			if params.kernel:
				pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"])
				pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.csv")
				fig = plot_PD(pd_1, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.png")
				matplotlib.pyplot.close()
				pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"])
				pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.csv")
				fig = plot_PD(pd, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.png")
				matplotlib.pyplot.close()
			if params.image:
				pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"])
				pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.csv")
				fig = plot_PD(pd_1, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.png")
				matplotlib.pyplot.close()
				pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"])
				pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.csv")
				fig = plot_PD(pd, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.png")
				matplotlib.pyplot.close()
			if params.cokernel:
				pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"])
				pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.csv")
				fig = plot_PD(pd_1, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.png")
				matplotlib.pyplot.close()
				pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth","death"])
				pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.csv")
				fig = plot_PD(pd, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.png")
				matplotlib.pyplot.close()
				APF = calculate_APF(pd)
				APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.csv")
				fig = plot_APF(APF, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.png")
				matplotlib.pyplot.close()
		else:
			dcmp, filt, dgm_1, dgm_2 = oineus_process(points, params)
			dgm_1.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
			dgm_2.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
			fig = plot_PD(dgm_1, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.png")
			matplotlib.pyplot.close()
			fig = plot_PD(dgm_2, structure_name)
			matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_2.png")
			matplotlib.pyplot.close()
			APF_1 = calculate_APF(dgm_1)
			APF_1.to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.csv")
			fig = plot_APF(APF_1, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.png")
			matplotlib.pyplot.close()
			APF_2 = calculate_APF(dgm_2)
			APF_2.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.csv")
			fig = plot_APF(APF_2, structure_name)
			matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.png")
			matplotlib.pyplot.close()


def batch_mode(parent_dir : str, file_ext : str, file_format : str,configuration_file : str, configuration : str, sample_start, sample_end, sample_step,  n_threads : int, kernel : bool, image : bool, cokernel : bool,upper_threshold = "Auto", lower_threshold = "Auto"):
	"""! Batch mode for CLI usage. 
	
	@param parent_dir 	base directory from which to initialise the search
	@param file_ext		extension used for the files containing the initial configurations, files with any other extension will be ignored
	@param configuration_file 	file which contains the structure (and maybe others) you want to use, should follow INI formatting
	@param configuration 	specific structure to load from the configuration file
	@param sample_start		time at which to start samples
	@param sample_end		time at which to end samples
	@param sample_step		step between samples
	@param n_threads		number of threads for Oineus
	@param kernel			True to calculate kernel persistence
	@param image			True to calculate image persistence
	@param cokernel			True to calculate cokernel persistence
	@param upper_threshold	the height above which points are in the subcomplex, set to 'Auto' for automatic selection
	@param lower_threshold	the height below which points are in the subcomplex, set to 'Auto' for automatic selection
	"""
	
	if not os.path.exists(os.path.join(parent_dir, "PD1")):
		os.mkdir(os.path.join(parent_dir, "PD1"))
	if not os.path.exists(os.path.join(dir, "PD2")):
		os.mkdir(os.path.join(parent_dir, "PD2"))
	if not os.path.exists(os.path.join(dir, "APF1")):
		os.mkdir(os.path.join(parent_dir, "APF1"))
	if not os.path.exists(os.path.join(dir, "APF2")):
		os.mkdir(os.path.join(parent_dir, "APF2"))
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(configuration_file, configuration)
	#sample_settings = read_sample(configuration_file, sample_range)
	#sample_every = list(range(sample_settings[0], sample_settings[1], sample_settings[2]))
	structure_files = []
	for root, dirs, files in os.walk(parent_dir):
		for f in files:
			if f.endswith(file_ext):
				structure_files.append(os.path.join(root, f))
	print("Have found the following configuration files:")
	for f in structure_files:
		print(f)
	for structure_file in structure_files:
		print("Looking at {}".format(structure_file))
		atom_locations = load_atom_file(structure_file)
		structure_name = os.path.splitext(os.path.split(structure_file)[1])[0]
		for sample_time in range(sample_start, sample_end, sample_step):
			print("looking at sample {}".format(s))
			points = sample_at(atom_locations, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
			if params.kernel or params.image or params.cokernel:
				if upper_threshold == "Auto":
					upper_threshold = math.floor(max(points["z"]))
				if lower_threshold == "Auto":
					lower_threshold =  math.ceil(min(points["z"]))
				kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
				if params.kernel:
					pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"])
					pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.csv")
					fig = plot_PD(pd_1, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_1.png")
					matplotlib.pyplot.close()
					pd = pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"])
					pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.csv")
					fig = plot_PD(pd, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_APF_2.png")
					matplotlib.pyplot.close()
				if params.image:
					pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"])
					pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.csv")
					fig = plot_PD(pd_1, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_1.png")
					matplotlib.pyplot.close()
					pd = pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"])
					pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.csv")
					fig = plot_PD(pd, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_image_APF_2.png")
					matplotlib.pyplot.close()
				if params.cokernel:
					pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"])
					pd.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.csv")
					fig = plot_PD(pd_1, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					pandas.DataFrame(APF).to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_1.png")
					matplotlib.pyplot.close()
					pd = pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth","death"])
					pd.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.csv")
					fig = plot_PD(pd, structure_name)
					matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2.png")
					matplotlib.pyplot.close()
					APF = calculate_APF(pd)
					APF.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.csv")
					fig = plot_APF(APF, structure_name)
					matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_APF_2.png")
					matplotlib.pyplot.close()
			else:
				dcmp, filt, dgm_1, dgm_2 = oineus_process(points, params)
				dgm_1.to_csv(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
				dgm_2.to_csv(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.csv")
				fig = plot_PD(dgm_1, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1.png")
				matplotlib.pyplot.close()
				fig = plot_PD(dgm_2, structure_name)
				matplotlib.pyplot.savefig(dir+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_2.png")
				matplotlib.pyplot.close()
				APF_1 = calculate_APF(dgm_1)
				APF_1.to_csv(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.csv")
				fig = plot_APF(APF_1, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF1/"+structure_name+"_sample_"+str(sample_time)+"_APF_1.png")
				matplotlib.pyplot.close()
				APF_2 = calculate_APF(dgm_2)
				APF_2.to_csv(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.csv")
				fig = plot_APF(APF_2, structure_name)
				matplotlib.pyplot.savefig(dir+"/APF2/"+structure_name+"_sample_"+str(sample_time)+"_APF_2.png")
				matplotlib.pyplot.close()