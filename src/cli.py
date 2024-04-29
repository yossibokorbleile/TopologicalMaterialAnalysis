##
# @internal
# @file cli.py
# contains the functions for the command line interface
from process import *
from plots import *
import os
import matplotlib.pyplot


#def process_sample(atom_locations, sample_time : int,  repeat_x, repeat_y, repeat_z, params, upper_threshold = "Auto", lower_threshold = "Auto"):


def single_mode(structure_file : str, file_format : str, configuration_file : str, configuration : str, sample_time : int,  n_threads : int, save_plots : bool, kernel : bool, image : bool, cokernel : bool, thickness = "Auto"): #upper_threshold, lower_threshold
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
	dir_path= ""
	for d in os.path.abspath(structure_file).split("/")[0:-1]:
		print(d)
		dir_path+= d+"/"
	print(dir_path)
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
	atom_locations = sample_at(structure_file, file_format, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
	structure_name = os.path.splitext(os.path.split(structure_file)[1])[0]
	if params.kernel or params.image or params.cokernel:
		top_pt = max(atom_locations["z"])
		bot_pt = min(atom_locations["z"])
		height = abs(top_pt - bot_pt)
		if thickness == "Auto":
			ut= top_point - 0.1*height
			lt = bot_pt + 0.1*height
		else:
			ut= top_point - thickness*height
			lt = bot_pt + thickness*height
		# if upper_threshold == "Auto":
		# 	ut = top_pt - 0.01*height
		# else: 
		# 	ut = top_pt - upper_threshold*height
		# if lower_threshold == "Auto":
		# 	lt =  bot_pt + 0.01*height
		# else:
		# 	lt =  bot_pt + lower_threshold*height
		try: 
			kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(atom_locations, params, ut, lt)
		except:
			print("issues with {}".format(structure_file))
			return False
		if params.kernel:
			write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 kernel PD sample "+str(sample_time))
			write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 kernel PD sample "+str(sample_time))
		if params.image:
			write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 image PD sample "+str(sample_time))
			write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 image PD sample "+str(sample_time))
		if params.cokernel:
			write_files(dgm=pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 cokernel PD sample "+str(sample_time))
			write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 cokernel PD sample "+str(sample_time))
		write_files(dgm=dgm_1, file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
		write_files(dgm=dgm_2, file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
	else:
		dcmp, filt, dgm_1, dgm_2 = oineus_process(points, params)
		write_files(dgm=dgm_1, file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
		write_files(dgm=dgm_2, file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
	return True

def multi_mode(structure_file : str, file_format : str, configuration_file : str, configuration : str, sample_start, sample_end, sample_step,  n_threads : int, kernel : bool, image : bool, cokernel : bool,  thickness = "Auto", save_plots = False): # upper_threshold = "Auto", lower_threshold = "Auto",
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
	@param upper_threshold	the height above which atom_locations are in the subcomplex, set to 'Auto' for automatic selection
	@param lower_threshold	the height below which atom_locations are in the subcomplex, set to 'Auto' for automatic selection
	"""

	
	dir_path= os.path.dirname(structure_file)
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
	structure_name = os.path.splitext(os.path.split(structure_file)[1])[0]
	for sample_time in range(sample_start, sample_end, sample_step):
		print("looking at sample {}".format(sample_time))
		try: 
				atom_locations = sample_at(structure_file, file_format, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
		except:
			print("issues with loading {}".format(structure_file))
			return False
		if params.kernel or params.image or params.cokernel:
			top_pt = max(atom_locations["z"])
			bot_pt = min(atom_locations["z"])
			height = abs(top_pt - bot_pt)
			if thickness == "Auto":
				ut= top_point - 0.1*height
				lt = bot_pt + 0.1*height
			else:
				ut= top_point - thickness*height
				lt = bot_pt + thickness*height
			# if upper_threshold == "Auto":
			# 	ut = top_pt - 0.01*height
			# else: 
			# 	ut = top_pt - upper_threshold*height
			# if lower_threshold == "Auto":
			# 	lt =  bot_pt + 0.01*height
			# else:
			# 	lt =  bot_pt + lower_threshold*height
			try: 
				kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(atom_locations, params, ut, lt)
			except:
				print("issues with {}".format(structure_file))
				return False
			if params.kernel:
				write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 1 kernel PD sample "+str(sample_time))
				write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 2 kernel PD sample "+str(sample_time))
			if params.image:
				write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 1 image PD sample "+str(sample_time))
				write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 2 image PD sample "+str(sample_time))
			if params.cokernel:
				write_files(dgm=pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 1 cokernel PD sample "+str(sample_time))
				write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2_upper_"+str(upper_threshold)+"_lower_"+str(lower_threshold), save_plots=save_plots, plot_name=structure_name+" dimension 2 cokernel PD sample "+str(sample_time))
			write_files(dgm=dgm_1, file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
			write_files(dgm=dgm_2, file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_2", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
		else:
			dcmp, filt, dgm_1, dgm_2 = oineus_process(points, params)
			write_files(dgm=dgm_1, file_path=dir_path+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
			write_files(dgm=dgm_2, file_path=dir_path+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_2", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
	return True

def batch_mode(parent_dir : str, file_ext : str, file_format : str,configuration_file : str, configuration : str, sample_start, sample_end, sample_step,  n_threads : int, kernel : bool, image : bool, cokernel : bool, thickness = "Auto", save_plots = False):#upper_threshold = "Auto", lower_threshold = "Auto"
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
	print(parent_dir)
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
		structure_name = os.path.splitext(os.path.split(structure_file)[1])[0]
		structure_folder = os.path.split(structure_file)[0]
		if not os.path.exists(os.path.join(structure_folder, "PD1")):
			os.mkdir(os.path.join(structure_folder, "PD1"))
		if not os.path.exists(os.path.join(structure_folder, "PD2")):
			os.mkdir(os.path.join(structure_folder, "PD2"))
		if not os.path.exists(os.path.join(structure_folder, "APF1")):
			os.mkdir(os.path.join(structure_folder, "APF1"))
		if not os.path.exists(os.path.join(structure_folder, "APF2")):
			os.mkdir(os.path.join(structure_folder, "APF2"))
		
		for sample_time in range(sample_start, sample_end, sample_step):
			print("looking at sample {}".format(sample_time))
			# try: 
			atom_locations = sample_at(structure_file, file_format, sample_time, repeat_x, repeat_y, repeat_z, atoms, radii)
			print("succesfully got atom locations")
			if params.kernel or params.image or params.cokernel:
				top_pt = max(atom_locations["z"])
				bot_pt = min(atom_locations["z"])
				height = abs(top_pt - bot_pt)
				if thickness == "Auto":
					ut= top_point - 0.1*height
					lt = bot_pt + 0.1*height
				else:
					ut = top_point - thickness*height
					lt = bot_pt + thickness*height
				# if upper_threshold == "Auto":
				# 	ut = top_pt - 0.01*height
				# else: 
				# 	ut = top_pt - upper_threshold*height
				# if lower_threshold == "Auto":
				# 	lt =  bot_pt + 0.01*height
				# else:
				# 	lt =  bot_pt + lower_threshold*height
				kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(atom_locations, params, ut, lt)
				print("got kicr for sample {} of structure {}".format(sample, structure_file))
				if params.kernel:
					write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=structure_folder+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 kernel PD sample "+str(sample_time))
					write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=structure_folder+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_kernel_PD_2", save_plots=save_plots, plot_name=structure_name+" dimension 2 kernel PD sample "+str(sample_time))
				if params.image:
					write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth","death"]), file_path=structure_folder+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 image PD sample "+str(sample_time))
					write_files(dgm=pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth","death"]), file_path=structure_folder+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_image_PD_2", save_plots=save_plots, plot_name=structure_name+" dimension 2 image PD sample "+str(sample_time))
				if params.cokernel:
					write_files(dgm=pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth","death"]), file_path=structure_folder+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 cokernel PD sample "+str(sample_time))
					write_files(dgm=pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth","death"]), file_path=structure_folder+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_cokernel_PD_2", save_plots=save_plots, plot_name=structure_name+" dimension 2 cokernel PD sample "+str(sample_time))
				write_files(dgm=dgm_1, file_path=structure_folder+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
				write_files(dgm=dgm_2, file_path=structure_folder+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
			else:
				dcmp, filt, dgm_1, dgm_2 = oineus_process(atom_locations, params)
				write_files(dgm=dgm_1, file_path=structure_folder+"/PD1/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 1 PD sample "+str(sample_time))
				write_files(dgm=dgm_2, file_path=structure_folder+"/PD2/"+structure_name+"_sample_"+str(sample_time)+"_PD_1", save_plots=save_plots, plot_name=structure_name+" dimension 2 PD sample "+str(sample_time))
			# except:
			# 	print("issues with loading {}".format(structure_file))
	return True