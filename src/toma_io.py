##
# @internal
# @file io.py
# @brief functions to read configuration and computation settings, load atoms as well as save outputs.
# @version 0.1
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

def load_configuration_settings():
	"""! Load the configuration settings
	@brief Load configuration settings from the session state
	
	This function reads configuration settings from st.session_state and loads them into the session state.
	It uses read_configuration() to parse the config file and extract:
	- atoms: ASE Atoms object containing atomic positions and types
	- radii: Dictionary mapping atom types to their radii
	- repeat_x/y/z: Number of times to repeat the cell in each direction
	
	The settings are read from:
	- st.session_state["config_file"]: Path to configuration file
	- st.session_state["config_name"]: Name of configuration section to use
	
	The parsed settings are stored back in st.session_state as:
	- st.session_state.atoms
	- st.session_state.radii  
	- st.session_state.repeat_x/y/z
	"""	
	st.session_state.atoms, st.session_state.radii, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z = read_configuration(st.session_state["config_file"], st.session_state["config_name"])

def load_computation_settings():
	"""! Load the computation settings
	@brief Load computation settings from the session state
	
	This function reads computation settings from st.session_state and loads them into the session state.
	It uses read_computation_settings() to parse the settings file and extract:
	- n_threads: Number of threads for Oineus to use
	- save_plots: Whether to save plots
	- kernel: Whether to compute kernel persistence
	- image: Whether to compute image persistence
	- cokernel: Whether to compute cokernel persistence
	- thickness: Thickness settings
	"""
	st.session_state.n_threads, st.session_state.save_plots, st.session_state.kernel, image, st.session_state.cokernel, st.session_state.thickness = process.load_computation_settings(st.session_state["comp_file"], st.session_state["comp_name"])


def write_dgm_csv(dgm, file_path,  plot_name = ""):
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
		# pd.DataFrame(APF, columns=["mean age", "lifetime"]).to_csv(file_path.replace("PD","APF")+".csv")
	return True

# function to save plots
def save_plots():
	"""! Save plots to a directory
	@brief Save plots to a directory
	
	This function saves plots to a specified directory.
	It uses the following from st.session_state:
	- st.session_state.file_path: Path to the file containing the data
	- st.session_state.sample_indices: Indices of samples to process
	- st.session_state.params: Computation settings
	"""	
	dir_name = os.path.dirname(st.session_state.file_path)
	file_name = os.path.splitext(os.path.split(st.session_state.file_path)[1])[0]
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

