##
# @internal
# @file Single_Configuration.py
# @brief Streamlit page for analysing a single configuration.
# @version 0.6
# @date December 2024
# @author: Yossi Bokor Bleile

import streamlit as st
# import oineus
import numpy as np
import pandas
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

from toma_functions import *

check_params()



#function to display plots
def display_plots():
	"""! Display plots for a single configuration
	@brief Display plots for a single configuration
	
	This function displays persistence diagrams and APF plots for a single configuration.
	It checks if plots have been generated, and if not, calls generate_plots().
	
	For each sample index, it displays:
	- Persistence diagrams in dimensions 0,1,2 if enabled
	- Kernel/image/cokernel persistence diagrams if those computations are enabled
	- APF plots in dimensions 0,1 if enabled
	- Kernel/image/cokernel APF plots if those computations are enabled
	
	The plots are displayed using Streamlit's plotly_chart function.
	If certain plots are not available (e.g. kernel APFs), error messages are printed.
	"""
	if "plots_generated" not in st.session_state or not st.session_state["plots_generated"]:
		generate_plots()
	plot_tab.write(st.session_state.sample_end)
	plot_tab.write(st.session_state.sample_indices)
	for i,s in enumerate(st.session_state.sample_indices):
		print(i)
		if "pd0" in st.session_state and st.session_state["pd0"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_0[i])
			# else: 
			plot_tab.plotly_chart(st.session_state.fig_pds_0[i])
		if "pd1" in st.session_state and st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_1[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_1[i])
		if "pd2" in st.session_state and st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_2[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_2[i])
		if "apf0" in st.session_state and st.session_state["apf0"] == True:
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
		if "apf1" in st.session_state and st.session_state["apf1"] == True:
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
		if "apf2" in st.session_state and st.session_state["apf2"] == True:
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

st.header("Single Configuration Mode")
st.checkbox("Specify directory and file prefix for saving plots", key="custom_save")
if st.session_state.custom_save:
	st.session_state.save_directory = st.text_input("Save directory (this should be plain text containing the path to the directory to save the files):", key="save_directory_input", placeholder="path/to/save/directory")
	st.session_state.save_file_name = st.text_input("File name to use when saving the files (this should be plain text containing the prefix for the files):", key="save_file_name_input", placeholder="save_file_name")
comp_tab, plot_tab, vis_tab = st.tabs(["Computation", "Plots", "Visuatlisation"]) #create tabs for the various parts
st.session_state.mode = "multi"


# else:
# 	st.session_state.save_directory = os.path.dirname(st.session_state.file_path)
# 	st.session_state.save_file_name = os.path.splitext(os.path.split(st.session_state.file_path)[1])[0]

#initialize the reduction parameters


### lets set up the computation tab
file_path = comp_tab.text_input("Initial structure file (this should be plain text containing the path to the initial structure file):",key="file_path", placeholder="path/to/file.xyz") #specify initial structure file 
file_format = comp_tab.text_input("File format (this should be plain text containing the format of the initial structure file):", key="file_format",placeholder="Auto") #specify format of the initial strutcure file
if "processed_file" in st.session_state:
	if st.session_state.file_path != st.session_state["processed_file"]:
		st.session_state.processed = False

comp_tab.header("Configuration settings")
manual_config = comp_tab.checkbox("Manually specify configuration", key="manual_config")#manually set configuration
if not manual_config:
	st.session_state.config_file = comp_tab.text_input("Configuration file (this should be plain text containing the path to the configuration ini file):", key="configuration_file_input", placeholder="path/to/config.ini")
	st.session_state.config_name = comp_tab.text_input("Section name (this should be plain text containing the name of the section in the file):", key="configuration_name", placeholder="config")
else:
	try:
		st.session_state.atoms = [str(a).strip() for a in comp_tab.text_input("Atoms (this should be plain text of the atomic symbols in the structure separated by commas):", key="atoms_input", placeholder="H,C,N,O").split(",")]
	except:
		st.session_state.atoms = ["H","C","N","O"]
	try:
		st.session_state.radii = [float(r) for r in comp_tab.text_input("Radii (this should be plain text of the atomic radii to use for the computation, separated by commas, in the same order as the atoms):", key="radii_input", placeholder="0.3,0.7,1.0,1.2").split(",")]
	except:
		st.session_state.radii = [0.3,0.7,1.0,1.2]


comp_tab.markdown("Computation settings")
manual_compute = comp_tab.checkbox("Manually specify settings for the computations (i.e number of threds, and if you want to compute kernel/image/cokernel)", key="maual_comp_config")
if not manual_compute:
	not_same_config_file = comp_tab.checkbox("The computation settings are in a different configuration file.", key="not_same_config_file")
	if not_same_config_file:
		st.session_state.comp_file = comp_tab.text_input("Configuration file (this should be plain text containing the path to the configuration ini file):", key="comp_config_file", placeholder="path/to/config.ini")
	else:
		st.session_state.comp_file = st.session_state.config_file
	st.session_state.comp_name = comp_tab.text_input("Section name (this should be plain text containing the name of the section in the file):", key="comp_config_name", placeholder="comp")
else:
	try:
		st.session_state.sample_start = int(comp_tab.text_input("Sample start (this should be plain text containing the start index of the samples to process):", key="sample_start_input", placeholder="0"))
	except:
		st.session_state.sample_start = 0
	try:
		st.session_state.sample_end = int(comp_tab.text_input("Sample end (this should be plain text containing the end index of the samples to process):", key="sample_end_input", placeholder="10"))
	except:
		st.session_state.sample_end = 1
	try:
		st.session_state.sample_step = int(comp_tab.text_input("Sample step (this should be plain text containing the step between samples to process):", key="sample_step_input", placeholder="1"))
	except:
		st.session_state.sample_step = 1
	try:
		st.session_state.repeat_x = int(comp_tab.text_input("Repitition in x-axis (this should be plain text containing the number of times to repeat the structure in the x-axis):", key="repeat_x_input", placeholder="1"))
	except:
		st.session_state.repeat_x = 1
	try:
		st.session_state.repeat_y = int(comp_tab.text_input("Repitition in y-axis (this should be plain text containing the number of times to repeat the structure in the y-axis):", key="repeat_y_input", placeholder="1"))
	except:
		st.session_state.repeat_y = 1
	try:
		st.session_state.repeat_z = int(comp_tab.text_input("Repitition in z-axis (this should be plain text containing the number of times to repeat the structure in the z-axis):", key="repeat_z_input", placeholder="1"))
	except:
		st.session_state.repeat_z = 1
	st.session_state.kernel = comp_tab.checkbox("Compute kernel persistence", key="kernel_check")
	st.session_state.image = comp_tab.checkbox("Compute image persistence", key="image_check")
	st.session_state.cokernel = comp_tab.checkbox("Compute cokernel persistence", key="cokernel_check")
	st.session_state.thickness = comp_tab.text_input("Select thickness of top and bottom layer (this should be plain text containing the thickness of the top and bottom layer):", key="thickness_input", placeholder="Automatic detection")
	st.session_state.n_threads = comp_tab.text_input("Select number of threads to use (this should be plain text containing the number of threads to use):", key="n_threads_input", placeholder="4")
	if st.session_state.n_threads == "":
		st.session_state.params.n_threads = 4
		# comp_tab.write("Number of threads is "+str(st.session_state.params.n_threads))
	else:
		st.session_state.params.n_threads = int(st.session_state.n_threads)
	comp_tab.write("Number of threads is "+str(st.session_state.params.n_threads))

	if st.session_state["thickness_input"] == "":
		st.session_state.thickness = "Auto"
	else:
		st.session_state.thickness= float(st.session_state["thickness_input"])

comp_tab.markdown(f"The file selected to analyse is "+file_path)	
settings_into = "<span style='color:#2E86C1'>You have selected to process the file:</span> " + str(st.session_state.file_path) + "<br><span style='color:#28B463'>You have selected to use the `.ini` file:</span> " + str(st.session_state.config_file) + "<br><span style='color:#D35400'>You have selected to use the configuration:</span> " + str(st.session_state.config_name) + "<br><span style='color:#8E44AD'>You have selected to use the computation settings file:</span> " + str(st.session_state.comp_file) + "<br><span style='color:#C0392B'>You have selected to use the computation settings section:</span> " + str(st.session_state.comp_name)
# comp_tab.write("You have selected to use the computation settings file: ", st.session_state.comp_file)
# comp_tab.write("You have selected to use the computation settings section: ", st.session_state.comp_name)
if file_format == "":
	comp_tab.html(settings_into + "<br><span style='color:#2E86C1'>File format will automatically detected.</span>")
else:
	comp_tab.html(settings_into + "<br><span style='color:#2E86C1'>File format is:</span> " + file_format)

if not st.session_state.custom_save:
	st.session_state.save_directory = os.path.dirname(st.session_state.file_path)
	st.session_state.save_file_name = os.path.splitext(os.path.split(st.session_state.file_path)[1])[0]

comp_buttons = comp_tab.columns(2)
with comp_buttons[0]:
	st.button("Process", key="process", on_click=compute)

with comp_buttons[1]:
	st.button("Save .csv files", key="save_csv_files", on_click=save_dgms_as_csv)


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
if "processed" in st.session_state and st.session_state["processed"] and st.session_state["visualisation"]:# and not st.session_state.kicr_params.kernel and not st.session_state.kicr_params.image and not st.session_state.kicr_params.cokernel:
	vis_tab.write(st.session_state.params)
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
	st.session_state.dfVis = generate_visulisation_df(st.session_state.vis_dgm, st.session_state.dcmps[st.session_state.selected_sample_index].r_data, st.session_state.filts[st.session_state.selected_sample_index], st.session_state.atom_locations_list[st.session_state.selected_sample_index], st.session_state.atoms)
	to_display = ["birth", "death", "lifetime"]
	for a in st.session_state.atoms:
		to_display.append(a)
	viz = vis_tab.dataframe(st.session_state.dfVis[to_display], on_select="rerun")
	if viz.selection.rows == []:
		st.session_state.selected_row = None
	else:
		st.session_state.selected_row = viz.selection.rows
	vis_tab.write("Selected Row:")
	vis_tab.write(st.session_state.selected_row)
	if st.session_state.selected_row != None:
		for cycle_id in st.session_state.selected_row:
			vis_tab.plotly_chart(generate_display(st.session_state.atom_locations, st.session_state.dfVis, cycle_id, st.session_state.filts[st.session_state.selected_sample_index], neighbours = st.session_state["neighbours"]))
elif st.session_state.kicr_params.kernel or st.session_state.kicr_params.image or st.session_state.kicr_params.cokernel:
	vis_tab.markdown("Visulation of kernel/image/cokernel persistent homology is not yet available.")
else:
	vis_tab.markdown("Persistent homology has not been computed, so representative cycles are unknow. Please proces the file, and then return to this tab.")

st.button("test", key="test", on_click=test)

if "params" in st.session_state:
	print(st.session_state.params)

