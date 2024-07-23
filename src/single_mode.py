import streamlit as st
import oineus
import numpy as np
import pandas as pd
import os

import streamlit_functions


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
	config_name = st.text_input("Configuration name:", key="comp_config_name")
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
	if st.session_state["manual_config"]:
		load_configuration_file()
	if st.session_state["maual_comp_config"]:
		load_computation_settings()
	st.session_state.atom_locations  = sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, rst.session_state.epeat_z, st.session_state.atoms, st.session_state.radii)
	st.session_state.processed=True
	top_pt = max(st.session_state.atom_locations["z"])
	bot_pt = min(st.session_state.atom_locations["z"])
	height = abs(top_pt - bot_pt)
	if st.session_state.thickness == "Auto":
		ut= top_pt - 0.1*height
		lt = bot_pt + 0.1*height
	else:
		ut = top_pt - thickness*height
	st.session_state.kicr, st.session_state.dgm_1, st.session_state.dgm_2 =  oineus_kernel_image_cokernel(st.session_state.atom_locations, st.session_state.params, st.session_state.thickness)
		
st.button("Process", key="process", on_click=compute)
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
# 					pandas.DataFrame(numpy.column_stack([births[1], deaths[1]]), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_PD_1.csv", header = None)
# 					pandas.DataFrame(APF_1, columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_APF_1.csv", header=  None)
# 				except:
# 					print("Persistence diagram in dimension 1 is empty.")
# 				try:
# 					pandas.DataFrame(numpy.column_stack([births[2], deaths[2]]), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_PD_2.csv", header = None)
# 					pandas.DataFrame(APF_2, columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_APF_2.csv", header = None)
# 				except:
# 					print("Persistence diagram in dimension 2 is empty.")
# 				if values_main["kernel"]:
# 					pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_1.csv", header = None)
# 					pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_kernel_PD_2.csv", header = None)
# 				if values_main["image"]:
# 					pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_image_PD_1.csv", header = None)
# 					pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_image_PD_2.csv", header = None)
# 				if values_main["cokernel"]:
# 					pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_1.csv", header = None)
# 					pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+file_name+"_sample_"+str(s)+"_cokernel_PD_2.csv", header = None)
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