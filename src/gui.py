##
# @internal
# @file gui.py
# @brief Create the graphical user interfaces. 
# There are 3 options:
# -# single_mode
# -# multi_mode
# -# batch_mode

from process import *
from plots import *
from visualisation import *
import os
import oineus

import PySimpleGUI as sg	

#sg.Print('Re-routing the stdout', do_not_reroute_stdout=False)
#sg.Print('Re-routing the stdout', do_not_reroute_stderr=False)

# this is clobbering the print command, and replacing it with sg's Print()
#print = sg.Print

# this will now output to the sg display.
#print('This is a normal print that has been re-routed.')


TitleFont = 'Times 28'
HeaderFont1 = 'Times 22'
HeaderFont2 = 'Times 20'
TextFont = 'Times 16'
InputFont = 'Times 14'
ButtonFont = 'Times 16'
sg.theme('LightGrey')
menu_bar = [["AMA", ["Quit"]], ["Mode", ["Single", "Multi", "Batch"]]]


def entry_window():
	"""! Entry point for the GUI
	"""
	right_click_menu_entry = ["Unused", ["Single", "Multi", "Batch", "Quit"]]
	
	
	layout_entry = [[sg.Menu(menu_bar, key="menu")], [sg.Button("Single", font=ButtonFont), sg.Button("Multi", font=ButtonFont), sg.Button("Batch", font=ButtonFont), sg.Button("Quit", font=ButtonFont)]]
	window_entry = sg.Window("Topological Amorphous Material Analysis", layout_entry, resizable = True, size=(500,500), right_click_menu=right_click_menu_entry)
	while True:
		event_entry, values_entry = window_entry.read()
		if event_entry == "Single":
			window_entry.close()
			single_mode()
		if event_entry == "Multi":
			window_entry.close()
			multi_mode()
		if event_entry == "Batch":
			window_entry.close()
			batch_mode()
		if event_entry == "Quit" or event_entry == sg.WIN_CLOSED:
			break
	window_entry.close()
 
def single_mode():
	"""! GUI for single mode
 
	This mode is designed to investigate the properties of a single structure at a specific time sample. 
  	It should not be used for generating persistence diagrams and accumulated persistence functions plots for several structures or samples. 
	For several samples of the same structure, you can use multi_mode, and for various structures you should batch mode.
  	"""
	right_click_menu = ['Unused', ["Multi", "Batch", "Exit", "Quit"]]
	#menu_def = [["AMA", ["Quit"]], ["Mode", ["Multi", "Batch"]]]
	layout_main = [[sg.Menu(menu_bar, key="menu")],[sg.Text("Single Mode", font=TitleFont)],  [sg.Text("Intial configuration file:", font=TextFont), sg.FileBrowse(key="file_path", font=ButtonFont), sg.Text("Enter file format:", font=TextFont), sg.Input("xyz", key="file_format", font=InputFont)], [sg.Text("Configuration file settings", font=HeaderFont1)], [sg.Text("File:", font=TextFont), sg.FileBrowse(key="config_file", font=ButtonFont), sg.Text("Structure name:", font=TextFont), sg.Input(key="structure", font=InputFont)], [sg.Text("Manual structure configuration", font=HeaderFont1)], [sg.Text("Atoms:", font=TextFont), sg.Input("H, C, N, Zn", key="atoms", font=InputFont)], [sg.Text("Radii:", font=TextFont), sg.Input("0.389, 0.718, 0.635, 1.491", key="radii", font=InputFont)], [sg.Text("Please select index of the sample you want:",font=TextFont), sg.Input("0", key="sample_at")], [sg.Text("Repeating the configuration:", font=HeaderFont2)], [sg.Text("Repeation in x-axis:", font=TextFont), sg.Input("1", key="repeat_x", font=InputFont)], [sg.Text("Repeation in y-axis:", font=TextFont), sg.Input("1", key="repeat_y", font=InputFont)], [sg.Text("Repeation in z-axis:", font=TextFont), sg.Input("1", key="repeat_z", font=InputFont)], [sg.Text("Kernel/Image/Cokernel Settings:", font=HeaderFont2)], [sg.Text("Number of threads:", font=TextFont), sg.Input("4", font=InputFont, key="n_threads")], [sg.Text("Upper Threshold:", font=TextFont), sg.Input("Automatic", key="upper_threshold")], [sg.Text("Lower Threshold:", font=TextFont), sg.Input("Automatic", key="lower_threshold")], [sg.Checkbox("Kernel", key="kernel", font=ButtonFont), sg.Checkbox("Image", key="image", font=ButtonFont), sg.Checkbox("Cokernel", key="cokernel", font=ButtonFont)], [sg.Text("Please select which plots you would like to generate:", font=HeaderFont1)],[sg.Text("Persistence Diagram", font=TextFont), sg.Checkbox("Dimension 1", key="PD1", font=ButtonFont), sg.Checkbox("Dimension 2", key="PD2", font=ButtonFont)], [sg.Text("Accumulated Persistence Function", font=TextFont), sg.Checkbox("Dimension 1", key="APF1", font=ButtonFont), sg.Checkbox("Dimension 2", key="APF2", font=ButtonFont)], [sg.Button("Process", font=ButtonFont),  sg.Button("Plot", font=ButtonFont), sg.Button("Visualisation",font=ButtonFont),sg.Button("Save", font=ButtonFont), sg.Button("Exit", font=ButtonFont)]]
	
	window_main = sg.Window("Topological Amorphous Material Analysis", layout_main, resizable = True, size=(700,700), right_click_menu=right_click_menu)

	loaded = False
	processed = False
	plotted = False
	params = oineus.ReductionParams()

	while True:
		event_main, values_main = window_main.read()
		if event_main == "Process":
			print("starting")
			if values_main["config_file"] != "":
				atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(values_main["config_file"], values_main["structure"])
			else:	
				atoms = [str(a).strip() for a in values_main["atoms"].split(",")]
				radii = [float(r) for r in values_main['radii'].split(",")]
				repeat_x = int(values_main["repeat_x"])
				repeat_y = int(values_main["repeat_y"])
				repeat_z = int(values_main["repeat_z"])
			file_path = values_main["file_path"]
			points = load_atom_file(file_path, values_main["file_format"])
			print("loaded")
			s = int(values_main["sample_at"])
			points = sample_at(points, s, repeat_x, repeat_y, repeat_z, atoms, radii)
			simplices = weighted_alpha_diode(points)
			params.n_threads = int(values_main["n_threads"])
			if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
				if values_main["upper_threshold"] == "Automatic":
					upper_threshold = math.floor(max(points["z"])) 
				else:
					upper_threshold = float(values_main["upper_threshold"])
				if values_main["lower_threshold"] == "Automatic":
						lower_threshold = math.ceil(min(points["z"]))
				else:
					lower_threshold = float(values_main["lower_threshold"])
				params.kernel = values_main["kernel"]
				params.image = values_main["image"]
				params.cokernel = values_main["cokernel"]
				kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
				if values_main["APF1"]:
					APF_1 = calculate_APF(dgm_1)
				if values_main["APF2"]:
					APF_2 = calculate_APF(dgm_2)	
			else:
				dgm_1, dgm_2 = oineus_process(points, params)
				if values_main["APF1"]:
					APF_1 = calculate_APF(dgm_1)
				if values_main["APF2"]:
					APF_2 = calculate_APF(dgm_2)
			processed = True
   
		if event_main == "Plot":
			if processed == False:
				sg.popup_error("File not processed, please Process it.")
			else:
				if values_main['PD1'] == True:
					fig_pd_1 = plot_PD(dgm_1, file_path+" PD1 sample "+str(s))
					fig_pd_1.show()
					if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
						try:
							fig_kic_pd_1 = plot_kernel_image_cokernel_PD(kicr, 1, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" kernel/image/cokernel dimension 1 sample "+str(s))
							fig_kic_pd_1.show()
						except:
							print("Kernel/image/cokernel persistence has not been calculated")
				if values_main['PD2'] == True:
					fig_pd_2 = plot_PD(dgm_2, file_path+" PD2 sample "+str(s))
					fig_pd_2.show()
					if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
						try:
							fig_kic_pd_2 = plot_kernel_image_cokernel_PD(kicr, 2, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" kernel/image/cokernel dimension 2 sample "+str(s))
							fig_kic_pd_2.show()
						except:
							print("Kernel/image/cokernel persistence has not been calculated")
				if values_main['APF1'] == True:
					try:
						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 sample "+str(s))
						fig_apf_1.show()
					except:
						sg.popup_error("APF1 has not been calculated, will calculate it now. Did you change settins after processig the structure?")
						APF_1 = calculate_APF(dgm_1)
						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 sample "+str(s))
						fig_apf_1.show()
				if values_main['APF2'] == True:
					try:
						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 sample "+str(s))
						fig_apf_2.show()	
					except:
						sg.popup_error("APF2 has not been calculated, will calculate it now. Did you change settins after processig the structure?")
						APF_2 = calculate_APF(dgm_2)
						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 sample "+str(s))
						fig_apf_2.show()
				plotted = True
				
		
		if event_main == "Save":
			if processed == False:
				sg.popup_error("File not processed, please process it.")
			else:
				dir = os.path.dirname(file_path)
				config_name = os.path.splitext(os.path.split(file_path)[1])[0]
				try:
					pandas.DataFrame(numpy.column_stack([births[1], deaths[1]]), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_PD_1.csv", header = None)
					pandas.DataFrame(APF_1, columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_APF_1.csv", header=  None)
				except:
					print("Persistence diagram in dimension 1 is empty.")
				try:
					pandas.DataFrame(numpy.column_stack([births[2], deaths[2]]), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_PD_2.csv", header = None)
					pandas.DataFrame(APF_2, columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_APF_2.csv", header = None)
				except:
					print("Persistence diagram in dimension 2 is empty.")
				if values_main["kernel"]:
					pandas.DataFrame(kicr.kernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_kernel_PD_1.csv", header = None)
					pandas.DataFrame(kicr.kernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_kernel_PD_2.csv", header = None)
				if values_main["image"]:
					pandas.DataFrame(kicr.image_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_image_PD_1.csv", header = None)
					pandas.DataFrame(kicr.image_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_image_PD_2.csv", header = None)
				if values_main["cokernel"]:
					pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(1), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_cokernel_PD_1.csv", header = None)
					pandas.DataFrame(kicr.cokernel_diagrams().in_dimension(2), columns=["birth", "death"]).to_csv(dir+"/"+config_name+"_sample_"+str(s)+"_cokernel_PD_2.csv", header = None)
				if plotted == False:
					if values_main['APF1'] == True:
						fig_apf_1 = plot_APF(APF_1, file_path+" APF1 "+str(s))
						fig_apf_1.show()
					if values_main['APF2'] == True:
						fig_apf_2 = plot_APF(APF_2, file_path+" APF2 "+str(s))
						fig_apf_2.show()			
					if values_main['PD1'] == True:
						fig_pd_1 = plot_PD(births[1], deaths[1], file_path+" PD1 "+str(s))
						fig_pd_1.show()
						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
							fig_kic_pd_1 = plot_kernel_image_cokernel_PD(kicr, 1, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" kernel/image/cokernel dimension 1"+str(s))
							fig_kic_pd_1.show()
					if values_main['PD2'] == True:
						fig_pd_2 = plot_PD(births[2], deaths[2], file_path+" PD2 "+str(s))
						fig_pd_2.show()
						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
							fig_kic_pd_2 = plot_kernel_image_cokernel_PD(kicr, 2, values_main["kernel"], values_main["image"], values_main["cokernel"], file_path+" kernel/image/cokernel dimension 2"+str(s))
							fig_kic_pd_2.show()
				else:
					if values_main["PD1"]:
						fig_pd_1.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_PD_1.png")
						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
							fig_kic_pd_1.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_kic_PD_1.png")
					if values_main["PD2"]:
						fig_pd_2.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_PD_2.png")
						if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
							fig_kic_pd_2.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_kic_PD_2.png")
					if values_main["APF1"]:
						fig_apf_1.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_APF_1.png")
					if values_main["APF2"]:
						fig_apf_2.write_image(dir+"/"+config_name+"_sample_"+str(s)+"_APF_2.png")
		
		if event_main == "Visualisation":
			if processed == False:
				sg.popup_error("File not processed, please process it.")
			else:
				dfPD = get_representative_loops(points, atoms, filt, m, dionysus_diagrams)
				visualisation_table_layout = [[sg.Table(values=dfPD.values.tolist(), headings=dfPD.columns.values.tolist(), auto_size_columns=True, num_rows = 50, display_row_numbers=True, selected_row_colors="red on yellow", enable_events=True)], [sg.Button("Display selected", key="Display")]]
				visualisation_table_window = sg.Window("AMA: 1-cycle representatives table", visualisation_table_layout, resizable=True)
				while True:
					event_visualisation, value_visualisation = visualisation_table_window.read()
					if event_visualisation == "Display":
						vis = generate_display(points, dfPD, value_visualisation[0][0], filt)
						vis.show()
					if event_visualisation == "Exit" or event_visualisation == sg.WIN_CLOSED:
						break
		if event_main == "Exit" or event_main == sg.WIN_CLOSED:
			window_main.close()
			entry_window()
			break
		
		if event_main == "Quit":
			break
	window_main.close()
 
def multi_mode():
	"""! GUI for multi mode"""
	right_click_menu = ['Unused', ["Single", "Batch", "Exit", "Quit"]]
	#menu_def = [["AMA", ["Quit"]], ["Mode", ["Single", "Batch"]]]
	layout_main = [[sg.Menu(menu_bar, key="menu")],[sg.Text("Multi Mode", font=TitleFont)],  [sg.Text("Intial configuration file:", font=TextFont), sg.FileBrowse(key="file_path", font=ButtonFont), sg.Text("Enter file format:", font=TextFont), sg.Input("xyz", key="file_format", font=InputFont)], [sg.Text("Configuration file settings", font=HeaderFont1)], [sg.Text("File:", font=TextFont), sg.FileBrowse(key="config_file", font=ButtonFont), sg.Text("Structure name:", font=TextFont), sg.Input(key="structure", font=InputFont)], [sg.Text("Manual structure configuration", font=HeaderFont1)], [sg.Text("Atoms:", font=TextFont), sg.Input("H, C, N, Zn", key="atoms", font=InputFont)], [sg.Text("Radii:", font=TextFont), sg.Input("0.389, 0.718, 0.635, 1.491", key="radii", font=InputFont)], [sg.Text("Enter the start of the sample range:", font=ButtonFont), sg.Input("4005", key="range-start", font=InputFont)], [sg.Text("Enter the end of the sample range:", font=ButtonFont), sg.Input("4405",key="range-end", font=InputFont)], [sg.Text("Enter the sample step range:", font=ButtonFont), sg.Input("200",key="range-step", font=InputFont)], [sg.Text("Repeating the configuration:", font=HeaderFont2)], [sg.Text("Repeation in x-axis:", font=TextFont), sg.Input("1", key="repeat_x", font=InputFont)], [sg.Text("Repeation in y-axis:", font=TextFont), sg.Input("1", key="repeat_y", font=InputFont)], [sg.Text("Repeation in z-axis:", font=TextFont), sg.Input("1", key="repeat_z", font=InputFont)], [sg.Text("Kernel/Image/Cokernel Settings:", font=HeaderFont2)], [sg.Text("Number of threads:", font=TextFont), sg.Input("4", font=InputFont, key="n_threads")], [sg.Text("Upper Threshold:", font=TextFont), sg.Input("Automatic", key="upper_threshold")], [sg.Text("Lower Threshold:", font=TextFont), sg.Input("Automatic", key="lower_threshold")], [sg.Checkbox("Kernel", key="kernel", font=ButtonFont), sg.Checkbox("Image", key="image", font=ButtonFont), sg.Checkbox("Cokernel", key="cokernel", font=ButtonFont)], [sg.Text("Please select which plots you would like to generate:", font=HeaderFont1)],[sg.Text("Persistence Diagram", font=TextFont), sg.Checkbox("Dimension 1", key="PD1", font=ButtonFont), sg.Checkbox("Dimension 2", key="PD2", font=ButtonFont)], [sg.Text("Accumulated Persistence Function", font=TextFont), sg.Checkbox("Dimension 1", key="APF1", font=ButtonFont), sg.Checkbox("Dimension 2", key="APF2", font=ButtonFont)], [sg.Text("Status messages and progress:", key="label", enable_events=True, font=TextFont, justification='left', expand_x=True), sg.MLine(key='messages'+sg.WRITE_ONLY_KEY, size=(200,8))], [sg.Button("Process", font=ButtonFont),  sg.Button("Plot", font=ButtonFont), sg.Button("Visualisation",font=ButtonFont),sg.Button("Save", font=ButtonFont), sg.Button("Exit", font=ButtonFont), sg.Button("Quit", font=ButtonFont)]] 
	#[sg.Text("Progress Tracking", key='-OUT-', enable_events=True, font=TextFont, justification='center', expand_x=True),sg.ProgressBar(100, orientation='h', expand_x=True, size=(20, 20),  key='progress')],

	window_main = sg.Window("Topological Amorphous Material Analysis", layout_main, resizable = True, size=(800,1000), right_click_menu=right_click_menu)
	params = oineus.ReductionParams()
	loaded = False
	processed = False
	while True:
		event_main, values_main = window_main.read()
		if (event_main == "Process"):
			window_main['messages'+sg.WRITE_ONLY_KEY].print("EVENT: Process")
			window_main['messages'+sg.WRITE_ONLY_KEY].print("begun")
			if values_main["config_file"] != "":
				atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(values_main["config_file"], values_main["structure"])
			else:	
				atoms = [str(a).strip() for a in values_main['atoms'].split(",")]
				radii = [float(r) for r in values_main['radii'].split(",")]
				repeat_x = int(values_main["repeat_x"])
				repeat_y = int(values_main["repeat_y"])
				repeat_z = int(values_main["repeat_z"])
			window_main['messages'+sg.WRITE_ONLY_KEY].print("configuration loaded")
			file_path = values_main['file_path']
			dgms_1 = []
			dgms_2 = []
			if values_main["APF1"]:
				APFs_1 = []
			if values_main["APF2"]:
				APFs_2 = []
			window_main['messages'+sg.WRITE_ONLY_KEY].print("settings loaded")
			coords = load_atom_file(file_path, values_main["file_format"])
			window_main['messages'+sg.WRITE_ONLY_KEY].print("coordinates loaded")
			sample_every = list(range(int(values_main["range-start"]), int(values_main["range-end"]), int(values_main["range-step"])))
			params.n_threads = int(values_main["n_threads"])
			for i,s in enumerate(sample_every):
				window_main['messages'+sg.WRITE_ONLY_KEY].print("sample "+str(i)+" out of "+str(len(sample_every))+"\r")
				points = sample_at(coords, s, repeat_x, repeat_y, repeat_z, atoms, radii)
				simplices = weighted_alpha_diode(points)
				#window_main['progress'].update_bar(math.ceil(100*i/len(sample_every)))
				window_main['messages'+sg.WRITE_ONLY_KEY].print("looking at sample "+str(s))
				if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
					if values_main["upper_threshold"] == "Automatic":
						upper_threshold = math.floor(max(points["z"])) 
					else:
						upper_threshold = float(values_main["upper_threshold"])
					if values_main["lower_threshold"] == "Automatic":
							lower_threshold = math.ceil(min(points["z"]))
					else:
						lower_threshold = float(values_main["lower_threshold"])
					params.kernel = values_main["kernel"]
					params.image = values_main["image"]
					params.cokernel = values_main["cokernel"]
					kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
					if values_main["APF1"]:
						APF_1 = calculate_APF(dgm_1)
					if values_main["APF2"]:
						APF_2 = calculate_APF(dgm_2)	
				else:
					dgm_1, dgm_2 = oineus_process(points, params)
					if values_main["APF1"]:
						APF_1 = calculate_APF(dgm_1)
					if values_main["APF2"]:
						APF_2 = calculate_APF(dgm_2)
			#window_main['progress'].update_bar(100)
			window_main['messages'+sg.WRITE_ONLY_KEY].print("finished")
			processed = True



		if event_main == "Plot":
			if processed == False:
				sg.popup_error("File not processed, please Process.")
			else:
				if values_main['APF1'] == True:
					fig_APFs_1 = plot_APFs(APFs_1)
				if values_main['APF2'] == True:
					fig_APFs_2 = plot_APFs(APFs_2)
				if values_main['PD1'] == True:
					fig_PDs_1 = plot_PDs(dgms_1)
				if values_main['PD2'] == True:	
					fig_PDs_2 = plot_PDs(dgms_2)
		
		if event_main == "Save":
			save_path = sg.popup_get_text("Please enter the path to the directory in which you want to save the PDs and APFs.")
			for i, s in enumerate(sample_every):
				pandas.DataFrame.to_csv(save_path+"_sample_"+str(s)+"_PD_1.csv", pandas.DataFrame(numpy.column_stack(births_1[i], deaths_1[i])))
				pandas.DataFrame.to_csv(save_path+"_sample_"+str(s)+"_PD_2.csv", pandas.DataFrame(numpy.column_stack(births_2[i], deaths_2[i])))
				pandas.DataFrame.to_csv(save_path+"_sample_"+str(s)+"_APF_1.csv", pandas.DataFrame(APFs_1[i]))
				pandas.DataFrame.to_csv(save_path+"_sample_"+str(s)+"_APF_2.csv", pandas.DataFrame(APFs_2[i]))
		
		if event_main == "Exit" or event_main == sg.WIN_CLOSED:
			window_main.close()
			entry_window()
			break
		if event_main == "Quit":
			break

	window_main.close()
 
def batch_mode():
	"""! GUI for batch mode"""
	right_click_menu = ['Unused', ["Single", "Multi", "Exit", "Quit"]]
	#menu_def = [["AMA", ["Quit"]], ["Mode", ["Single", "Multi"]]]
	layout_main = [[sg.Menu(menu_bar, key="menu")],[sg.Text("Batch Mode", font=TitleFont)], [ sg.Text("Select parent directory:", font=ButtonFont), sg.FolderBrowse(key="file_dir", font=ButtonFont), sg.Text("Enter file extension:", font=ButtonFont), sg.Input(".xyz", key="file_ext", font=ButtonFont)], [sg.Text("Enter file format:", font=TextFont), sg.Input("xyz", key="file_format", font=InputFont)], [sg.Text("Configuration file settings", font=HeaderFont1)], [sg.Text("File:", font=TextFont), sg.FileBrowse(key="config_file", font=ButtonFont), sg.Text("Structure name:", font=TextFont), sg.Input(key="structure", font=InputFont)],[sg.Text("Manual structure configuration", font=HeaderFont1)], [sg.Text("Atoms:", font=TextFont), sg.Input("H, C, N, Zn", key="atoms", font=InputFont)], [sg.Text("Radii:", font=TextFont), sg.Input("0.389, 0.718, 0.635, 1.491", key="radii", font=InputFont)], [sg.Text("Sample range:", font=HeaderFont2)],[sg.Text("Sample range start:", font=TextFont), sg.Input("4005", key="range-start", font=InputFont)], [sg.Text("Sample range end:", font=TextFont), sg.Input("4405",key="range-end", font=InputFont)], [sg.Text("Sample range step:", font=TextFont), sg.Input("200",key="range-step", font=InputFont)], [sg.Text("Repeating the configuration:", font=HeaderFont2)], [sg.Text("Repeation in x-axis:", font=TextFont), sg.Input("1", key="repeat_x", font=InputFont)], [sg.Text("Repeation in y-axis:", font=TextFont), sg.Input("1", key="repeat_y", font=InputFont)], [sg.Text("Repeation in z-axis:", font=TextFont), sg.Input("1", key="repeat_z", font=InputFont)], [sg.Text("Kernel/Image/Cokernel Settings:", font=HeaderFont2)], [sg.Text("Number of threads:", font=TextFont), sg.Input("4", font=InputFont, key="n_threads")], [sg.Text("Upper Threshold:", font=TextFont), sg.Input("Automatic", key="upper_threshold")], [sg.Text("Lower Threshold:", font=TextFont), sg.Input("Automatic", key="lower_threshold")], [sg.Checkbox("Kernel", key="kernel", font=ButtonFont), sg.Checkbox("Image", key="image", font=ButtonFont), sg.Checkbox("Cokernel", key="cokernel", font=ButtonFont)], [sg.Text("Please select which plots you would like to generate:", font=HeaderFont1)],[sg.Text("Persistence Diagram", font=TextFont), sg.Checkbox("Dimension 1", key="PD1", font=ButtonFont), sg.Checkbox("Dimension 2", key="PD2", font=ButtonFont)], [sg.Text("Accumulated Persistence Function", font=TextFont), sg.Checkbox("Dimension 1", key="APF1", font=ButtonFont), sg.Checkbox("Dimension 2", key="APF2", font=ButtonFont)],[sg.Button("Process", font=ButtonFont), sg.Button("Exit", font=ButtonFont)]]
 
	window_main = sg.Window("Topological Amorphous Material Analysis", layout_main, resizable = True, size=(600,750), right_click_menu=right_click_menu)
	loaded = False
	processed = False
	while True:	
		event_main, values_main = window_main.read()
		if (event_main == "Process"):
			if values_main["config_file"] != "":
				atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration(values_main["config_file"], values_main["structure"])
			else:	
				atoms = [str(a) for a in values_main['atoms'].split(",")]
				radii = [float(r) for r in values_main['radii'].split(",")]
				repeat_x = int(values_main["repeat_x"])
				repeat_y = int(values_main["repeat_y"])
				repeat_z = int(values_main["repeat_z"])
			births_1 = []
			deaths_1 = []
			births_2 = []
			deaths_2 = []
			APFs_1 = []
			APFs_2 = []
			sample_every = list(range(int(values_main["range-start"]), int(values_main["range-end"])+1, int(values_main["range-step"])))
			config_files = []
			for root, dirs, files in os.walk(values_main["file_dir"]):
				for f in files:
					print(f)
					if f.endswith(values_main["file_ext"]):
						config_files.append(os.path.join(root, f))
			print("Have found the following configuration files:")
			for f in config_files:
				print(f)
			births_1 = []
			deaths_1 = []
			births_2 = []
			deaths_2 = []
			APFs_1 = []
			APFs_2 = []
			for config in config_files:
				print("Looking at {}".format(config))
				atoms = load_atom_file(config, values_main["file_format"])
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
				for s in sample_every:
					points = load_atom_file(file_path, values_main["file_format"])
					s = int(values_main["sample_at"])
					points = sample_at(points, s, repeat_x, repeat_y, repeat_z, atoms, radii)
					simplices = weighted_alpha_diode(points)
					params.n_threads = int(values_main["n_threads"])
					if values_main["kernel"] or values_main["image"] or values_main["cokernel"]:
						if values_main["upper_threshold"] == "Automatic":
							upper_threshold = math.floor(max(points["z"])) 
						else:
							upper_threshold = float(values_main["upper_threshold"])
						if values_main["lower_threshold"] == "Automatic":
								lower_threshold = math.ceil(min(points["z"]))
						else:
							lower_threshold = float(values_main["lower_threshold"])
						params.kernel = values_main["kernel"]
						params.image = values_main["image"]
						params.cokernel = values_main["cokernel"]
						kicr, dgm_1, dgm_2 =  oineus_kernel_image_cokernel(points, params, upper_threshold, lower_threshold)
						dgm_1.to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.csv")
						dgm_2.to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_2.csv")
						if values_main["APF1"]:
							APF_1 = calculate_APF(dgm_1)
						if values_main["APF2"]:
							APF_2 = calculate_APF(dgm_2)	
					else:
						dgm_1, dgm_2 = oineus_process(points, params)
						dgm_1.to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_1.csv")	
						dgm_2.to_csv(dir+"/PD1/"+config_name+"_sample_"+str(s)+"_PD_2.csv")
						if values_main["APF1"]:
							APF_1 = calculate_APF(dgm_1)
							pandas.DataFrame(APF_1, columns=["birth", "death"]).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_1.csv")
						if values_main["APF2"]:
							APF_2 = calculate_APF(dgm_2)
							pandas.DataFrame(APF_2, columns=["birth", "death"]).to_csv(dir+"/APF1/"+config_name+"_sample_"+str(s)+"_APF_2.csv")
			print("All done!")
			processed = True

		if event_main == "Exit" or event_main == sg.WIN_CLOSED:
			window_main.close()
			entry_window()
			break
		
		if event_main == "Quit":
			break
	window_main.close()
 
 