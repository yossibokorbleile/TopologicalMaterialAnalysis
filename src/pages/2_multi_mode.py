##
# @internal
# @file multi_mode.py
# @brief Streamlit page for analysing a multiple samples.
# @version 0.1
# @date July 2024


import sys
sys.path.insert(0, '..')
import streamlit as st
import oineus
import numpy
import pandas
import os
import process
import streamlit_functions
# import plots
import visualisation
from ase import io
import sys

st.header("Multi Mode")
comp_tab, plot_tab, vis_tab = st.tabs(["Computation", "Plots", "Visuatlisation"]) #create tabs for the various parts
st.session_state.mode = "multi"


###define various functions needed for later
def test():
	st.session_state["maual_comp_config"] = True
	st.session_state.processed=True
	st.session_state.config_file = "../examples/structure-types.ini"
	st.session_state.file_path = "../examples/ZIF_test.xyz"
	st.session_state.config_name = "ZIF-TEST"
	st.session_state.comp_name = "ZIF-test"
	st.session_state.sample_start = 0
	st.session_state.sample_end = 2
	st.session_state.sample_step = 1
	st.session_state.repeat_x = 1
	st.session_state.repeat_y = 1
	st.session_state.repeat_z = 1
	st.session_state.thickness=0.1
	st.session_state.kernel = True
	st.session_state.image = True
	st.session_state.cokernel = True
	compute()

# Function to compute the persistent homology 
def compute():
	st.session_state.params = oineus.ReductionParams()
	if not st.session_state["manual_config"]:
		streamlit_functions.load_configuration_settings()
	if not st.session_state["maual_comp_config"]:
		streamlit_functions.load_computation_settings()
	st.session_state.sample_indices = []
	st.session_state.dgms_0 = []
	st.session_state.dgms_1 = []
	st.session_state.dgms_2 = []
	st.session_state.APFs_0 = []
	st.session_state.APFs_1 = []
	st.session_state.APFs_2 = []
	st.session_state.atom_locations_list = []
	st.session_state.dcmps = [] 
	st.session_state.filts = []
	if "kernel" not in st.session_state:
		st.session_state.kernel = False
		st.session_state.params.kernel = False
	else:
		print("KERNEL is ", st.session_state.kernel)
		st.session_state.params.kernel = st.session_state.kernel
	if "image" not in st.session_state:
		st.session_state.image = False
		st.session_state.params.image = False
	else:
		st.session_state.params.image= st.session_state.image
	if "cokernel" not in st.session_state:
		st.session_state.cokernel = False
		st.session_state.params.cokernel = False
	else:
		st.session_state.params.cokernel = st.session_state.cokernel
	if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
		st.session_state.kicrs = []
	if st.session_state["kernel"]:
		st.session_state.kernel_dgms_0 = []
		st.session_state.kernel_dgms_1 = []
		st.session_state.kernel_dgms_2 = []
		st.session_state.kernel_APFs_0 = []
		st.session_state.kernel_APFs_1 = []
		st.session_state.kernel_APFs_2 = []
	if st.session_state["image"]:
		st.session_state.image_dgms_0 = []
		st.session_state.image_dgms_1 = []
		st.session_state.image_dgms_2 = []
		st.session_state.image_APFs_0 = []
		st.session_state.image_APFs_1 = []
		st.session_state.image_APFs_2 = []
	if st.session_state["cokernel"]:
		st.session_state.cokernel_dgms_0 = []
		st.session_state.cokernel_dgms_1 = []
		st.session_state.cokernel_dgms_2 = []
		st.session_state.cokernel_APFs_0 = []
		st.session_state.cokernel_APFs_1 = []
		st.session_state.cokernel_APFs_2 = []
	if st.session_state.sample_end == "Auto":
		st.session_state.sample_end = 1
	for s in range(st.session_state.sample_start, st.session_state.sample_end, st.session_state.sample_step):
		st.session_state.sample_index = s
		st.session_state.sample_indices.append(s)
		st.session_state.atom_locations = process.sample_at(st.session_state.file_path, st.session_state.file_format, st.session_state.sample_index, st.session_state.repeat_x, st.session_state.repeat_y, st.session_state.repeat_z, st.session_state.atoms, st.session_state.radii)
		st.session_state.atom_locations_list.append(st.session_state.atom_locations)
		if st.session_state.params.kernel or st.session_state.params.image or st.session_state.params.cokernel:
			top_pt = max(st.session_state.atom_locations["z"])
			bot_pt = min(st.session_state.atom_locations["z"])
			height = abs(top_pt - bot_pt)
			if st.session_state.thickness == "Auto":
				ut= top_pt - 0.1*height
				lt = bot_pt + 0.1*height
			else:
				ut = top_pt - st.session_state.thickness*height
				lt = bot_pt + st.session_state.thickness*height
			st.session_state.kicr, st.session_state.dgm_0, st.session_state.dgm_1, st.session_state.dgm_2 = process.oineus_kernel_image_cokernel(st.session_state.atom_locations, st.session_state.params, ut, lt)
			st.session_state.kicrs.append(st.session_state.kicr)
			st.session_state.dgms_0.append(st.session_state.dgm_0)
			st.session_state.dgms_1.append(st.session_state.dgm_1)
			st.session_state.dgms_2.append(st.session_state.dgm_2)
			st.session_state.APFs_0.append(process.calculate_APF(st.session_state.dgm_0))
			st.session_state.APFs_1.append(process.calculate_APF(st.session_state.dgm_1))
			st.session_state.APFs_2.append(process.calculate_APF(st.session_state.dgm_2))
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				st.session_state.kicrs.append(st.session_state.kicr)
			if st.session_state["kernel"]:
				kernel_dgm_0 = pandas.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(0), columns=["birth", "death"])
				kernel_dgm_1 = pandas.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(1), columns=["birth", "death"])
				kernel_dgm_2 = pandas.DataFrame(st.session_state.kicr.kernel_diagrams().in_dimension(2), columns=["birth", "death"])
				print("kernel_dgm_2 is:")
				print(kernel_dgm_2)
				st.session_state.kernel_dgms_0.append(kernel_dgm_0)
				st.session_state.kernel_dgms_1.append(kernel_dgm_1)
				st.session_state.kernel_dgms_2.append(kernel_dgm_2)
				st.session_state.kernel_APFs_0.append(process.calculate_APF(kernel_dgm_0))
				st.session_state.kernel_APFs_1.append(process.calculate_APF(kernel_dgm_1))
				st.session_state.kernel_APFs_2.append(process.calculate_APF(kernel_dgm_2))
			if st.session_state["image"]:
				image_dgm_0 = pandas.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(0), columns=["birth", "death"])
				image_dgm_1 = pandas.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(1), columns=["birth", "death"])
				image_dgm_2 = pandas.DataFrame(st.session_state.kicr.image_diagrams().in_dimension(2), columns=["birth", "death"])
				st.session_state.image_dgms_0.append(image_dgm_0)
				st.session_state.image_dgms_1.append(image_dgm_1)
				st.session_state.image_dgms_2.append(image_dgm_2)
				st.session_state.image_APFs_0.append(process.calculate_APF(image_dgm_0))
				st.session_state.image_APFs_1.append(process.calculate_APF(image_dgm_1))
				st.session_state.image_APFs_2.append(process.calculate_APF(image_dgm_2))
			if st.session_state["cokernel"]:
				cokernel_dgm_0 = pandas.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				cokernel_dgm_1 = pandas.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				cokernel_dgm_2 = pandas.DataFrame(st.session_state.kicr.cokernel_diagrams().in_dimension(0), columns=["birth", "death"])
				st.session_state.cokernel_dgms_0.append(cokernel_dgm_0)
				st.session_state.cokernel_dgms_1.append(cokernel_dgm_1)
				st.session_state.cokernel_dgms_2.append(cokernel_dgm_2)
				st.session_state.cokernel_APFs_0.append(process.calculate_APF(cokernel_dgm_0))
				st.session_state.cokernel_APFs_1.append(process.calculate_APF(cokernel_dgm_1))
				st.session_state.cokernel_APFs_2.append(process.calculate_APF(cokernel_dgm_2))
		else:
			st.session_state.dcmp, st.session_state.filt, st.session_state.dgm_0, st.session_state.dgm_1, st.session_state.dgm_2 = process.oineus_process(st.session_state.atom_locations, st.session_state.params)
			st.session_state.dgms_0.append(st.session_state.dgm_0)
			st.session_state.dgms_1.append(st.session_state.dgm_1)
			st.session_state.dgms_2.append(st.session_state.dgm_2)
			st.session_state.APFs_0.append(process.calculate_APF(st.session_state.dgm_0))
			st.session_state.APFs_1.append(process.calculate_APF(st.session_state.dgm_1))
			st.session_state.APFs_2.append(process.calculate_APF(st.session_state.dgm_2))
			st.session_state.dcmps.append(st.session_state.dcmp) 
			st.session_state.filts.append(st.session_state.filt)

	if len(st.session_state.sample_indices) != len(st.session_state.dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.APFs_2):
		st.markdown("*WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
	if "kernel_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.kernel_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.kernel_APFs_2):
			st.markdown("KERNEL *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
			st.write(len(st.session_state.kernel_dgms_0))
			st.write(len(st.session_state.sample_indices))
	if "image_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.image_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.image_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.image_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.image_APFs_2):
			st.markdown("IMAGE *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.") 
	if "cokernel_dgms_0" in st.session_state:
		if len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_0) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_1) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_dgms_2) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_0) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_1) or len(st.session_state.sample_indices) != len(st.session_state.cokernel_APFs_2):
			st.markdown("COKERNEL *WARNING* something went wrong, the number of diagrams/APFs calcualted does not match the number expected.")
	st.session_state["processed_file"] =  st.session_state.file_path
	st.session_state.processed = True
	# print("Finished!")

# function to generate plots
def generate_plots():
	st.session_state.fig_pds_0 = []
	st.session_state.fig_pds_1 = []
	st.session_state.fig_pds_2 = []
	st.session_state.fig_apfs_0 = []
	st.session_state.fig_apfs_1 = []
	st.session_state.fig_apfs_2 = []
	if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
		st.session_state.fig_kic_pds_0 = []
		st.session_state.fig_kic_pds_1 = []
		st.session_state.fig_kic_pds_2 = []
		st.session_state.fig_pds_1 = []
		st.session_state.fig_pds_2 = []
		st.session_state.fig_kernel_apfs_0 = []
		st.session_state.fig_kernel_apfs_1 = []
		st.session_state.fig_kernel_apfs_2 = []
		st.session_state.fig_image_apfs_0 = []
		st.session_state.fig_image_apfs_1 = []
		st.session_state.fig_image_apfs_2 = []
		st.session_state.fig_cokernel_apfs_0 = []
		st.session_state.fig_cokernel_apfs_1 = []
		st.session_state.fig_cokernel_apfs_2 = []
		st.session_state.fig_apfs_0 = []
		st.session_state.fig_apfs_1 = []
		st.session_state.fig_apfs_2 = []
	if "kernel" not in st.session_state:
		st.session_state["kernel"] = False
	if "image" not in st.session_state:
		st.session_state["image"] = False
	if "cokernel" not in st.session_state:
		st.session_state["cokernel"] = False
	for i,s in enumerate(st.session_state.sample_indices):
		if st.session_state["pd0"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_0.append(plots.plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 0, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 0 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 0.")
			st.session_state.fig_pds_0.append(plots.plot_PD(st.session_state.dgms_0[i], st.session_state.file_path+" PD0 sample "+str(s)))
		if st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_1.append(plots.plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 1, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 1 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 1.")
			st.session_state.fig_pds_1.append(plots.plot_PD(st.session_state.dgms_1[i], st.session_state.file_path+" PD1 sample "+str(s)))
		if st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_2.append(plots.plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 2, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 2 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 2.")
			st.session_state.fig_pds_2.append(plots.plot_PD(st.session_state.dgms_2[i], st.session_state.file_path+" PD2 sample "+str(s)))
		if st.session_state["apf0"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_0.append(plots.plot_APF(st.session_state.kernel_APFs_0[i], file_path+" kernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 0.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_0.append(plots.plot_APF(st.session_state.image_APFs_0[i], file_path+" image APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 0.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_0.append(plots.plot_APF(st.session_state.cokernel_APFs_0[i], file_path+" cokernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 0.")
			st.session_state.fig_apfs_0.append(plots.plot_APF(st.session_state.APFs_0[i], file_path+" APF1 sample "+str(s)))
		if st.session_state["apf1"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_1.append(plots.plot_APF(st.session_state.kernel_APFs_1[i], file_path+" kernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 1.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_1.append(plots.plot_APF(st.session_state.image_APFs_1[i], file_path+" image APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 1.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_1.append(plots.plot_APF(st.session_state.cokernel_APFs_1[i], file_path+" cokernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 1.")
			st.session_state.fig_apfs_1.append(plots.plot_APF(st.session_state.APFs_1[i], file_path+" APF1 sample "+str(s)))
		if st.session_state["apf2"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_2.append(plots.plot_APF(st.session_state.kernel_APFs_2[i], file_path+" kernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 2.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_2.append(plots.plot_APF(st.session_state.image_APFs_2[i], file_path+" image APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 2.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_2.append(plots.plot_APF(st.session_state.cokernel_APFs_2[i], file_path+" cokernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 2.")
			st.session_state.fig_apfs_2.append(plots.plot_APF(st.session_state.APFs_2[i], file_path+" APF2 sample "+str(s)))
	st.session_state.plots_generated = True

#function to display plots
def display_plots():
	if "plots_generated" not in st.session_state or not st.session_state["plots_generated"]:
		generate_plots()
	plot_tab.write(st.session_state.sample_end)
	plot_tab.write(st.session_state.sample_indices)
	for i,s in enumerate(st.session_state.sample_indices):
		if st.session_state["pd0"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_0[i])
			# else: 
			plot_tab.plotly_chart(st.session_state.fig_pds_0[i])
		if st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_1[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_1[i])
		if st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				plot_tab.plotly_chart(st.session_state.fig_kic_pds_2[i])
			# else:
			plot_tab.plotly_chart(st.session_state.fig_pds_2[i])
		if st.session_state["apf0"]:
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
		if st.session_state["apf1"]:
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
		if st.session_state["apf2"]:
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


# function to save plots
def save_plots():
	dir_name = os.path.dirname(file_path)
	file_name = os.path.splitext(os.path.split(file_path)[1])[0]
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


### lets set up the computation tab
file_path = comp_tab.text_input("Initial structure file:",key="file_path") #specify initial structure file 
file_format = comp_tab.text_input("File format:", key="file_format",placeholder="Auto") #specify format of the initial strutcure file
if "processed_file" in st.session_state:
	if file_path != st.session_state["processed_file"]:
		st.session_state.processed = False

comp_tab.header("Configuration settings")
manual_config = comp_tab.checkbox("Manually specify configuration", key="manual_config")#manually set configuration
if not manual_config:
	st.session_state.config_file = comp_tab.text_input("Configuration file:", key="configuration_file")
	st.session_state.config_name = comp_tab.text_input("Configuration name:", key="configuration_name")
else:
	st.session_state.atoms = comp_tab.text_input("Atoms:", key="atoms_input")
	st.session_state.radii = comp_tab.text_input("Radii:", key="radii_input")


comp_tab.markdown("Computation settings")
manual_compute = comp_tab.checkbox("Manually specify settings for the computations (i.e number of threds, and if you want to compute kernel/image/cokernel)", key="maual_comp_config")
if not manual_compute:
	same_config_file = comp_tab.checkbox("The computation settings are in the same configuration file.", key="same_config_file")
	if not same_config_file:
		st.session_state.comp_file = comp_tab.text_input("Configuration file:", key="comp_config_file")
	else:
		st.session_state.comp_file = st.session_state.config_file
	st.session_state.comp_name = comp_tab.text_input("Configuration name:", key="comp_config_name")
else:
	st.session_state.sample_start = comp_tab.text_input("Sample start", key="sample_start_input")
	st.session_state.sample_end = comp_tab.text_input("Sample end", key="sample_end_input")
	st.session_state.sample_step = comp_tab.text_input("Sample step", key="sample_step_input")
	st.session_state.repeat_x = comp_tab.text_input("Repitition in x-axis:", key="repeat_x_input")
	st.session_state.repeat_y = comp_tab.text_input("Repitition in y-axis:", key="repeat_y_input")
	st.session_state.repeat_z = comp_tab.text_input("Repitition in z-axis:", key="repeat_z_input")
	st.session_state.kernel = comp_tab.checkbox("Compute kernel persistence", key="kernel_check")
	st.session_state.image = comp_tab.checkbox("Compute image persistence", key="image_check")
	st.session_state.cokernel = comp_tab.checkbox("Compute cokernel persistence", key="cokernel_check")
	st.session_state.thickness = comp_tab.text_input("Select thickness of top and bottom layer:", key="thickness_input", placeholder="Automatic detection")
	st.session_state.n_threads = comp_tab.text_input("Select number of threads to use:", key="n_threads_input", placeholder="4")
	if st.session_state.n_threads == "":
		st.session_state.params.n_threads = 4
	else:
		st.session_state.params.n_threads = int(st.session_state.n_threads)
	comp_tab.markdown(f"Number of threads is".format(st.session_state.params.n_threads))

	if st.session_state["thickness_input"] == "":
		st.session_state.thickness = "Auto"
	else:
		st.session_state.thickness= float(st.session_state["thickness_input"])

comp_tab.markdown(f"The file selected to analyse is "+file_path)	
if file_format == "":
	comp_tab.markdown("File format will automatically detected.")
else:
	comp_tab.markdown("File format is", file_format)


comp_tab.button("Process", key="process", on_click=compute)


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
if "processed" in st.session_state and st.session_state["processed"] and st.session_state["visualisation"] and not st.session_state.params.kernel and not st.session_state.params.image and not st.session_state.params.cokernel:
	vis_tab.write(st.session_state.params)
	vis_tab.write(st.session_state.dcmps)
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
	st.session_state.dfVis = visualisation.generate_visulisation_df(st.session_state.vis_dgm, st.session_state.dcmps[st.session_state.selected_sample_index].r_data, st.session_state.filts[st.session_state.selected_sample_index], st.session_state.atom_locations_list[st.session_state.selected_sample_index], st.session_state.atoms)
	to_display = ["birth", "death", "lifetime"]
	for a in st.session_state.atoms:
		to_display.append(a)
	viz = vis_tab.dataframe(st.session_state.dfVis[to_display], on_select="rerun")
	if viz.selection.rows == []:
		st.session_state.selected_row = None
	else:
		st.session_state.selected_row = viz.selection.rows
	# vis_tab.write("Selected Row:")
	# vis_tab.write(st.session_state.selected_row)
	if st.session_state.selected_row != None:
		for cycle_id in st.session_state.selected_row:
			vis_tab.plotly_chart(visualisation.generate_display(st.session_state.atom_locations, st.session_state.dgm_1, cycle_id, st.session_state.filt, neighbours = st.session_state["neighbours"]))
elif st.session_state.params.kernel or st.session_state.params.image or st.session_state.params.cokernel:
	vis_tab.markdown("Visulation of kernel/image/cokernel persistent homology is not yet available.")
else:
	vis_tab.markdown("Persistent homology has not been computed, so representative cycles are unknow. Please proces the file, and then return to this tab.")

st.button("test", key="test", on_click=test)

if "params" in st.session_state:
	print(st.session_state.params)

