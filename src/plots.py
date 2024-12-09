##
# @internal
# @file plots.py
# @brief functions to generate plots of the summaries
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

def plot_APF(APF : np.array, name : str):
	"""! Plot an accumulated persistence function
	
	@param APF - np.array with 2 columns of coordinates which define the APF
	@param name - title for the plot
	
	@result a plotly.express figure
	"""
	
	fig = px.line(x=APF["mean age"], y=APF["lifetime"], labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_APFs(APFs : list, APF_names : list, fig_name : str):#, APF_colour, APF_label):
	"""! Plot a set accumulated persistence function, with automatic colour differentiation.
	
	@param APFs - accumlated persistence functions to plot
	
	@result a matplotlib figure
	"""
	assert len(APFs) == len(APF_names)
	fig = go.Figure(labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=fig_name)
	last_pt = math.ceil(max([APFs[i][-1,0] for i in range(len(APFs))])*1.1)
	for i in range(len(APFs)):
		APFs[i] = np.vstack([APFs[i], [last_pt, APFs[i][-1,1]]])
	for i in range(len(APFs)):
		fig.add_trace(go.Scatter(x=APFs[i][:,0], y=APFs[i][:,1], mode="lines", name=APF_names[i]))
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_PD(dgm, name : str):
	"""! Plot a persistence diagram, with a specific colour
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.
	
	@param dgm 	- 	pandas.DataFrame of the diagram
	@param name    - name to use as title of the plot
	
	@result a plotly.express figure
	"""
	try:
		max_val = max(dgm["death"][dgm["death"] != math.inf])
	except:
		max_val = max(dgm["birth"])
	birth = []
	death = []
	inf_fin = []
	fig = go.Figure()
	for i in range(dgm.shape[0]):
		if dgm["death"].iloc[i] == math.inf:
			birth.append(dgm["birth"].iloc[i])
			death.append(max_val*1.1)
			inf_fin.append("inf")
		else:
			birth.append(dgm["birth"].iloc[i])
			death.append(dgm["death"].iloc[i])
			inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="inf_fin", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig


def plot_PDs(dgms, name : str):
	"""! Plot several persistence diagrams, with  automatic colour choices
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.

	@param dgms - list of diagrams
	@param name - title to use for the plot

	@results a plotly.express figure
	"""
	birth = []
	death = []
	samp = []
	inf_fin = []
	vals = []
	for dgm in dgms:
		dgm_vals = []
		for d in dgm["death"]:
			if d != math.inf:
				dgm_vals.append(d)
		if len(dgm_vals) !=0:
			vals.append(max(dgm_vals))
	if len(vals) != 0:
		max_val = max(vals)
	else:
		max_val = max([max(b) for b in dgm["birth"] for dgm in dmgs])
	fig = go.Figure()
	for i in range(len(dgms)):
		for j in range(len(dgms[i]["death"])):
			if dgms[i]["death"].iloc[j] == math.inf:
				birth.append(dgs[i]["birth"].iloc[j])
				death.append(max_val*1.1)
				samp.append(str(i))
				inf_fin.append("inf")
			else:
				birth.append(dgms[i]["birth"].iloc[j])
				death.append(dgms[i]["death"].iloc[j])
				samp.append(str(i))
				inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "sample":samp, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", color="sample", symbol="inf_fin", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig


def plot_kernel_image_cokernel_PD(kicr, d : int, codomain : bool, kernel : bool, image : bool, cokernel : bool, name : str):
	"""! Plot kernel, image, cokernel on same figure
	@param kicr 	oineus::KerImCokReduced 
	@param d	 	the dimension to extract (either 1 or 2)
	@param kernel	bool to plot kernel
	@param image	bool to plot image
	@param cokernel	bool to plot cokernel
	@return figu	figure with the chosen PD diagrams
	"""
	fig = go.Figure()
	max_val = -math.inf
	if codomain:
		codomain_pd = kicr.codomain_diagrams().in_dimension(d)
		print("codomain diagram has {} points".format(codomain_pd.shape[0]))
		if math.inf in codomain_pd[:,1] and max_val < max([d for d in codomain_pd[:,1] if d !=math.inf]):
			max_val = max(codomain_pd[:,1])
	if kernel:
		kernel_pd = kicr.kernel_diagrams().in_dimension(d)
		print("kernel diagram has {} points".format(kernel_pd.shape[0]))
		if math.inf in kernel_pd[:,1] and max_val < max([d for d in kernel_pd[:,1] if d !=math.inf]):
			max_val = max(kernel_pd[:,1])
	if image:
		image_pd = kicr.image_diagrams().in_dimension(d)
		print("image diagram has {} points".format(image_pd.shape[0]))
		if math.inf in image_pd[:,1]  and max_val < max([d for d in image_pd[:,1] if d !=math.inf]):
			max_val = max(image_pd[:,1])
	if cokernel:
		print("cokernel diagram has {} points".format(cokernel_pd.shape[0]))
		cokernel_pd = kicr.cokernel_diagrams().in_dimension(d)
		if math.inf in cokernel_pd[:,1] and max_val < max([d for d in cokernel_pd[:,1] if d !=math.inf]):
			max_val =  max(cokernel_pd[:,1])
	birth = []
	death = []
	pt_type = []
	inf_fin = []
	if codomain:
		for i in range(codomain_pd.shape[0]):
			if codomain_pd[i,1] == math.inf:
				birth.append(codomain_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("codomain")
				inf_fin.append("inf")
			else:
				birth.append(codomain_pd[i,0])
				death.append(codomain_pd[i,1])
				pt_type.append("codomain")
				inf_fin.append("fin")
	if kernel:
		for i in range(kernel_pd.shape[0]):
			if kernel_pd[i,1] == math.inf:
				birth.append(kernel_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("kernel")
				inf_fin.append("inf")
			else:
				birth.append(kernel_pd[i,0])
				death.append(kernel_pd[i,1])
				pt_type.append("kernel")
				inf_fin.append("fin")
	if image:
		for i in range(image_pd.shape[0]):
			if image_pd[i,1] == math.inf:
				birth.append(image_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("image")
				inf_fin.append("inf")
			else:
				birth.append(image_pd[i,0])
				death.append(image_pd[i,1])
				pt_type.append("image")
				inf_fin.append("fin")
	if cokernel:
		for i in range(cokernel_pd.shape[0]):
			if cokernel_pd[i,1] == math.inf:
				birth.append(cokernel_pd[i,0])
				death.append(max_val*1.1)
				pt_type.append("cokernel")
				inf_fin.append("inf")
			else:
				birth.append(cokernel_pd[i,0])
				death.append(cokernel_pd[i,1])
				pt_type.append("cokernel")
				inf_fin.append("fin")
	to_plot = pd.DataFrame({"birth":birth, "death":death, "pt_type":pt_type, "inf_fin":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="inf_fin", color="pt_type", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

# function to generate plots
def generate_plots():
	"""! Generate plots for a single configuration
	@brief Generate plots for a single configuration
	
	This function generates persistence diagrams and APF plots for a single configuration.
	It creates empty lists to store the plots, then populates them based on enabled options.
	
	The following plots can be generated:
	- Persistence diagrams in dimensions 0,1,2 if enabled via pd0, pd1, pd2 flags
	- Kernel/image/cokernel persistence diagrams if those computations are enabled
	- APF plots in dimensions 0,1,2 if enabled via apf0, apf1, apf2 flags
	- Kernel/image/cokernel APF plots if those computations are enabled
	
	The plots are stored in st.session_state for later display.
	
	@note Requires that persistent homology has already been computed and stored in session state
	"""
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
					st.session_state.fig_kic_pds_0.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 0, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 0 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 0.")
			st.session_state.fig_pds_0.append(plot_PD(st.session_state.dgms_0[i], st.session_state.file_path+" PD0 sample "+str(s)))
		if st.session_state["pd1"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_1.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 1, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 1 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 1.")
			st.session_state.fig_pds_1.append(plot_PD(st.session_state.dgms_1[i], st.session_state.file_path+" PD1 sample "+str(s)))
		if st.session_state["pd2"] == True:
			if st.session_state["kernel"] or st.session_state["image"] or st.session_state["cokernel"]:
				try:
					st.session_state.fig_kic_pds_2.append(plot_kernel_image_cokernel_PD(st.session_state.kicrs[i], 2, True, st.session_state["kernel"], st.session_state["image"], st.session_state["cokernel"], st.session_state.file_path+" codmain/kernel/image/cokernel dimension 2 sample "+str(s)))
				except:
					plot_tab.markdown("Encountered issues with kernel/image/cokernel diagram in dimension 2.")
			st.session_state.fig_pds_2.append(plot_PD(st.session_state.dgms_2[i], st.session_state.file_path+" PD2 sample "+str(s)))
		if st.session_state["apf0"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_0.append(plot_APF(st.session_state.kernel_APFs_0[i], st.session_state.file_path+" kernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 0.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_0.append(plot_APF(st.session_state.image_APFs_0[i], st.session_state.file_path+" image APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 0.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_0.append(plot_APF(st.session_state.cokernel_APFs_0[i], st.session_state.file_path+" cokernel APF0 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 0.")
			st.session_state.fig_apfs_0.append(plot_APF(st.session_state.APFs_0[i], st.session_state.file_path+" APF0 sample "+str(s)))
		if st.session_state["apf1"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_1.append(plot_APF(st.session_state.kernel_APFs_1[i], st.session_state.file_path+" kernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 1.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_1.append(plot_APF(st.session_state.image_APFs_1[i], st.session_state.file_path+" image APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 1.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_1.append(plot_APF(st.session_state.cokernel_APFs_1[i], st.session_state.file_path+" cokernel APF1 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 1.")
			st.session_state.fig_apfs_1.append(plot_APF(st.session_state.APFs_1[i], st.session_state.file_path+" APF1 sample "+str(s)))
		if st.session_state["apf2"]:
			if st.session_state["kernel"]:
				try:
					st.session_state.fig_kernel_apfs_2.append(plot_APF(st.session_state.kernel_APFs_2[i], st.session_state.file_path+" kernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute kernel APF in dimension 2.")
			if st.session_state["image"]:
				try:
					st.session_state.fig_image_apfs_2.append(plot_APF(st.session_state.image_APFs_2[i], st.session_state.file_path+" image APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute image APF in dimension 2.")
			if st.session_state["cokernel"]:
				try:
					st.session_state.fig_cokernel_apfs_2.append(plot_APF(st.session_state.cokernel_APFs_2[i], st.session_state.file_path+" cokernel APF2 sample "+str(s)))
				except:
					plot_tab.markdown("Can't compute cokernel APF in dimension 2.")
			st.session_state.fig_apfs_2.append(plot_APF(st.session_state.APFs_2[i], st.session_state.file_path+" APF2 sample "+str(s)))
	st.session_state.plots_generated = True