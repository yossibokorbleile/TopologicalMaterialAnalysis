##
# @internal
# @file plots.py
# @brief Functions for generating plots.
# Given persistence diagrams and accumulated persistence functions, there are functions to plot either a single persistence diagram (PD) or accumulated persistence function (APF), or plot several together.
#import matplotlib.pyplot as plt
#import matplotlib as mpl
import math
from colour import Color
#from matplotlib import cmfind
#from matplotlib.colors import LogNorm
from scipy.interpolate import interpn
#from matplotlib.ticker import MaxNLocator
#from matplotlib.ticker import FuncFormatter
#import matplotlib.colors as mcolors
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
#from matplotlib.widgets  import RectangleSelector
import plotly.express as px
import plotly.graph_objects as go

import PySimpleGUI as sg	

import numpy
import pandas

def plot_APF(APF : numpy.array, name : str):
	"""! Plot an accumulated persistence function
	
	@param APF - numpy.array with 2 columns of coordinates which define the APF
	@param name - title for the plot
	
	@result a plotly.express figure
	"""
	
	fig = px.line(x=APF[:,0], y=APF[:,1], labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_APFs(APFs : list, name : str):#, APF_colour, APF_label):
	"""! Plot a set accumulated persistence function, with automatic colour differentiation.
	
	@param APFs - accumlated persistence functions to plot
	
	@result a matplotlib figure
	"""
	fig = go.Figure(labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=name)
	last_pt = math.ceil(max([APFs[i][-1,0] for i in range(len(APFs))])*1.1)
	for i in range(len(APFs)):
		APFs[i] = numpy.vstack([APFs[i], [last_pt, APFs[i][-1,1]]])
	for i in range(len(APFs)):
		fig.add_trace(go.Scatter(x=APFs[i][:,0], y=APFs[i][:,1], mode="lines", name=str(i)))
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig

def plot_PD(dgm, name : str):
	"""! Plot a persistence diagram, with a specific colour
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.
	
	@param dgm 	- 	numpy.array of size nx2
	@param name    - name to use as title of the plot
	
	@result a plotly.express figure
	"""
	max_val = max(dgm["death"][dgm["death"] != math.inf])
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
	to_plot = pandas.DataFrame({"birth":birth, "death":death, "type":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="type", title=name)
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
	#assert len(births) == len(deaths), f"Different number of sets of points provided."
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
	for i in range(len(dgm)):
		for j in range(len(dgms[i]["death"])):
			if dgms[i]["death"].iloc[j] == math.inf:
				birth.append(dgs[i]["birth"].iloc[j])
				death.append(max_val*1.1)
				samp.append(str(i))
				inf_fin.append("inf")
			else:
				birth.append(dgs[i]["birth"].iloc[j])
				death.append(dgs[i]["death"].iloc[j])
				samp.append(str(i))
				inf_fin.append("fin")
	to_plot = pandas.DataFrame({"birth":birth, "death":death, "sample":samp, "type":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", color="sample", symbol="type", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig



def plot_kernel_image_cokernel_PD(kicr, d : int, kernel : bool, image : bool, cokernel : bool, name : str):
	"""! Plot kernel, image, cokernel on same figure
	@param kicr 	oineus::KerImCokReduced 
	@param d	 	the dimension to extract (either 1 or 2)
	@param kernel	bool to plot kernel
	@param image	bool to plot image
	@param cokernel	bool to plot cokernel
	@return figu	figure with the chosen PD diagrams
	"""
	print("settings are kernel {} image {} cokernel {}".format(kernel, image, cokernel))
	#fig, ax = plt.subplots()
	fig = go.Figure()
	max_val = -math.inf
	if kernel:
		kernel_pd = kicr.kernel_diagrams().in_dimension(d)
		if math.inf in kernel_pd[:,1] and max_val < max(kernel_pd[:,1]):
			max_val = max(kernel_pd[:,1])
	if image:
		image_pd = kicr.image_diagrams().in_dimension(d)
		if math.inf in image_pd[:,1]  and max_val < max(image_pd[:,1]):
			max_val = max(image_pd[:,1])
	if cokernel:
		cokernel_pd = kicr.cokernel_diagrams().in_dimension(d)
		if math.inf in cokernel_pd[:,1] and max_val < max(cokernel_pd[:,1]):
			max_val =  max(cokernel_pd[:,1])
	birth = []
	death = []
	ker_im_cok = []
	inf_fin = []
	if kernel:
		for i in range(kernel_pd.shape[0]):
			if kernel_pd[i,1] == math.inf:
				birth.append(kernel_pd[i,0])
				death.append(max_val*1.1)
				ker_im_cok.append("kernel")
				inf_fin.append("inf")
			else:
				birth.append(kernel_pd[i,0])
				death.append(kernel_pd[i,1])
				ker_im_cok.append("kernel")
				inf_fin.append("fin")
	if image:
		for i in range(image_pd.shape[0]):
			if image_pd[i,1] == math.inf:
				birth.append(image_pd[i,0])
				death.append(max_val*1.1)
				ker_im_cok.append("image")
				inf_fin.append("inf")
			else:
				birth.append(image_pd[i,0])
				death.append(image_pd[i,1])
				ker_im_cok.append("image")
				inf_fin.append("fin")
	if cokernel:
		for i in range(cokernel_pd.shape[0]):
			if cokernel_pd[i,1] == math.inf:
				birth.append(cokernel_pd[i,0])
				death.append(max_val*1.1)
				ker_im_cok.append("cokernel")
				inf_fin.append("inf")
			else:
				birth.append(cokernel_pd[i,0])
				death.append(cokernel_pd[i,1])
				ker_im_cok.append("cokernel")
				inf_fin.append("fin")
	#fig.tight_layout(pad=5.0)
	to_plot = pandas.DataFrame({"birth":birth, "death":death, "ker_im_cok":ker_im_cok, "type":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="type", color="ker_im_cok", title=name)
	fig.update_xaxes(rangemode="tozero")
	fig.update_yaxes(rangemode="tozero")
	return fig
