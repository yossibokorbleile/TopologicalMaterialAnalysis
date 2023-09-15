##
# @internal
# @file plots.py
# @brief Functions for generating plots.
# Given persistence diagrams and accumulated persistence functions, there are functions to plot either a single persistence diagram (PD) or accumulated persistence function (APF), or plot several together.
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from colour import Color
from matplotlib import cmfind $HOME -name "*activate" -type f
from matplotlib.colors import LogNorm
from scipy.interpolate import interpn
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets  import RectangleSelector
import plotly.express as px
import plotly.graph_objects as go

import numpy
import pandas

def plot_APF(APF : numpy.array, name : str):
	"""! Plot an accumulated persistence function
	
	@param APF - numpy.array with 2 columns of coordinates which define the APF
	@param name - title for the plot
	
	@result a plotly.express figure
	"""
	## Documentation for a function.
	# @param APF the apf to plot
	# @param APF_colour the colour to use
	
	fig = px.line(x=APF[:,0], y=APF[:,1], labels={'x':'m (Å$^2$)', 'y':'APF (Å$^2$)'}, title=name)
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
	return fig

def plot_PD(birth : list, death : list, name : str):
	"""! Plot a persistence diagram, with a specific colour
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.
	
	@param births - list of birth times
	@param deaths - list of death times
	@param name    - name to use as title of the plot
	
	@result a plotly.express figure
	"""
	assert len(birth) == len(death), f"Different number of sets of points provided."
	inf_fin = []
	vals = []
	plot_death = []
	for d in pts:
		if d != math.inf:
			vals.append(d)
	if len(vals) != 0:
		max_val = max(vals)
	else:
		max_val = max(births)
	fig = go.Figure()
	for i in range(len(births)):
		if death[i] == math.inf:
			plot_death.append(max_val*1.1)
			inf_fin.append("inf")
		else:
			plot_death.append(death[i])
			inf_fin.append("fin")
	to_plot = pandas.DataFrame({"birth":birth, "death":plot_death, "type":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", symbol="type", title=name)
	return fig


def plot_PDs(births : list, deaths : list, name : str):
	"""! Plot several persistence diagrams, with  automatic colour choices
	
	Points at infinity are plotted at a height of 1.1 times the last finite point to die.

	@param births - list of list of birth values
	@param deaths - list of list of death values
	@param name - title to use for the plot

	@results a plotly.express figure
	"""
	assert len(births) == len(deaths), f"Different number of sets of points provided."
	birth = []
	death = []
	samp = []
	inf_fin = []
	vals = []
	for pts in pds:
		dgm_vals = []
		for d in pts:
			if d != math.inf:
				dgm_vals.append(d)
		if len(dgm_vals) !=0:
			vals.append(max(dgm_vals))
	if len(vals) != 0:
		max_val = max(vals)
	else:
		max_val = max([max(pts) for pts in births])
	fig = go.Figure()
	for i in range(len(births)):
		for j in range(len(deaths[i])):
			if deaths[i][j] == math.inf:
				birth.append(births[i][j])
				death.append(max_val*1.1)
				samp.append(str(i))
				inf_fin.append("inf")
			else:
				birth.append(births[i][j])
				death.append(deaths[i][j])
				samp.append(str(i))
				inf_fin.append("fin")
	to_plot = pandas.DataFrame({"birth":birth, "death":death, "sample":samp, "type":inf_fin})
	fig = px.scatter(to_plot, x="birth", y="death", color="sample", symbol="type", title=name)
	return fig



def plot_kernel_image_cokernel_PD(kicr, d : int, kernel : bool, image : bool, cokernel : bool):
	"""! Plot kernel, image, cokernel on same figure
	@param kicr 	oineus::KerImCokReduced 
	@param d	 	the dimension to extract (either 1 or 2)
	@param kernel	bool to plot kernel
	@param image	bool to plot image
	@param cokernel	bool to plot cokernel
	@return figu	figure with the chosen PD diagrams
	"""
	print("settings are kerne {} image {} cokernel{}".format(kernel, image, cokernel))
	fig, ax = plt.subplots()
	max_val = -math.inf
	if kernel:
		kernel_pd = kicr.kernel_diagrams().in_dimension(d)
		if math.inf in kernel_pd[:,1] and max_val < max(kernel_pd[:,1]):
			max_val = max_val = max(kernel_pd[:,1])
	if image:
		image_pd = kicr.image_diagrams().in_dimension(d)
		if math.inf in image_pd[:,1]  and max_val < max(image_pd[:,1]):
			max_val  = max(image_pd[:,1])
	if cokernel:
		cokernel_pd = kicr.cokernel_diagrams().in_dimension(d)
		if math.inf in cokernel_pd[:,1] and max_val < max(cokernel_pd[:,1]):
			max_val =  max(cokernel_pd[:,1])
	if kernel:
		kernel_birth_fin = []
		kernel_birth_inf = []
		kernel_death_fin = []
		kernel_death_inf = []
		for i in range(kernel_pd.shape[0]):
			if kernel_pd[i,1] == math.inf:
				kernel_birth_inf.append(kernel_pd[i,0])
				kernel_death_inf.append(max_val*1.1)
			else:
				kernel_birth_fin.append(kernel_pd[i,0])
				kernel_death_fin.append(kernel_pd[i,1])
		ax.scatter(kernel_birth_fin, kernel_death_fin, color="#984ea3", s = 1, marker = ".", label="ker fin")	
		if len(kernel_birth_inf) !=0:
			ax.scatter(kernel_birth_inf, kernel_death_inf, color="#984ea3", s = 1, marker = "^", label="ker inf")
	if image:
		image_birth_fin = []
		image_birth_inf = []
		image_death_fin = []
		iamge_death_inf = []
		for i in range(image_pd.shape[0]):
			if image_pd[i,1] == math.inf:
				image_birth_inf.append(image_pd[i,0])
				image_death_inf.append(max_val*1.1)
			else:
				image_birth_fin.append(image_pd[i,0])
				image_death_fin.append(image_pd[i,1])
		ax.scatter(image_birth_fin, image_death_fin, color="#ff7f00", s = 1, marker = ".", label="im fin")	
		if len(image_birth_inf) !=0:
			ax.scatter(image_birth_inf, image_death_inf, color="#ff7f00", s = 1, marker = "^", label="im inf")
	if cokernel:
		cokernel_birth_fin = []
		cokernel_birth_inf = []
		cokernel_death_fin = []
		cokernel_death_inf = []
		for i in range(cokernel_pd.shape[0]):
			if cokernel_pd[i,1] == math.inf:
				cokernel_birth_inf.append(cokernel_pd[i,0])
				cokernel_death_inf.append(max_val*1.1)
			else:
				cokernel_birth_fin.append(cokernel_pd[i,0])
				cokernel_death_fin.append(cokernel_pd[i,1])
		ax.scatter(cokernel_birth_fin, cokernel_death_fin, color="#a65628", s = 1, marker = ".", label="cok fin")	
		if len(kernel_birth_inf) !=0:
			ax.scatter(cokernel_birth_inf, cokernel_death_inf, color="#a65628", s = 1, marker = "^", label="cok inf")
	ax.set_xlabel('birth')
	ax.set_ylabel('death')
	ax.legend()
	fig.tight_layout(pad=5.0)
	return fig
	

def draw_figure(canvas, figure):
	figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
	figure_canvas_agg.draw()
	figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
	return figure_canvas_agg

def layout_plot_sample_at(name : str, object : str, sample_at):
	layout = [[sg.Text("{} Plot: {} at sample {}".format(name, object,sample_at), font="Arial 20")],[sg.T('Controls:')], 
				[sg.Canvas(key='controls_cv')],
				[sg.Text('Figure:')],
				[sg.Column(layout=[
				[sg.Canvas(key='figCanvas',
					size=(700, 700)
					)]
				],
				background_color='#DAE0E6',
				pad=(0, 0)
				)]]
	return layout

def layout_plot(name : str, object : str):
	layout = [[sg.Text("{} Plot: {}".format(name, object), font="Arial 20")],[sg.T('Controls:')], 
				[sg.Canvas(key='controls_cv')],
				[sg.Text('Figure:')],
				[sg.Column(layout=[
				[sg.Canvas(key='figCanvas',
					size=(00, 700)
					)]
				],
				background_color='#DAE0E6',
				pad=(0, 0)
				)]]
	return layout

def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
	if canvas.children:
		for child in canvas.winfo_children():
			child.destroy()
	if canvas_toolbar.children:
		for child in canvas_toolbar.winfo_children():
			child.destroy()
	figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
	figure_canvas_agg.draw()
	toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
	toolbar.update()
	figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)


class Toolbar(NavigationToolbar2Tk):
	def __init__(self, *args, **kwargs):
		super(Toolbar, self).__init__(*args, **kwargs)
