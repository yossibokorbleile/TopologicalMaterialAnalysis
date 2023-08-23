##
# @internal
# @file plots.py
# @brief Functions for generating plots.
# Given persistence diagrams and accumulated persistence functions, there are functions to plot either a single persistence diagram (PD) or accumulated persistence function (APF), or plot several together.
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from colour import Color
from matplotlib import cm
from matplotlib.colors import LogNorm
from scipy.interpolate import interpn
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets  import RectangleSelector

import numpy


def plot_APF(APF : numpy.array, APF_colour : str):#, APF_colour, APF_label):
    """! Plot an accumulated persistence function
    
    @param APF - numpy.array with 2 columns of coordinates which define the APF
    @param APF_colour - colour to use
    
    @result a matplotlib figure
    """
    ## Documentation for a function.
    # @param APF the apf to plot
    # @param APF_colour the colour to use
    
    fig, ax = plt.subplots()
    ax.plot(APF[:,0], APF[:,1], color=APF_colour)
    ax.set_xlabel('m (Å$^2$)')
    ax.set_ylabel('APF (Å$^2$)')
    fig.tight_layout(pad=5.0)
    return fig

def plot_APFs(APFs):#, APF_colour, APF_label):
    """! Plot a set accumulated persistence function, with automatic colour differentiation.
    
    @param APFs - accumlated persistence functions to plot
    
    @result a matplotlib figure
    """
    last_pt = math.ceil(max([APFs[i][-1,0] for i in range(len(APFs))])*1.1)
    for i in range(len(APFs)):
        APFs[i] = numpy.vstack([APFs[i], [last_pt, APFs[i][-1,1]]])
    fig, ax = plt.subplots()
    for i in range(len(APFs)):
        c = Color("blue")
        c.red=i/len(APFs)
        c.blue =  1 - i/len(APFs)
        ax.plot(APFs[i][:,0], APFs[i][:,1], color=c.rgb)
    ax.set_xlabel('m (Å$^2$)')
    ax.set_ylabel('APF')
    fig.tight_layout(pad=5.0)
    return fig

def plot_PD(births : list, deaths : list, PD_colour : str):
    """! Plot a persistence diagram, with a specific colour
    
    Points at infinity are plotted at a height of 1.1 times the last finite point to die.
    
    @param births - list of birth times
    @param deaths - list of death times
    @param PD_colour    - colour to use
    
    @result a matplotlib figure
    """
    fig, ax = plt.subplots()
    if math.inf in deaths:
        max_val = max(deaths)
        birth_fin = []
        death_fin = []
        birth_inf = []
        death_inf = []
        for i in range(len(deaths)):
            if deaths[i] == math.inf:
                birth_inf.append(births[i])
                death_inf.append(max_val*1.1)
            else:
                birth_fin.append(births[i])
                death_fin.append(deaths[i])
        ax.scatter(birth_fin, death_fin, color=PD_colour, s = 1, marker = ".")
        ax.scatter(birth_inf, death_inf, color=PD_colour, s = 1, marker = "^")
    else:
        ax.scatter(births, deaths, color=PD_colour, s= 1, marker = ".")
    ax.set_xlabel('birth')
    ax.set_ylabel('death')
    fig.tight_layout(pad=5.0)
    return fig

def plot_PDs(births : list, deaths : list):
    """! Plot several persistence diagrams, with  automatic colour choices
    """
    fig, ax = plt.subplots()
    for i in range(len(births)):
        c = Color("blue")
        c.red=i/len(births)
        c.blue =  1 - i/len(births)
        if math.inf in deaths[i]:
            max_val = max(deaths[i])
            birth_fin = []
            death_fin = []
            birth_inf = []
            death_inf = []
            for j in range(len(deaths[i])):
                if deaths[i][j] == math.inf:
                    birth_inf.append(births[i][j])
                    death_inf.append(max_val*1.1)
                else:
                    birth_fin.append(births[i][j])
                    death_fin.append(deaths[i][j])
            ax.scatter(birth_fin, death_fin, color=c.rgb, s = 1, marker = ".")
            ax.scatter(birth_inf, death_inf, color=c.rgb, s = 1, marker = "^")
        else:
            ax.scatter(births[i], deaths[i], color=c.rgb, s=1, marker=".")
    ax.set_xlabel('birth')
    ax.set_ylabel('death')
    fig.tight_layout(pad=5.0)
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
