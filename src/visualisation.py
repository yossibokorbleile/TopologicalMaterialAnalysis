##
# @internal
# @file visualisation.py
# @brief Look for represetnatives of the cycles.

import numpy
import pandas
import math
from ase import Atoms
import plotly.express as px
import plotly.graph_objects as go

def get_representative_loops(dgm : pandas.DataFrame, V, filt):
	"""! Get representative of each homology class in dimension 1.
	@param dgm 		pandas.DataFrame containing the diagram and the birth and death simplex id
	@param filt		oineus filtration
	
	@return dgm		pandas.DataFrame with new column `cycle rep`
	"""
	cycle_reps = []
	print(dgm.shape[0])
	for i in range(dgm.shape[0]):
		idx = dgm["death simplex"].iloc[i]
		sorted_rep = V[filt.get_sorted_id_by_id(idx)]
		rep = []
		for v in sorted_rep:
			rep.append(filt.get_id_by_sorted_id(v))
		cycle_reps.append(rep)
		print("{} ({}, {}): unsorted birth_id {} sorted birth id {} rep sorted {} rep unsorted {}".format(i,dgm["birth"].iloc[i], dgm["death"].iloc[i],idx,filt.get_sorted_id_by_id(idx),V[filt.get_sorted_id_by_id(idx)], rep))
	dgm["cycle rep"] = cycle_reps
	return dgm

def loop_composition(loop, points, atom_types):
	"""! Get the composition of a given representative 

	@param loop		list containing the ids of the vertices in the loop
	@param points	pandas.DataFrame of the points (atoms) in the structure
	@param atom_types list of the atom types we are considering (restricting to)

	@return comp	dictionary listing the number of atoms of each type in the loop
	"""
	comp = dict([(a, 0) for a in atom_types])
	for x in loop:
		comp[points["Atom"].iloc[x]] += 1
	return comp

def generate_visulisation_df(dgm : pandas.DataFrame, V, filt, points, atom_types):  
	"""! generate the pandas.DataFrame containing the information about the points so we can display it 

	@param dgm		pandas.DataFrame of the diagram, with columns `birth`, `death`, `birth simplex`, `death simplex`, `cycle rep`
	"""
	dgm =  get_representative_loops(dgm, V, filt)
	#get the composition of the cycle representatives
	comps = []
	for i in range(dgm.shape[0]):
		comps.append(loop_composition(dgm["cycle rep"].iloc[i], points, atom_types))
	#for each atom type we are looking at, add a column with the number of atoms of this type in the cycle representative
	for a in atom_types:
		dgm[a] = [c[a] for c in comps]
	return dgm

def get_neighbour_cycles(points : pandas.DataFrame, cycle : list, filt):
	neighbours = []
	for s in filt:
		for v in s:
			if v in cycle:
				neighbours.append(s)
				break
	return neighbours

def generate_display(points : pandas.DataFrame, dgm : pandas.DataFrame, id : int, filt): #TODO: visualise a neighbourhood of the representative
	"""! Display a representative of a cycle.
	@param points 	pandas.DataFrame of the atoms
	@param dgm 	pandas.DataFrame of the representatives of cycles
	@param id		int corresponding to the id of the cycle you want to visualise

	@return fig		plotly.express figure displaying the ring
	"""
	#print(dgm["idPoint"].loc[id])
	cycle = points.iloc[dgm["idPointz"].loc[id]]
	fig_data = px.scatter_3d(cycle, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",cycle.index]).data
	for i in range(len(cycle)):
		fig_data = fig_data+px.line_3d(cycle.iloc[[i,(i+1)%len(cycle)]],x="x", y="y", z="z").data 
	neighbour_cycles = []
	cycle_set = set(cycle)
	neighbours = get_neighbour_cycles(points, cycle.index.values, filt)
	print(neighbours)
	neighbour_atoms = list(set(v for s in neighbours for v in s))
	print(points.loc[neighbour_atoms])
	neighbour_atoms = points.loc[neighbour_atoms]
	fig_data = fig_data+px.scatter_3d(neighbour_atoms, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",neighbour_atoms.index]).data 
	for s in neighbours:
		s_cycle = list(v for v in s)
		s_cycle = points.loc[s_cycle]
		fig_data = fig_data+go.Figure(go.Mesh3d(x=s_cycle["x"], y=s_cycle["y"],   z=s_cycle["z"],  color="blue",  opacity=.01, alphahull=0)).data
		#for i in range(len(s_cycle)):
		#	fig_data = fig_data+px.line_3d(points.loc[[s_cycle[i],s_cycle[(i+1)%len(s_cycle)]]],x="x", y="y", z="z", fill="toself")).data 
	fig_ring = go.Figure(data=fig_data)
	fig_ring.update_layout( title="Visualisation of a representative of loop with ID: {}".format(id))
	return fig_ring