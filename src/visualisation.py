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

def get_representative_loops(dgm : pandas.DataFrame, R, filt):
	"""! Get representative of each homology class in dimension 1.
	@param dgm 		pandas.DataFrame containing the diagram and the birth and death simplex id
	@param filt		oineus filtration
	
	@return dgm		pandas.DataFrame with new column `cycle rep`
	"""
	cycle_reps = []
	for i in range(dgm.shape[0]):
		idx = dgm["death simplex"].iloc[i]
		#print(idx)
		sorted_rep = R[idx]
		rep = []
		for v in sorted_rep:
			rep.append(filt.get_cell(v))
		cycle_reps.append(rep)
		#print("{} ({}, {}): unsorted birth_id {} sorted birth id {} rep sorted {} rep unsorted {}".format(i,dgm["birth"].iloc[i], dgm["death"].iloc[i],idx,filt.get_sorted_id_by_id(idx),R[filt.get_sorted_id_by_id(idx)], rep))
	dgm["cycle rep"] = cycle_reps
	return dgm

def get_vertices_and_edges(loop, filt):
	verts = []
	edges = []
	#print("gettin vertices and edges")
	for x in loop:
		#print(x)
		verts.append(x.vertices[0])
		verts.append(x.vertices[1])
		edges.append([x.vertices[0], x.vertices[1]])
	#print("got them")
	return verts, edges

def loop_composition(verts, filt, points, atom_types):
	"""! Get the composition of a given representative 

	@param loop		list containing the ids of the edges in the loop
	@param points	pandas.DataFrame of the points (atoms) in the structure
	@param atom_types list of the atom types we are considering (restricting to)

	@return comp	dictionary listing the number of atoms of each type in the loop
	"""
	comp = dict([(a, 0) for a in atom_types])
	for v in verts:
		comp[points["Atom"].iloc[filt.get_id_by_sorted_id(v)]] += 1
	return comp

def generate_visulisation_df(dgm : pandas.DataFrame, R, filt, points, atom_types):  
	"""! generate the pandas.DataFrame containing the information about the points so we can display it 

	@param dgm		pandas.DataFrame of the diagram, with columns `birth`, `death`, `birth simplex`, `death simplex`, `cycle rep`
	"""
	dgm =  get_representative_loops(dgm, R, filt)
	#print("next")
	#get the composition of the cycle representatives
	comps = []
	edges = []
	verts = []
	for i in range(dgm.shape[0]):
		verts_i, edges_i = get_vertices_and_edges(dgm["cycle rep"].iloc[i], filt)
		#print("here?")
		edges.append(edges_i)
		verts.append(verts_i)
		comps.append(loop_composition(verts_i, filt, points, atom_types))
	#for each atom type we are looking at, add a column with the number of atoms of this type in the cycle representative
	for a in atom_types:
		dgm[a] = [c[a] for c in comps]
	dgm["edges"]=edges
	dgm["vertices"]=verts
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
	cycle = points.iloc[[filt.get_id_by_sorted_id(v) for v in dgm["vertices"].loc[id]]]
	fig_data = px.scatter_3d(cycle, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",cycle.index]).data
	print(dgm["edges"].loc[id])
	for e in dgm["edges"].loc[id]:
		fig_data = fig_data+px.line_3d(points.iloc[[filt.get_id_by_sorted_id(e[0]),filt.get_id_by_sorted_id(e[1])]],x="x", y="y", z="z").data 
	#neighbour_cycles = []
	#cycle_set = set(cycle)
	#neighbours = get_neighbour_cycles(points, cycle.index.values, filt)
	#print(neighbours)
	#neighbour_atoms = list(set(v for s in neighbours for v in s))
	#print(points.loc[neighbour_atoms])
	#neighbour_atoms = points.loc[neighbour_atoms]
	#fig_data = fig_data+px.scatter_3d(neighbour_atoms, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",neighbour_atoms.index]).data 
	#for s in neighbours:
		# s_cycle = list(v for v in s)
		# s_cycle = points.loc[s_cycle]
		# fig_data = fig_data+go.Figure(go.Mesh3d(x=s_cycle["x"], y=s_cycle["y"],   z=s_cycle["z"],  color="blue",  opacity=.01, alphahull=0)).data
		#for i in range(len(s_cycle)):
		#	fig_data = fig_data+px.line_3d(points.loc[[s_cycle[i],s_cycle[(i+1)%len(s_cycle)]]],x="x", y="y", z="z", fill="toself")).data 
	fig_ring = go.Figure(data=fig_data)
	fig_ring.update_layout( title="Visualisation of a representative of loop with ID: {}".format(id))
	return fig_ring