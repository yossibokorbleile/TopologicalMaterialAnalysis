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
		rep = []
		for v in R[idx]:
			rep.append(filt.get_cell(v))
		cycle_reps.append(rep)
	dgm["cycle rep"] = cycle_reps
	return dgm

def get_vertices_and_edges(loop, filt):
	verts = []
	edges = []
	for x in loop:
		verts.append(x.vertices[0])
		verts.append(x.vertices[1])
		edges.append([x.vertices[0], x.vertices[1]])
	return list(set(verts)), edges

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
	dgm["vertices"]=verts
	dgm["edges"]=edges
	return dgm

def get_neighbour_cells(points : pandas.DataFrame, cycle_vertices : list, filt):
	neighbour_vertices = []
	neighbour_edges = []
	for s in filt.simplices():
		if len(s.vertices) == 2:#restrict to edges for the moment
			for v in s.vertices:
				if v in cycle_vertices:
					neighbour_vertices.append(s.vertices[0])
					neighbour_vertices.append(s.vertices[1])
					neighbour_edges.append(s.vertices)
					break
	return neighbour_vertices, neighbour_edges

def generate_display(points : pandas.DataFrame, dgm : pandas.DataFrame, id : int, filt, neighbours = False): #TODO: visualise a neighbourhood of the representative
	"""! Display a representative of a cycle.
	@param points 	pandas.DataFrame of the atoms
	@param dgm 	pandas.DataFrame of the representatives of cycles
	@param id		int corresponding to the id of the cycle you want to visualise

	@return fig		plotly.express figure displaying the ring
	"""
	cycle = points.iloc[[filt.get_id_by_sorted_id(v) for v in dgm["vertices"].loc[id]]]
	fig_data = px.scatter_3d(cycle, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",cycle.index]).data
	for e in dgm["edges"].loc[id]:
		fig_data = fig_data+px.line_3d(points.iloc[[filt.get_id_by_sorted_id(e[0]),filt.get_id_by_sorted_id(e[1])]],x="x", y="y", z="z").update_traces(line_color='red', line_width=5).data 
	if neighbours:
		neighbour_vertices, neighbour_edges = get_neighbour_cells(points, [v for v in dgm["vertices"].loc[id]], filt)
		neighbour_atoms = list(set(filt.get_id_by_sorted_id(v) for v in neighbour_vertices))
		neighbour_atoms = points.loc[neighbour_atoms]
		fig_data = fig_data+px.scatter_3d(neighbour_atoms, x="x", y="y", z="z", size="w", color="Atom").data 
		for e in neighbour_edges:
			fig_data = fig_data+px.line_3d(points.iloc[[filt.get_id_by_sorted_id(e[0]),filt.get_id_by_sorted_id(e[1])]],x="x", y="y", z="z").data 
			#s_cycle = list(v for v in s)
			#s_cycle = points.loc[s_cycle]
			#fig_data = fig_data+go.Figure(go.Mesh3d(x=s_cycle["x"], y=s_cycle["y"],   z=s_cycle["z"],  color="blue",  opacity=.01, alphahull=0)).data
			#for i in range(len(s_cycle)):
			#	fig_data = fig_data+px.line_3d(points.loc[[s_cycle[i],s_cycle[(i+1)%len(s_cycle)]]],x="x", y="y", z="z", fill="toself").data 
	fig_ring = go.Figure(data=fig_data)
	fig_ring.update_layout( title="Visualisation of a representative of loop with ID: {}".format(id))
	return fig_ring