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

def get_representative_loops(points, atoms, filt, m, d_dgms):
	"""! Get representative of each homology class in dimension 1.

	@param points
	@param atoms
	@param filt
	@param m
	@param d_dgms

	@return dfPD
	"""
	simps_birth = []
	cycle_comps = []
	cycle_reps = []
	for i, c in enumerate(m):
		if i % 1000 == 0:
			print(i)
		simps_birth.append(filt[i].data)
		cycle_comps.append([j for j in filt[i]])
		a = []
		for x in c:
			a = a+[j for j in filt[x.index]]
		a = list(set(a))
	simps_birth = numpy.array(simps_birth)
	pd_cycle = []
	for p in d_dgms[1]:
		# get indices of 2-simplices
		two_cycles = numpy.where(numpy.fromiter(map(len,cycle_comps), dtype="int")==3)[0]
		# indices of simplices born at p.death
		diff = numpy.absolute(p.death-simps_birth) 
		ids_same_death = [i for i in numpy.where(diff == min(diff))[0]]
		#intersec the two sets
		ids_simp_death = numpy.intersect1d(two_cycles, ids_same_death)
		#test if m[ids_simp_death] contains filt[p.data]
		birth_cycl_in = [p.data in [j.index for j in m[i]] for i in ids_simp_death]
		# arbitrarily choose the first one
		pd_cycle.append(cycle_comps[ids_simp_death[birth_cycl_in][0]])
	dfPD = pandas.DataFrame(data={
								"Dimension" : [1 for p in d_dgms[1]],
								"Birth" : [p.birth for p in d_dgms[1]],
								"Death" : [p.death for p in d_dgms[1]],  
								"idPoint" : pd_cycle,
								"Size" : [len(cycle) for cycle in pd_cycle]
	})   
	atom_count = [[] for a in atoms]
	for i in range(dfPD.shape[0]):
		list_atoms_in_cycle = points.iloc[dfPD.iloc[i]["idPoint"]]["Atom"].tolist()
		for i, a in enumerate(atoms):
			atom_count[i].append(list_atoms_in_cycle.count(a))
	#append the count to the dataframe
	print(atom_count)
	for i in range(len(atoms)):
		dfPD[atoms[i]+" count"] = atom_count[i]
	return dfPD

def generate_display(points : pandas.DataFrame, dfPD : pandas.DataFrame, id : int, filt): #TODO: visualise a neighbourhood of the representative
	"""! Display a representative of a cycle.
	@param points 	pandas.DataFrame of the atoms
	@param dfPD 	pandas.DataFrame of the representatives of cycles
	@param id		int corresponding to the id of the cycle you want to visualise

	@return fig		plotly.express figure displaying the ring
	"""
	print(dfPD["idPoint"].iloc[id])
	cycle = points.iloc[dfPD["idPoint"].iloc[id]]
	fig_data = px.scatter_3d(cycle, x="x", y="y", z="z", size="w", color="Atom", hover_data=["Atom",cycle.index]).data
	for i in range(len(cycle)):
		fig_data = fig_data+px.line_3d(cycle.iloc[[i,(i+1)%len(cycle)]],x="x", y="y", z="z").data 
	fig_ring = go.Figure(data=fig_data)
	neighbour_cycles = []
	cycle_set = set(cycle)
	for c in filt:
		print(c)
		c_list=list(c)
		print(c_list)
		if len(c_list) == 3 and len(cycle_set.intersection(set(c_list)))!=0:
			neighbour_cycles.append(c_list)
	print(neighbour_cycles)	
	return fig_ring