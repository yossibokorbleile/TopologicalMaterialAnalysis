from process import *
from plots import *
from visualisation import *
import os
import oineus

file_path = '/Users/jl65ai/Downloads/ZIF_300_wrapped.xyz'
atoms=["H","C","N","Zn"]
radii=[0.389, 0.718, 0.635, 1.491]
repeat_x=repeat_y=repeat_z=1
params=oineus.ReductionParams()
Params.n_threads=16
s=0
atom_locations = sample_at(file_path,"xyz", s, repeat_x, repeat_y, repeat_z, atoms, radii)
dcmp, filt, dgm_1, dgm_2 = oineus_process(atom_locations, params)

