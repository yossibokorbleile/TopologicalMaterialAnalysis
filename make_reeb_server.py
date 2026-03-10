import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing as mp
import pickle
import os

from scipy.spatial.distance import pdist, squareform, cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from itertools import chain, combinations, product

from reeb_graph import Reeb_Graph
from reeb_aux import *


if __name__ == '__main__':

    """
    Setup Data Folder
    """

    data_folder = './Data'
    compositions = ['/67-33','/70-30','/75-25']
    versions = ['/Quench_4','/Quench_5','/Quench_6']
    add_folder = '/MSD'
    data_set = '/md.statsis3'
    axes_names = ['_x_','_y_','_z_']
    results_folder = './Results'
    MP = False

    """
    Make Reeb Graphs
    """

    n_grid = 281 #the complexity scales as n_grid^3: lower it if you cannot run the computations
    fat = 1
    covering = np.array([-1,1])
    reeb_stride = 2
#    relax = [reeb_stride+1,-reeb_stride-2]
    relax = [n_grid//2,-(n_grid//2+1)]
    stride = 20
    t_step = 2
    radius = np.array([0])

    for comp in compositions: 

        print('Doing Glass: ', comp,'                 ')

        for quench_count, vers in enumerate(versions):
            
            print('Doing Version: ', vers,'                    ')

            inputfile = data_folder + comp + vers + add_folder + data_set #preprocessing
            fmean, pLi, radii = estimate_radius(inputfile, t_step)

            print('Radius: ', radii[:10])
            
            for j,V in enumerate(canonical_matrix_grid()): #selecting the axis along which max-flow is computed (and mapper graphs constructed)

                for k in radius: #this was done for robustness, as it perturbs the radii computed at line 55. As you can see in line 44, it is now set to 0. So, no perturbation.

                    if k==np.min(radius):
                        PREV_DATA = (fmean, pLi, None, None)
                    else:
                        PREV_DATA = (fmean, pLi, reeb.distances_to_balls, reeb.atoms_kinds) 

                    print('Doing Radius: ', comp, vers, j, k, '                    ')

                    r = radii + 0.1*k

                    reeb = Reeb_Graph(inputfile, radii = r,
                            grid_size = n_grid, 
                            t_step = t_step,
                            PREV_DATA = PREV_DATA,
                            periodic = True,
                            fat_radius = fat,
                            covering = covering,
                            reeb_stride = reeb_stride,
                            transform_points = V,
                            relax_z_axis = relax,
                            verbose = True, save_RAM = True, stride=stride, MP=False)

        
                    reeb.make_reeb_graph(plot=False)

                    print('Reeb Graph Done!')

#                    flow = circular_max_flow(reeb)                    
#                    reeb.flow = flow

#                    print('Flow Computed ', flow)

                    name = results_folder + comp + '_' + vers[1:] + '_reeb' + axes_names[j] + str(k) + '.pickle' 

                    with open(name, "wb") as fp:
                        pickle.dump(reeb, fp)

                    name_ = comp[1:] + '_' + vers[1:] + '_reeb' + axes_names[j] + str(k) + '.pickle' 
                    zip_name = comp[1:] + '_' + vers[1:] + '_reeb' + axes_names[j] + str(k) + '.zip'      

                    cmd = 'cd ' + results_folder + ' && zip '+ zip_name + ' ' + name_
                    os.system(cmd)

                    cmd = 'rm '+ name 
                    os.system(cmd)
                    


    
















