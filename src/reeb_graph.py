##
# @file reeb_graph.py
# @brief Reeb graph construction and analysis for void region connectivity in materials.
# @version 1.3.0
# @date March 2026

import numpy as np
import matplotlib.pyplot as plt
#import galois
import networkx as nx
import multiprocessing as mp

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.spatial.distance import pdist, squareform, cdist
from scipy.sparse import csr_matrix, vstack
from scipy.sparse.csgraph import dijkstra, connected_components, bellman_ford, maximum_flow, shortest_path
from scipy.sparse.linalg import lsqr

from itertools import chain, combinations, product

from diffusion_utils import *
from reeb_aux import *

class Reeb_Graph(object):
    
    """! @brief Reeb graph class for analysing void region connectivity in materials.

    The Reeb_Graph class constructs a Reeb graph from the void regions of a material
    structure, enabling the study of transport pathways and connectivity. It discretises
    the void space onto a grid, constructs a graph from connected components of level sets,
    and provides methods for maximum flow computation and tunnel identification.

    @section reeb_usage Usage
    @code{.py}
    from reeb_graph import Reeb_Graph
    import numpy as np

    reeb = Reeb_Graph(
        backbone=atom_positions,
        flow=flow_positions,
        radii=radii_array,
        M=box_max, m=box_min,
        grid_size=50,
        periodic=True
    )
    reeb.make_reeb_graph(axis=2, old=False, plot=True)
    max_flow = reeb.compute_max_flow()
    @endcode
    """
    
    def __init__(self, inputfile = None, backbone = None, flow = None, radii = None, M = None, m = None,
                                grid_size = 50, 
                                t_step = 50,
                                PREV_DATA = (None,None,None,None),
                                periodic = True,
                                fat_radius = 1,
                                covering = np.array([-1,1]),
                                reeb_stride = 1,
                                swap_res_grid_and_balls = False,
                                transform_points = None,
                                relax_z_axis = None,
                                verbose = False,
                                MP = False,
                                save_RAM = True,
                                stride = 1):
        """! @brief Initialise the Reeb_Graph with material structure data and grid parameters.

        @param inputfile Optional path to an input file containing atomic structure data.
        @param backbone Numpy array of shape (N, 3) with backbone atom positions.
        @param flow Numpy array of shape (N, 3) with flow particle positions.
        @param radii Vector or scalar determining the radius around backbone atoms, defining the void region.
        @param M Numpy array of shape (3,) with the upper bounds of the simulation box.
        @param m Numpy array of shape (3,) with the lower bounds of the simulation box.
        @param grid_size Number of grid points along each axis (default 50).
        @param t_step Time step parameter (default 50).
        @param PREV_DATA Tuple (backbone, flow, distances, atom_kinds) from previous calculations to reuse distances.
        @param periodic Whether to apply periodic boundary conditions (default True).
        @param fat_radius Multiplier for the cell diameter to define the graph neighbourhood radius (default 1).
        @param covering Array defining the level set thickness, e.g. [-1,1] for symmetric covering (default [-1,1]).
        @param reeb_stride Stride along the axis for Reeb graph construction (default 1).
        @param swap_res_grid_and_balls If True, swap forbidden region and residual grid (default False).
        @param transform_points Affine transformation (V, v) to apply to atom positions as V*pts + v.
        @param relax_z_axis Array [a, -b] to extend the z-axis by placing boundary layers (default None).
        @param verbose If True, print progress and diagnostic information (default False).
        @param MP If True, use multiprocessing for distance calculations (default False).
        @param save_RAM If True, use memory-efficient graph construction (default True).
        @param stride Stride for subdividing distance calculations between backbone and grid (default 1).
        """   
    
    
        """
        Parameters Setup
        """
        # self.inputfile = inputfile        
        self.radii = radii
        self.t_step = t_step
        self.M = M
        self.m = m
        self.grid_size = grid_size
        self.axes = np.array([0,1])
        self.dim = 3
        # self.backbone, self.flow, self.distances_to_balls_aux, self.atoms_kinds_aux = PREV_DATA
        self.backbone = backbone
        self.flow = flow
        self.distances_to_balls_aux = None
        self.atoms_kinds_aux = None
        self.periodic = periodic    
        self.fat_radius = fat_radius
        self.reeb_stride = reeb_stride
        self.covering = covering
        self.swap_res_grid_and_balls = swap_res_grid_and_balls
        self.transform_points = transform_points
        
        if self.periodic:
            self.simply_connected_top_bottom = False
            self.relax_z_axis = np.array([self.covering[-1]+1,-1])
        else:
            self.simply_connected_top_bottom = True
            self.relax_z_axis = relax_z_axis
        
        self.verbose = verbose
        self.save_RAM = save_RAM
        self.stride = stride
        self.MP = MP
        
        if self.verbose and self.periodic:
            print('This number should be ONE for accurate periodic graphs!', self.grid_size%self.reeb_stride)
            print('Coverings should be symmetric for accurate periodic graphs!', np.min(self.covering),np.max(self.covering))
            print('Coverings and stride should be coherent: stride = covering + 1 ! Stride: ', 
                  self.reeb_stride, 'Covering Step: ', np.max(self.covering))
        self.setup_grid(balls_centres = self.backbone, balls_radii = self.radii, M = self.M, m = self.m)
           
        
    """
    Auxiliary Functions
    """    
    def setup_grid(self, balls_centres = None, balls_radii = None,
                                M = None, m = None):
        """! @brief Set up the computational grid and residual region.

        Optionally applies an affine transformation to backbone positions, then
        generates the residual grid by excluding forbidden regions around atoms.

        @param balls_centres Numpy array of shape (N, 3) with atom centre positions.
        @param balls_radii Numpy array of atom radii defining forbidden regions.
        @param M Upper bounds of the simulation box.
        @param m Lower bounds of the simulation box.
        """
        if not self.transform_points is None:
            balls_centres, balls_radii = self.transform_backbone(balls_centres = balls_centres ,balls_radii = balls_radii,
                                M = M, m = m)

        self.produce_res_grid(balls_centres = balls_centres ,balls_radii = balls_radii,
                                        M = M, m = m)  
                 

    def transform_backbone(self, balls_centres = None, balls_radii = None,
                                M = None, m = None):
        """! @brief Apply affine transformation to backbone atom positions.

        Extracts the transformation matrix V and translation vector v from
        self.transform_points and applies them to atom centres.

        @param balls_centres Numpy array of shape (N, 3) with atom positions.
        @param balls_radii Numpy array of atom radii.
        @param M Upper bounds of the simulation box.
        @param m Lower bounds of the simulation box.
        @return Tuple (balls_centres, balls_radii) with transformed positions and filtered radii.
        """
        a=0        
        try:
            self.transform_points.shape
        except:
            a=1
        
        if a==1:
            V,v = self.transform_points
        else:
            V = self.transform_points
            v = 0

        if not self.inputfile is None: 
            
            P = balls_centres[:self.n_P,:]
            r_P = balls_radii[:self.n_P]
            S = balls_centres[self.n_P:,:]
            r_S = balls_radii[self.n_P:]

            P, r_P = self.affine_transform(P, r_P, V, v, M, m)
            S, r_S = self.affine_transform(S, r_S, V, v, M, m)

            balls_centres = np.concatenate([P, S])    
            balls_radii = np.concatenate([r_P,r_S])
            
        else:
            balls_centres, balls_radii = self.affine_transform(balls_centres, balls_radii, V, v, M, m)
       
        return balls_centres, balls_radii
                
    
    def affine_transform(self, points, radii, V, v, M, m):
        """! @brief Apply affine coordinate transformation with periodic boundary handling.

        Replicates points across periodic boundaries, applies the affine
        transformation V*points + v, and filters to retain only those within
        the simulation box [m, M].

        @param points Numpy array of shape (N, 3) with point coordinates.
        @param radii Numpy array of radii associated with each point.
        @param V Transformation matrix.
        @param v Translation vector.
        @param M Upper bounds of the simulation box.
        @param m Lower bounds of the simulation box.
        @return Tuple (points, radii) after transformation and filtering.
        """
        DELTAS = M-m
        vectors = list([[0,-DELTAS[i],DELTAS[i]] for i in range(len(DELTAS))])
        
        points_tmp = np.copy(points)
        radii_tmp = np.copy(radii)

        for v_ in set(product(*vectors)):
            
            if np.max(np.abs(v_))>0:
                points = np.concatenate([points, points_tmp + v_])
                radii = np.concatenate([radii, radii_tmp])
            
        points = points@V+v
        tmp_u = np.max(points-M,axis=-1)<0
        tmp_l = np.min(points-m,axis=-1)>0
        tmp = np.multiply(tmp_u,tmp_l)

        points = points[tmp]
        radii = radii[tmp]

        return points, radii
           
        
    def produce_res_grid(self, balls_centres = None ,balls_radii = None,
                                M = None, m = None):
        """! @brief Generate the residual grid (void region) by excluding forbidden regions around atoms.

        Creates a 3D grid over the simulation box, computes distances from each grid
        point to the nearest atom, and marks points outside all atomic radii as the
        residual (void) region. Also builds a sparse adjacency graph over the residual grid.

        @param balls_centres Numpy array of shape (N, 3) with atom centre positions.
        @param balls_radii Numpy array of atom radii defining forbidden regions.
        @param M Upper bounds of the simulation box.
        @param m Lower bounds of the simulation box.
        """
        if self.verbose:
            print('The filtration radii are: ', balls_radii[0],
                                                balls_radii[-1])
          
        tmp = 1

        self.nx, self.ny, self.nz = (self.grid_size, self.grid_size, self.grid_size)
        
        self.x = np.linspace(self.m[0], self.M[0], self.nx+tmp)[:-1]
        self.y = np.linspace(self.m[1], self.M[1], self.ny+tmp)[:-1]
        self.z = np.linspace(self.m[2], self.M[2], self.nz+tmp)[:-1]
            
        self.d_x = self.x[1]-self.x[0]
        self.d_y = self.y[1]-self.y[0]
        self.d_z = self.z[1]-self.z[0]
        
        self.r_graph = np.sqrt(self.d_x**2 + self.d_y**2 + self.d_z**2)*1.001*self.fat_radius   

        if not self.relax_z_axis is None: 
            
            idx_l,idx_u = self.relax_z_axis
            
            idxs_low = np.sum(self.z<self.z[idx_l])
            idxs_up = np.sum(self.z>self.z[idx_u])
            
            up = self.z[-1]+self.d_z*(idxs_low)
            low = self.z[0]-self.d_z*idxs_up
            
            z = np.linspace(low, up, self.nz+np.abs(idx_u)-1+np.abs(idx_l))
            self.z = np.arange(low,self.z[-1]+self.d_z*(idxs_low+1),self.d_z)
            
            if len(z)<len(self.z):
                self.z = z
                self.d_z = self.z[1]-self.z[0]
            
            self.nz = len(self.z)
            
        xv, yv, zv = np.meshgrid(self.x, self.y, self.z)

        self.grid = np.zeros((self.nx*self.ny*self.nz,self.dim))
        self.N = len(self.grid)
        
        self.grid[:,0] = xv.flatten()
        self.grid[:,1] = yv.flatten()
        self.grid[:,2] = zv.flatten()
        
        self.unit_3d = np.prod(self.M-self.m)/self.N
        self.unit_2d = np.prod(self.M[:-1]-self.m[:-1])/(self.nx*self.ny)                  
                        
        k = (1-0.99)
        z = np.unique(self.grid[:,-1])
        self.z_0 = z[0]+self.d_z*k
        self.z_1 = z[-1]-self.d_z*k
        axes_aux = np.array([0,1,2])
        
        if self.simply_connected_top_bottom:   
            aux_top = self.grid[:,-1]>=self.z_1
            aux_bottom = self.grid[:,-1]<=self.z_0
        else:
            aux_top = np.zeros_like(self.grid.shape[0])
            aux_bottom = np.zeros_like(self.grid.shape[0])

        if self.distances_to_balls_aux is None:
            grid_aux = self.m + np.remainder(self.grid-self.m, self.M-self.m) 
            self.distances_to_balls = np.inf*np.ones_like(self.grid[:,0]) 
            self.atoms_kinds = -1*np.ones_like(self.grid[:,0]).astype(int) 
        else:
            self.distances_to_balls = self.distances_to_balls_aux  
            self.atoms_kinds = self.atoms_kinds_aux.astype(int)
            
        if self.distances_to_balls_aux is None:
            print("saving RAM is: ", self.save_RAM)
            print("MP is: ", self.MP)
            if self.save_RAM:

                idxs_grid = np.arange(0,len(balls_centres),self.stride)
                
                if self.MP:
                    pool = mp.Pool(processes=4)

                    D = pool.map(dist_from_pts_periodic_boundaries_pool,
                                     ([balls_centres[i:i+self.stride,:],grid_aux,self.M,self.m,axes_aux,self.dim] 
                                             for i in idxs_grid))
                    pool.close()
                    
                    for i,d in enumerate(D):
                        self.distances_to_balls = np.vstack([self.distances_to_balls,d])
                        r_aux = np.argmin(self.distances_to_balls,axis=0).astype(int)
                        self.distances_to_balls = np.min(self.distances_to_balls,axis=0)

                        self.atoms_kinds[r_aux>0] = r_aux[r_aux>0]-1+idxs_grid[i] 
                        
                else:
                    for i in idxs_grid:
                        if self.verbose and (i%1000)==0:
                            print('Doing atom ',i,len(balls_centres), end='\r')
                        d = dist_from_pts_periodic_boundaries_numba(balls_centres[i:i+self.stride,:], grid_aux,self.M,self.m,axes_aux)#safe_dist_from_pts_periodic_boundaries_numba
 
                        self.distances_to_balls = np.vstack([self.distances_to_balls,d])
                        r_aux = np.argmin(self.distances_to_balls,axis=0).astype(int)
                        self.distances_to_balls = np.min(self.distances_to_balls,axis=0)
                        # print("distances to balls: ", self.distances_to_balls, " with shape ", self.distances_to_balls.shape)
                        # print("r_aux: ", r_aux, " with shape ", r_aux.shape)
                        # print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)    
                        self.atoms_kinds[r_aux>0] = r_aux[r_aux>0]-1+i 
            else:
                D = dist_from_pts_periodic_boundaries_numba(balls_centres,grid_aux,self.M,self.m,axes_aux)#safe_dist_from_pts_periodic_boundaries_numba
                self.atoms_kinds = np.argmin(D,axis=0)
                self.distances_to_balls = np.min(D,axis=0)
                
        # print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)
        # print("balls radii: ", balls_radii, " with shape ", balls_radii.shape)
        # print("distances to balls: ", self.distances_to_balls, " with shape ", self.distances_to_balls.shape)
        print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)
        print("distances to balls: ", self.distances_to_balls)
        thresh_radii = [balls_radii[i] for i in self.atoms_kinds]
        thresh_aux = np.vstack([self.distances_to_balls, thresh_radii])
        idxs = np.argmin(thresh_aux, axis=0)==1  
        
        self.idxs = (idxs + aux_top + aux_bottom)>0     
        
        if self.swap_res_grid_and_balls:
            if self.simply_connected_top_bottom:
                self.idxs = ((1-self.idxs) + aux_top + aux_bottom)>0
            else:
                self.idxs = (1-self.idxs)
        
            r = self.d_x*0
            self.balls_centres = self.grid[self.idxs]  
            self.balls_radii = np.ones_like(self.idxs)*r
        else:
            self.balls_centres = balls_centres
            self.balls_radii = balls_radii
            
        N=sum(self.idxs)  
        self.res_grid=self.grid[self.idxs]          
            
        if self.verbose:
            print('The Residual Grid has: ', N,' points.')
            print('Building Graph for the Point Cloud.')

        if self.save_RAM:
            
            res_to_grid = np.where(self.idxs>0)[0]
            grid_to_res = np.cumsum(self.idxs,dtype=int)-1

            grid_to_cubic = np.unravel_index(np.arange(len(self.grid)),(self.nx,self.ny,self.nz))
            grid_to_cubic = np.array(grid_to_cubic).T[:,[1,0,2]]

            xv, yv, zv = np.meshgrid(np.arange(self.nx), np.arange(self.ny), np.arange(self.nz))

            cubic_to_grid = np.ravel_multi_index((xv,yv,zv),(self.nx,self.ny,self.nz))

            if self.verbose:
                        print('Building Graph with Numba.')

            A, cnt = make_graph_fast(self.grid,self.idxs,
                                     res_to_grid,grid_to_res,
                                     grid_to_cubic, cubic_to_grid,
                                     self.nx,self.ny,self.nz,
                                     self.dim,self.M,self.m)
            
            A = A[:cnt,:]
            self.graph=csr_matrix((np.ones_like(A[:,0]), (A[:,0],A[:,1])), shape=(N, N)) 
                
        else:
            D = dist_periodic_boundaries(self.res_grid, self.M, self.m, self.axes, self.dim)        
            self.graph = csr_matrix(D<self.r_graph)
            
        if self.verbose:
            print('Grid Finished.')


    def make_reeb_graph(self, axis=-1, old=False, plot = False):
        """! @brief Construct the Reeb graph from the residual grid.

        Computes the Reeb graph along the specified axis using either the current
        or legacy algorithm, builds a NetworkX graph, and optionally plots the result.

        @param axis The coordinate axis along which to slice the residual grid (default -1, i.e. z-axis).
        @param old If True, use the legacy Reeb graph algorithm (default False).
        @param plot If True, display a 3D plot of the Reeb graph after construction (default False).
        """
        grid = self.res_grid

        if old:
            E, weights, fun, emb = self.reeb_graph_old(grid, axis=axis, D_ = self.graph) 
        else:
            E, weights, fun, emb = self.reeb_graph(grid, axis=axis, D_ = self.graph) 
                   
        self.G = self.build_graph(E, weights, fun, emb)       
        self.measure_path_comp = fun
        
        if plot:
            self.plot_reeb()

               
    def reeb_graph(self, grid, axis=-1, D_ = None):
        """! @brief Core Reeb graph computation using connected components of level sets.

        Sweeps through the grid along the given axis, computing connected components
        at each level set and connecting components between consecutive levels based
        on their intersection.

        @param grid Numpy array of shape (N, 3) with residual grid points.
        @param axis The coordinate axis for the sweep (default -1, i.e. z-axis).
        @param D_ Precomputed sparse adjacency matrix for the grid (default None, computed internally).
        @return Tuple (E, weights, fun, emb) containing edges, edge weight matrix,
                node measures, and node embedding coordinates.
        """
        self.axis = axis
        g = lambda x: len(x)
        
        f = grid[:,axis]
        values = np.unique(f)
        idxs = np.arange(0,len(f))
        l,u = self.covering
        
        if D_ is None:
            D_graph = safe_dist_from_pts_periodic_boundaries_numba(grid,grid,self.M,self.m,self.axes, self.dim)
            D_graph = (D_graph<self.r_graph)
        else:
            D_graph = D_
                
#        delta = make_neigh(grid, D_graph.indptr, D_graph.indices, f)
        delta = [self.d_x,self.d_y,self.d_z][axis]

        """
        delta -> di quanto cambia la funzione in un intorno di un punto: devo essere sicuro che tutti i punti di un intorno 
                 di un punto di un aperto del ricoprimento, stiamo o nell'aperto precedente, o nel successivo. 
                 Quindi delta e' il raggio degli aperti del ricoprimento.
        norm ->  il diametro di un aperto del ricoprimento é DELTA*sum(covering). Calcolando quante volte ci sta un passo della
                 griglia lungo z in un aperto del ricoprimento ottengo da quanti "PIANI" é composta la "TORTA" dell intersez di 
                 f^-1. 
                 Devo aggiungere 1 perché il piano 0 conta. Quindi prendo il numero medio di punti per piano. Con minimo 1.
        """
        
        self.delta = 1.001*delta        
        self.norm = ((self.delta*np.max(self.covering))//(self.d_z*self.reeb_stride)) + 1
        
        count = 0
        
        """
        Tengo conto dello stride
        """
        reeb_size = ((values.shape[0]-1)//self.reeb_stride)*self.reeb_stride
        reeb_idxs = np.arange(0,reeb_size+1,self.reeb_stride)        
        reeb_idxs = np.concatenate([reeb_idxs, 
                                    np.arange(reeb_idxs[-1]+1,values.shape[0],1)])
        reeb_size = len(reeb_idxs)
        
        """
        Faccio il reeb
        """
        reeb_values = values[reeb_idxs]
        
        self.REEB_GRAPH = -1*np.ones((grid.shape[0],reeb_size)).astype(int)
        fun = []
        emb = []
        weights = np.array([[0,0,0]], dtype=int)
        
        fat_lvl = (f<=values[0]+u*self.delta)
        lvl_to_grid = np.where(fat_lvl>0)[0]

        n_c,C = connected_components(D_graph[:,fat_lvl][fat_lvl,:])

        for c in range(n_c):
            tmp_ = np.where(C==c)[0]
            tmp = lvl_to_grid[tmp_]
            self.REEB_GRAPH[tmp,0] = count 
            fun.append(np.max([np.rint(g(tmp_)/self.norm),1]))
            emb.append(np.mean(grid[fat_lvl],axis=0))
            count += 1

        fun = np.array(fun)
        emb = np.array(emb)
    
        for j,v in enumerate(reeb_values[1:]):
            
            i=j+1
            
            sub_lvl = (f<=v+u*self.delta)
            sup_lvl = (f>=v+l*self.delta)
            fat_lvl = (sub_lvl*sup_lvl)

            lvl_to_grid = np.where(fat_lvl>0)[0]

            n_c,C = connected_components(D_graph[:,fat_lvl][fat_lvl,:])

            self.REEB_GRAPH[:,i], AUX_w, AUX_f, AUX_e, count = connect_lvl_sets_aux(i, n_c, C, grid, lvl_to_grid, self.REEB_GRAPH, self.norm, count)

            weights= np.concatenate([weights, AUX_w],axis=0)
            fun = np.hstack([fun, AUX_f]) 
            emb = np.vstack([emb, AUX_e]) 

        weights = csr_matrix((weights[0:,-1],(weights[0:,0].astype(int),weights[0:,1].astype(int))), shape=(len(fun),len(fun)) )
        E = np.argwhere(weights>0)
        self.n = len(emb)-1
            
        self.E = E
            
#            pts_bottom = np.unique(self.REEB_GRAPH[:,0])
#            pts_bottom = pts_bottom[pts_bottom>-1]
#            pts_top = np.unique(self.REEB_GRAPH[:,-1])
#            pts_top = pts_top[pts_top>-1]
            
#            weights[pts_top,pts_bottom] = np.max(weights) #assuming corresponding components are returned with the same index
            
#            E_aux = np.zeros((len(pts_top),2),dtype=int)
#            E_aux[:,0] = pts_top
#            E_aux[:,1] = pts_bottom
            
#            E = np.concatenate([E,E_aux],axis=0)
            
        return E, weights, fun, emb
    
    def reeb_graph_old(self, grid, axis=-1, D_ = None):
        """! @brief Legacy Reeb graph computation method.

        An older implementation of the Reeb graph algorithm that uses dictionary-based
        edge weights and a list-based edge collection. Retained for backward compatibility.

        @param grid Numpy array of shape (N, 3) with residual grid points.
        @param axis The coordinate axis for the sweep (default -1, i.e. z-axis).
        @param D_ Precomputed sparse adjacency matrix for the grid (default None, computed internally).
        @return Tuple (E, weights, fun, emb) containing edge list, weight dictionary,
                node measures, and node embedding coordinates.
        """
        self.axis = axis
        g = lambda x: len(x)
        
        f = grid[:,axis]
        values = np.unique(f)
        idxs = np.arange(0,len(f))
        l,u = self.covering
        
        if D_ is None:
            D_graph = safe_dist_from_pts_periodic_boundaries_numba(grid,grid,self.M,self.m,self.axes, self.dim)
            D_graph = (D_graph<self.r_graph)
        else:
            D_graph = D_
                
#        delta = make_neigh(grid, D_graph.indptr, D_graph.indices, f)
        delta = [self.d_x,self.d_y,self.d_z][axis]
        """
        delta -> di quanto cambia la funzione in un intorno di un punto: devo essere sicuro che tutti i punti di un intorno 
                 di un punto di un aperto del ricoprimento, stiamo o nell'aperto precedente, o nel successivo. 
                 Quindi delta e' il raggio degli aperti del ricoprimento.
        norm ->  il diametro di un aperto del ricoprimento é DELTA*sum(covering). Calcolando quante volte ci sta un passo della
                 griglia lungo z in un aperto del ricoprimento ottengo da quanti "PIANI" é composta la "TORTA" dell intersez di 
                 f^-1. 
                 Devo aggiungere 1 perché il piano 0 conta. Quindi prendo il numero medio di punti per piano. Con minimo 1.
        """
        
        self.delta = 1.001*delta        
        self.norm = ((self.delta*np.max(self.covering))//(self.d_z*self.reeb_stride)) + 1
        
        count = 0
        
        """
        Tengo conto dello stride
        """
        reeb_size = ((values.shape[0]-1)//self.reeb_stride)*self.reeb_stride
        reeb_idxs = np.arange(0,reeb_size+1,self.reeb_stride)        
        reeb_idxs = np.concatenate([reeb_idxs, 
                                    np.arange(reeb_idxs[-1]+1,values.shape[0],1)])
        reeb_size = len(reeb_idxs)
        
        """
        Faccio il reeb
        """
        reeb_values = values[reeb_idxs]
        
        self.REEB_GRAPH = -1*np.ones((grid.shape[0],reeb_size)).astype(int)
        E = []
        fun = []
        emb = []
        weights = {}
        
        fat_lvl = (f<=values[0]+u*self.delta)
        lvl_to_grid = np.where(fat_lvl>0)[0]

        n_c,C = connected_components(D_graph[:,fat_lvl][fat_lvl,:])

        for c in range(n_c):
            tmp_ = np.where(C==c)[0]
            tmp = lvl_to_grid[tmp_]
            self.REEB_GRAPH[tmp,0] = count 
            fun.append(np.max([np.rint(g(tmp_)/self.norm),1]))
            emb.append(np.mean(grid[fat_lvl],axis=0))
            count += 1

        for j,v in enumerate(reeb_values[1:]):
            
            i=j+1
            
            sub_lvl = (f<=v+u*self.delta)
            sup_lvl = (f>=v+l*self.delta)
            fat_lvl = (sub_lvl*sup_lvl)

            lvl_to_grid = np.where(fat_lvl>0)[0]

            n_c,C = connected_components(D_graph[:,fat_lvl][fat_lvl,:])

            for c in range(n_c):
                tmp_ = np.where(C==c)[0]
                tmp = lvl_to_grid[tmp_]

                self.REEB_GRAPH[tmp,i] = count 
                idxs_aux = np.unique(self.REEB_GRAPH[tmp,i-1])        
                idxs_aux=idxs_aux[idxs_aux>-1]

                for old_c in idxs_aux:
                    now = (self.REEB_GRAPH[:,i] == count)
                    old = (self.REEB_GRAPH[:,i-1] == old_c)
                    inters = np.sum(np.multiply(now,old))
                
                    if inters>0:
                        E.append((old_c,count))
                        weights[(old_c,count)]= np.max([inters//self.norm,1])

                fun.append(np.max([np.rint(g(tmp_)/self.norm),1]))
                emb.append(np.mean(grid[tmp],axis=0))

                count += 1

        fun = np.array(fun)
        emb = np.array(emb)
        self.n = len(emb)-1
        
        return E, weights, fun, emb

    
    def build_graph(self,E,weights,fun,emb):
        """! @brief Build a NetworkX graph from Reeb graph edges, weights, and embeddings.

        Creates an undirected graph with node positions set from the embedding coordinates
        and edge capacities set from the weight matrix. Also computes the reachability
        mask between source (node 0) and sink (last node).

        @param E Edge array or list of (source, target) pairs.
        @param weights Edge weight matrix (sparse or dict) representing intersection sizes.
        @param fun Numpy array of node measures (normalised component sizes).
        @param emb Numpy array of shape (N, 3) with node embedding positions.
        @return NetworkX Graph with node positions and edge capacities.
        """
        # create networkx graph
        G=nx.Graph()
        
        fun_ = {}
        
        # add vertices
        for p in range(len(emb)):
            G.add_node(p, pos=(emb[p][0],emb[p][1],emb[p][-1]))
            fun_[p] = fun[p]
            
        nx.set_node_attributes(G,fun_,'measure')


        if type(weights)==type({}):

            self.D_reeb_w = np.zeros((len(emb),len(emb)))
                
            # add edges
            for edge in E:
                G.add_edge(edge[0], edge[1])
                self.D_reeb_w[edge[0], edge[1]] = weights[edge[0], edge[1]]

            self.D_reeb_w = csr_matrix(self.D_reeb_w)
            nx.set_edge_attributes(G,weights,'capacity')
        else:
            self.D_reeb_w = csr_matrix(weights)
            capacities = {}

            for edge in E:
                G.add_edge(edge[0], edge[1])
                capacities[tuple(edge)] = weights[edge[0], edge[1]]

            nx.set_edge_attributes(G,capacities,'capacity')
            

        N = len(emb)-1
        D_0 = shortest_path(self.D_reeb_w, directed=True, indices=[0], return_predecessors=False)[0]<np.inf                
        D_1 = shortest_path(self.D_reeb_w.T, directed=True, indices=[N], return_predecessors=False)[0]<np.inf 
                
        self.reeb_connected = np.multiply(D_0,D_1)
 
        return G
    
    
    def compute_max_flow(self,D_ = None):
        """! @brief Compute the maximum flow through the Reeb graph from source to sink.

        Uses scipy's maximum_flow algorithm on the weighted adjacency matrix to
        find the maximum flow from node 0 (source) to the last node (sink).
        Results are stored in self.max_flow and self.D_flow.

        @param D_ Optional sparse weight matrix to use instead of self.D_reeb_w.
        """
        if D_ is None:
            D = self.D_reeb_w.astype(int)
        else:
            D = np.copy(D_.astype(int))
 
        graph = csr_matrix(D)
        N = D.shape[0]-1
        
        if N == 0:
            D = np.zeros((2,2), dtype=int)
            graph = csr_matrix(D)
            N=1
                        
        self.max_flow = maximum_flow(graph,0,N)
        self.D_flow = self.max_flow.flow
        mask = (self.max_flow.flow>0).astype(int)
        self.D_flow = self.D_flow.multiply(mask)
        
    
    def compute_tunnels(self, n=10, D_= None, independent = False, plot=False):
        """! @brief Find transport tunnels through the Reeb graph using shortest path algorithms.

        Iteratively finds paths from source (node 0) to sink (last node) using the
        Bellman-Ford algorithm on negated weights, extracting the widest paths first.
        After each path is found, edge capacities are reduced (or removed if independent).

        @param n Maximum number of tunnels to find (default 10).
        @param D_ Optional sparse weight matrix to use instead of self.D_reeb_w.
        @param independent If True, remove edges entirely after use; if False, subtract the
               minimum capacity along the path (default False).
        @param plot If True, plot each tunnel as it is found (default False).
        """
        dummy = 1-independent
        self.PATHS = []
        self.MIN = [] 
        self.VOL = []
        
        if D_ is None:
            M = np.max(self.D_reeb_w)
            D = csr_matrix(self.D_reeb_w/M)
        else:
            M = 1
            D = csr_matrix(D_)
            
        N = D.shape[0]-1
        mask = (D>0).astype(int)
        
        for i in range(n):
                    
            graph = (-D).multiply(mask)
            d, pred = bellman_ford(graph, directed=True, indices=[0], return_predecessors=True)    
            
            d = d[0][N]
            pred = pred[0]
                            
            if d < np.inf:
                idx = N
                path = [idx]
                aux = []

                while idx!=0:
                    idx = pred[idx]
                    path.append(idx)
                    aux.append(D[path[-1],path[-2]])

                path = np.array(path)[::-1]
                aux = np.array(aux)[::-1]
                m = np.min(aux)
                v = np.sum(aux)
                self.PATHS.append(path)
                self.MIN.append(m) 
                self.VOL.append(v)
                for i in range(len(path[:-1])):
                    D[path[i],path[i+1]] = (D[path[i],path[i+1]]-m)*dummy
                mask = (D>0).astype(int)
                
                if self.verbose:
                    print('Computed Tunnel Number ', len(self.PATHS),v*M,m*M)
                    
                if plot:
                    self.plot_tunnel(path,D)
                
                
            else:
                self.MIN = np.array(self.MIN)*M
                self.VOL = np.array(self.VOL)*M

                if self.verbose:
                    print('No More Tunnels Available!')
                break
                

    def plot_reeb(self,):
        """! @brief Visualise the Reeb graph as a 3D network plot.

        Retrieves node measures from the graph attributes and uses them as
        colour values for the 3D scatter plot.
        """
        try:
            fun_ = nx.get_node_attributes(self.G, 'measure')
            fun = np.array([fun_[v] for v in self.G.nodes])
        except:
            fun='blue'
        
        self.network_plot_3D(fun)

        
    def plot_tunnel(self,path, D):
        """! @brief Visualise a specific tunnel path on the Reeb graph.

        Highlights the nodes belonging to the tunnel path in red and all other
        nodes in blue, then renders the 3D network plot.

        @param path Array of node indices forming the tunnel path.
        @param D Sparse or dense weight matrix used to colour depleted edges.
        """
        fun = np.array(['blue' for v in self.G.nodes])
        
        for v in path:
            fun[v] = 'red'
            
        self.network_plot_3D(fun, D)

                
    def network_plot_3D(self, fun, D = None):
        """! @brief Create a 3D network visualisation of the Reeb graph using matplotlib.

        Plots nodes as a 3D scatter plot coloured by fun, and draws edges as lines.
        Edges with zero weight in D (if provided) are drawn in red.

        @param fun Node colour values: a numpy array of scalars or colour strings.
        @param D Optional weight matrix; edges with D[i,j]==0 are coloured red (default None).
        """
        # Get node positions
        pos = nx.get_node_attributes(self.G, 'pos')
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        V = np.array([pos[key] for key in pos.keys()])    

        v_x = np.min(V[:,0])
        V_x = np.max(V[:,0])
        v_y = np.min(V[:,1])
        V_y = np.max(V[:,1])
        v_z = np.min(V[:,-1])
        V_z = np.max(V[:,-1])        
        
        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(projection='3d')
        plt.xlim([v_x, V_x])
        plt.ylim([v_y, V_y])
        ax.set_zlim([v_z, V_z])
        
        ax.scatter(V[:,0],V[:,1],V[:,2], c=fun)

        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,j in enumerate(self.G.edges()):
            
            if not D is None and D[j[0],j[1]]==0:
                c = 'r'
            else:
                c='k'

            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))

        # Plot the connecting lines
            ax.plot(x, y, z, c=c, alpha=0.5)  

        plt.tight_layout()
        plt.show() 
        
        
    def plot_res_grid(self,three_d = False, heights_=None, axis=-1):
        """! @brief Visualise the residual grid points.

        In 3D mode, renders all residual grid points as a scatter plot. In 2D mode,
        plots cross-sectional slices of the residual grid at specified heights along
        the given axis.

        @param three_d If True, produce a single 3D scatter plot; otherwise plot 2D slices (default False).
        @param heights_ Optional array of height values at which to plot 2D slices (default None, uses all unique values).
        @param axis The coordinate axis perpendicular to the slicing planes (default -1, i.e. z-axis).
        """
        if three_d:
            C = self.res_grid

            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(projection='3d')

            v_x = np.min(C[:,0])
            V_x = np.max(C[:,0])
            v_y = np.min(C[:,1])
            V_y = np.max(C[:,1])
            v_z = np.min(C[:,-1])
            V_z = np.max(C[:,-1])


            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(projection='3d')
            plt.xlim([v_x, V_x])
            plt.ylim([v_y, V_y])
            ax.set_zlim([v_z, V_z])

            ax.scatter(C[:,0],C[:,1],C[:,2], c='b', alpha=0.3)

            plt.tight_layout()
            plt.show() 
        else:
            
            if heights_ is None:
                heights = np.unique(self.res_grid[:,axis])
            else:
                heights = heights_
                
            for h in heights:

       #         print(h)
                
                fig = plt.figure()
                ax1 = fig.add_subplot(111)

                ax1.set_xlim([self.m_x, self.M_x])
                ax1.set_ylim([self.m_y, self.M_y])

                tmp = np.where(self.res_grid[:,axis] == h)[0]

                ax1.scatter(self.res_grid[tmp][:,(axis+1)%3],self.res_grid[tmp][:,(axis+2)%3], c='b')


    def plot_path_connected_comp(self, comps, three_d = False, boundaries=False, plot_reeb=False, plot_sorrounding=True):
        """! @brief Plot path-connected components of the Reeb graph in the residual grid.

        Visualises the residual grid points belonging to the specified Reeb graph
        components, optionally overlaying the Reeb graph edges and surrounding
        backbone atoms.

        @param comps Array or list of Reeb graph vertex indices identifying the components to plot.
        @param three_d If True, produce a 3D scatter plot; otherwise plot 2D slices (default False).
        @param boundaries If True, include boundary nodes (first and last) in the selection (default False).
        @param plot_reeb If True, overlay the Reeb graph edges and nodes on the plot (default False).
        @param plot_sorrounding If True, overlay surrounding backbone atoms on the plot (default True).
        """
        
        comp = []
        
        for v in comps:
            if (v>1 and v<self.n) or boundaries:
                tmp = list(np.argwhere(self.REEB_GRAPH==v)[:,0])
                comp += tmp
        
        
        comp = np.unique(comp)        
        C = self.res_grid[comp]
        
        if not three_d:
            axis=-1
            for h in np.unique(self.res_grid[:,axis]):

                fig = plt.figure()
                ax1 = fig.add_subplot(111)

                ax1.set_xlim([self.m_x, self.M_x])
                ax1.set_ylim([self.m_y, self.M_y])

                tmp = np.where(self.res_grid[:,axis] == h)[0]
                ax1.scatter(self.res_grid[tmp][:,(axis+1)%3],self.res_grid[tmp][:,(axis+2)%3], c='b')

                tmp = np.where(C[:,axis] == h)[0]
                ax1.scatter(C[tmp][:,(axis+1)%3],C[tmp][:,(axis+2)%3], c='g')

        else:
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot(projection='3d')
            plt.xlim([self.m_x, self.M_x])
            plt.ylim([self.m_y, self.M_y])
            ax.set_zlim([self.m_z, self.M_z])

            #ax.scatter(C[:,0],C[:,1],C[:,2], c='b', alpha=0.5)
            ax.scatter(C[:,0],C[:,1],C[:,2], c=C[:,-1], alpha=0.5)

            if plot_reeb: 
                try:
                    fun_ = nx.get_node_attributes(self.G, 'measure')
                    fun = np.array([fun_[v] for v in self.G.nodes])
                except:
                    fun='blue'

                pos = nx.get_node_attributes(self.G, 'pos')
                V = np.array([pos[key] for key in pos.keys()])    
                ax.scatter(V[:,0],V[:,1],V[:,2], c=fun)

                # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
                # Those two points are the extrema of the line to be plotted
                for i,j in enumerate(self.G.edges()):

                    x = np.array((pos[j[0]][0], pos[j[1]][0]))
                    y = np.array((pos[j[0]][1], pos[j[1]][1]))
                    z = np.array((pos[j[0]][2], pos[j[1]][2]))

                # Plot the connecting lines
                    ax.plot(x, y, z, c='black', alpha=0.3) 

            if plot_sorrounding:
                _,P,S = self.compute_sorrounding_structure(comps,boundaries=boundaries) 

                if len(P)>0:
                    ax.scatter(P[:,0],P[:,1],P[:,2], c='g')
                if len(S)>0:
                    ax.scatter(S[:,0],S[:,1],S[:,2], c='r')

            plt.tight_layout()
            plt.show()

        
    def make_tunnel_matrix(self,):
        """! @brief Create a tunnel capacity matrix scaled by area and particle size.

        Converts the Reeb graph weight matrix into a tunnel capacity matrix by
        scaling edge weights by the ratio of the grid cell area to the particle
        cross-sectional area, then rounding up.

        @return Numpy array with the scaled tunnel capacity matrix.
        """
        D = np.copy(self.D_reeb_w)
        area_unit = self.d_x*self.d_y
        D = D*(area_unit/(np.pi*np.power(self.d_Li/2,2)))
        D = np.ceil(D)
        return D
    
    
    def compute_bottlenecks(self, D_ = None, plot=False, col_thresh = -0.001, plot_thresh = 400, include_deadlocks = False):
        """! @brief Identify bottleneck nodes in the Reeb graph flow network.

        Computes a divergence-like quantity at each node by examining the difference
        between incoming and outgoing weights. Nodes with significant net flow imbalance
        indicate bottlenecks. Results are stored in self.DIV.

        @param D_ Optional weight matrix to use instead of self.D_reeb_w (default None).
        @param plot If True, plot the divergence profile and a 3D network coloured by bottleneck status (default False).
        @param col_thresh Threshold for colouring nodes as bottlenecks in the plot (default -0.001).
        @param plot_thresh Threshold for filtering divergence values in the line plot (default 400).
        @param include_deadlocks If True, include dead-end nodes in the divergence computation (default False).
        """
        if D_ is None:
            D_w = np.copy(self.D_reeb_w)
        else:
            D_w = np.copy(D_)
                    
        D_w = D_w - D_w.T
        check = np.array([(np.min(D_w[i,:])<0)*(np.max(D_w[i,:])>0) for i in range(self.n+1)])
    
        self.DIV = np.zeros((self.n+1,))
        
        if np.max(np.abs(D_w))>0:
            
            if include_deadlocks:
                check = np.ones_like(D_w[0,:])
            elif self.simply_connected_top_bottom:
                check[0]=1
                check[-1]=1
            
            for i in range(self.n+1):
                self.DIV[i]=np.dot(D_w[i,:],check)

        if plot:
            
            if self.simply_connected_top_bottom:
                aux = np.abs(self.DIV)<plot_thresh
            else:
                aux = np.abs(self.DIV)>-1
                
            print(aux)
            
            plt.plot(np.arange(self.n+1)[aux],self.DIV[aux])
            plt.show()
            
            self.network_plot_3D(self.DIV>=col_thresh, D = self.D_reeb_w)
            
            
    def compute_sorrounding_structure(self,path,boundaries=False):
        """! @brief Find backbone atoms surrounding a given Reeb graph path.

        For each vertex in the path, collects the corresponding residual grid points,
        then finds the nearest backbone atoms to those grid points using periodic
        distance calculations.

        @param path Array or list of Reeb graph vertex indices (not grid indices).
        @param boundaries If True, include boundary vertices in the selection (default False).
        @return Tuple (cloud, cloud_P, cloud_S) with all nearby atoms, primary atoms,
                and secondary atoms respectively.
        """
        cloud_P = []
        cloud_S = []
        cloud = []

        comp = []
        
        for v in path:
            if (v>1 and v<self.n) or boundaries:
#            n_comp = self.REEB_GRAPH[v,:]
            
#            for j,n in enumerate(n_comp):
#                if n>0:
                tmp = list(np.argwhere(self.REEB_GRAPH==v)[:,0])
                comp += tmp
        
        
        comp = np.unique(comp)
        
        if self.z_axis_boundaries:
            axes_aux = [0,1,2]
        else:
            axes_aux = np.copy(self.axes)
            
        D = dist_from_pts_periodic_boundaries(self.res_grid[comp],self.balls_centres,self.M,self.m,axes_aux)

        idxs = np.argmin(D,axis=-1)
        idxs = np.unique(idxs)
        cloud = self.balls_centres[idxs]
        
        if self.inputfile:             
            idxs_P = idxs[idxs<self.P.shape[1]]
            cloud_P = self.balls_centres[idxs_P]

            idxs_S = idxs[idxs>=self.P.shape[1]]
            cloud_S = self.balls_centres[idxs_S]
        else:
            cloud_P = cloud
            cloud_S = []
            
        return cloud, cloud_P, cloud_S

    
    
    
    