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
    
    """Reeb Graph Class"""
    
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
        """
        LEGEND
        radii  -> vector or real value determining the radius around the atoms of the backbone, inducing the void region;
        axes -> the axis along which keep the boundary conditions. The remaining axis is "unrolled";
        PREV_DATA -> use to pass from previous calculations the distances between grid and backbone;
        simply_connected_top_bottom -> consider top and bottom layer as simply connected, 
                                       ignoring the forbidden region in such sets;
        swap_res_grid_and_balls -> swap forbidden region and residual grid;
        fat_radius -> scalar saying how many times the diamater of a cell is multiplied to obtain a neighborhood of the graph;
        reeb_stride -> the stride along the axis to make the Reeb graph;
        covering -> e.g. [-1,0] or [-1,1] determine the thickness of the "level sets": (v-delta,v) or (v-delta,v+delta). 
                    Other values different from -1,1 can be tried but can cause topological errors in the graph;
        relax_z_axis -> [a,-b] puts the first a layers of the cube on top of it and the last b on the bottom: less or equal than a, 
                 greater then -b;
        transform_points -> change atoms location with affine transformation Vxpts+v
        stride -> the stride to subdivide the calculations of the distances between backbone and the grid
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
        if not self.transform_points is None:
            balls_centres, balls_radii = self.transform_backbone(balls_centres = balls_centres ,balls_radii = balls_radii,
                                M = M, m = m)

        self.produce_res_grid(balls_centres = balls_centres ,balls_radii = balls_radii,
                                        M = M, m = m)  
                 

    def transform_backbone(self, balls_centres = None, balls_radii = None,
                                M = None, m = None): 
        
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
                        d = dist_from_pts_periodic_boundaries_numba(balls_centres[i:i+self.stride,:],
                                                                    grid_aux,self.M,self.m,axes_aux)
 
                        self.distances_to_balls = np.vstack([self.distances_to_balls,d])
                        r_aux = np.argmin(self.distances_to_balls,axis=0).astype(int)
                        self.distances_to_balls = np.min(self.distances_to_balls,axis=0)
                        # print("distances to balls: ", self.distances_to_balls, " with shape ", self.distances_to_balls.shape)
                        # print("r_aux: ", r_aux, " with shape ", r_aux.shape)
                        # print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)    
                        self.atoms_kinds[r_aux>0] = r_aux[r_aux>0]-1+i 
            else:
                D = dist_from_pts_periodic_boundaries_numba(balls_centres,grid_aux,self.M,self.m,axes_aux)
                self.atoms_kinds = np.argmin(D,axis=0)
                self.distances_to_balls = np.min(D,axis=0)
                
        # print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)
        # print("balls radii: ", balls_radii, " with shape ", balls_radii.shape)
        # print("distances to balls: ", self.distances_to_balls, " with shape ", self.distances_to_balls.shape)
        print("atoms kinds: ", self.atoms_kinds, " with shape ", self.atoms_kinds.shape)
        print("distances to balls: ", self.distances_to_balls)
        thresh_radii = balls_radii[self.atoms_kinds]
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
        
        self.axis = axis
        g = lambda x: len(x)
        
        f = grid[:,axis]
        values = np.unique(f)
        idxs = np.arange(0,len(f))
        l,u = self.covering
        
        if D_ is None:
            D_graph = dist_from_pts_periodic_boundaries_numba(grid,grid,self.M,self.m,self.axes, self.dim)
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
        
        self.axis = axis
        g = lambda x: len(x)
        
        f = grid[:,axis]
        values = np.unique(f)
        idxs = np.arange(0,len(f))
        l,u = self.covering
        
        if D_ is None:
            D_graph = dist_from_pts_periodic_boundaries_numba(grid,grid,self.M,self.m,self.axes, self.dim)
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
        
        try:
            fun_ = nx.get_node_attributes(self.G, 'measure')
            fun = np.array([fun_[v] for v in self.G.nodes])
        except:
            fun='blue'
        
        self.network_plot_3D(fun)

        
    def plot_tunnel(self,path, D):
       
        fun = np.array(['blue' for v in self.G.nodes])
        
        for v in path:
            fun[v] = 'red'
            
        self.network_plot_3D(fun, D)

                
    def network_plot_3D(self, fun, D = None):
        
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
        """
        comps -> a set of path connected components, indexed as the vertices of the reeb graph
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
        D = np.copy(self.D_reeb_w)
        area_unit = self.d_x*self.d_y
        D = D*(area_unit/(np.pi*np.power(self.d_Li/2,2)))
        D = np.ceil(D)
        return D
    
    
    def compute_bottlenecks(self, D_ = None, plot=False, col_thresh = -0.001, plot_thresh = 400, include_deadlocks = False):

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
        """
        path -> a subset of the vertices of the reeb graph (not the indexes! the actual vertices!)
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

    
    
    
    