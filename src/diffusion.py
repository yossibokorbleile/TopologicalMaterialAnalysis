import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.spatial.distance import pdist, squareform, cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components
from scipy.sparse.linalg import lsqr

from itertools import chain, combinations, product

from diffusion_utils import *
import galois

import networkx as nx

import multiprocessing as mp


class Diffusion(object):
    
    """Diffusion Class"""
    
    def __init__(self, inputfile = None, box = True, r_Li = 0, d_P = 2.5, d_S = 2, sign_var = -1,
                                grid_size = 22, axes = [0,1], dim = 3,  
                                simply_connected_top_bottom = False,
                                swap_res_grid_and_balls = False,
                                close_loops_with_boundary_conditions = False,
                                z_axis_boundaries = True,
                                simplify_grid_with_triangles = True, 
                                cut_paths_until_A_0 = True,
                                a = 1, b = 0.5, c=0.1,
                                K_graph = 1,
                                MP = True,
                                verbose = False):
        """
        LEGEND
        
        simply_connected_top_bottom -> consider A_0 and A_1 as simply connected, ignoring the forbidden region in such sets.
        swap_res_grid_and_balls -> swap forbidden region and residual grid.
        close_loops_with_boundary_conditions -> consider all possible combinations for closing loops, 
                                                as made available by boundary conditions. 
                                                Requires A_0 and A_1 simply connected and geodesically convex.
        z_axis_boundaries -> consider periodicity also on the z-axis to setup the forbidden region.
        simplify_grid_with_triangles -> make sure that every edge is in a triangle. 
        cut_paths_until_A_0 -> when simplifying a path, cuts the part of the path which is in A_0, but for the latest point.
        homological_simplification -> when simplifying a path, uses homology to shorten it.
        K_graph -> coefficient for building a graph on the grid when check path connected comp. Keep it equal to 1. 
        MP -> using multiprocessing when checking if to paths are homologous.
        """   
    
    
        """
        Parameters Setup
        """
        self.inputfile = inputfile
        
        self.r_Li = r_Li
        self.d_P = d_P
        self.d_S = d_S
        self.sign_var = sign_var
        
        self.grid_size = grid_size
        self.axes = axes
        self.dim = dim
        self.MP = MP
        self.swap_res_grid_and_balls = swap_res_grid_and_balls
        self.close_loops_with_boundary_conditions = close_loops_with_boundary_conditions
        
        if self.close_loops_with_boundary_conditions:
            self.simply_connected_top_bottom = True
        else:
            self.simply_connected_top_bottom = simply_connected_top_bottom
            
        
        self.z_axis_boundaries = z_axis_boundaries
        self.simplify_grid_with_triangles = simplify_grid_with_triangles
        self.cut_paths_until_A_0 = cut_paths_until_A_0
        
        self.a = a
        self.b = b
        self.c = c
        self.K_graph = K_graph
        
        self.verbose = verbose

        if not self.inputfile is None:
            self.Li,self.P,self.S,self.N,self.M,self.m = out_to_atoms(inputfile, box=box)
            
            if not box:
                
                self.Li = self.m + np.remainder(self.Li-self.m, self.M-self.m)    
                self.P = self.m + np.remainder(self.P-self.m, self.M-self.m)    
                self.S = self.m + np.remainder(self.S-self.m, self.M-self.m)                
                
            
            """
            Functions to Setup the Grid and Triangulation
            """
            
            self.setup_triangulation()

            if self.verbose:
                print('Triangulation Done.')
            
        
    """
    Auxiliary Functions
    """    
    def produce_fobidden_region(self, balls_centres_ = None ,balls_radii_ = None,
                                M_ = None, m_ = None): 
        
        if not self.inputfile is None: 
            
            if len(self.M)>0:
                self.M_x,self.M_y,self.M_z = self.M
                self.m_x,self.m_y,self.m_z = self.m

            n_Li = self.Li.shape[0]
            self.Li = clean_Li(self.Li,self.M,self.m)
            new_n_Li = self.Li.shape[0]

            if self.verbose:
                print('Before we had ', n_Li, ' Lithium Atoms.')
                print('Now we have ', new_n_Li, ' Lithium Atoms.')


            """
            Setup Forbidden Region
            """
            P_centres = np.mean(self.P,axis=0)
            aux_P = np.argwhere(np.max(np.abs(np.linalg.norm(self.P-P_centres,axis=-1)), axis=0)<(self.M_x-self.m_x)/2)       
            self.r_P = self.d_P + self.sign_var*np.var(np.linalg.norm(self.P-P_centres,axis=-1)[:,aux_P])

            S_centres = np.mean(self.S,axis=0)
            aux_S = np.argwhere(np.max(np.abs(np.linalg.norm(self.S-S_centres,axis=-1)), axis=0)<(self.M_x-self.m_x)/2)       
            self.r_S = self.d_S + self.sign_var*np.var(np.linalg.norm(self.S-S_centres,axis=-1)[:,aux_S])

            balls_centres = np.concatenate([P_centres, S_centres])
    #        balls_radii = np.concatenate([self.r_P, self.r_S]) + self.r_Li
            balls_radii = np.array([self.r_P for _ in P_centres] + [self.r_S for _ in S_centres])  + self.r_Li
        else:
            self.M = M_
            self.m = m_
            
            if len(self.M)>0:
                self.M_x,self.M_y,self.M_z = self.M
                self.m_x,self.m_y,self.m_z = self.m
                
            balls_centres = balls_centres_
            balls_radii = balls_radii_
               
        return balls_centres, balls_radii
    
    
    def setup_triangulation(self, balls_centres_ = None ,balls_radii_ = None,
                                M_ = None, m_ = None):
    
        if (self.d_P is None or self.d_S is None) and not self.inputfile is None:
            self.estimate_atoms_min_distances(n=10)

        self.produce_res_grid(balls_centres_ = balls_centres_ ,balls_radii_ = balls_radii_,
                            M_ = M_, m_ = m_)
        self.setup_geodesics_and_boundaries()
        self.make_edges()
        self.make_triangulation()
    
    def estimate_atoms_min_distances(self,n=10):
        
        tmp = np.arange(0,len(self.Li),n)

        if self.z_axis_boundaries:
            axes_aux = [0,1,2]
        else:
            axes_aux = np.copy(self.axes)

        d_Li = []
        d_P = []
        d_S = []
        
        n_Li = self.Li.shape[1]
        n_P = self.P.shape[1]
        n_S = self.S.shape[1]
        
        for t in tmp:
            L = self.Li[t,:,:] 
            P = self.P[t,:,:]
            S = self.S[t,:,:] 
            
            B = np.concatenate([L,P,S])
            
            D = dist_from_pts_periodic_boundaries(B,L,self.M,self.m,axes_aux)
            D[:n_Li,:] = D[:n_Li,:] + np.identity(n_Li)*100
            
            d_Li = np.concatenate([d_Li,np.min(D[:n_Li,:],axis=-1)])
            d_P = np.concatenate([d_P,np.min(D[n_Li:n_Li+n_P,:],axis=-1)])
            d_S = np.concatenate([d_S,np.min(D[n_Li+n_P:n_Li+n_P+n_S,:],axis=-1)])
            
        
        fig = plt.figure()
        
        bp = plt.boxplot(d_Li,showmeans=True)
        self.d_Li = bp['caps'][0].get_ydata()[0] 
        
        bp = plt.boxplot(d_P,showmeans=True)
        self.d_P = bp['caps'][0].get_ydata()[0] 
        
        bp = plt.boxplot(d_S,showmeans=True)
        self.d_S = bp['caps'][0].get_ydata()[0]    

        
        plt.close(fig)
        
        if self.verbose:
            print('Estimated Radii (Li, P, S): ',self.d_Li,self.d_P, self.d_S)
            
    def produce_res_grid(self, balls_centres_ = None ,balls_radii_ = None,
                                M_ = None, m_ = None):
        
        balls_centres, balls_radii = self.produce_fobidden_region(balls_centres_ = balls_centres_,
                                                                  balls_radii_ = balls_radii_,
                                                                  M_ = M_, m_ = m_)

        if self.verbose and not self.inputfile is None:
            print('The ranges of the radii are: ', self.r_P, self.r_S)
        
        tmp = 1

        self.nx, self.ny, self.nz = (self.grid_size, self.grid_size, self.grid_size)
        
        self.x = np.linspace(self.m_x, self.M_x, self.nx+tmp)[:-1]
        self.y = np.linspace(self.m_y, self.M_y, self.ny+tmp)[:-1]
        self.z = np.linspace(self.m_z, self.M_z, self.nz+tmp)[:-1]
            
        self.d_x = self.x[1]-self.x[0]
        self.d_y = self.y[1]-self.y[0]
        self.d_z = self.z[1]-self.z[0]

        k = 1-0.999

        self.z_0 = self.z[0]+self.d_z*k
        self.z_1 = self.z[-1]-self.d_z*k

        xv, yv, zv = np.meshgrid(self.x, self.y, self.z)

        self.grid = np.zeros((self.nx*self.ny*self.nz,self.dim))
        self.N = len(self.grid)
        
        self.grid[:,0] = xv.flatten()
        self.grid[:,1] = yv.flatten()
        self.grid[:,2] = zv.flatten()

        labels = np.ones_like(self.grid[:,0], dtype=int)
        aux = []


        if self.z_axis_boundaries:
            axes_aux = [0,1,2]
        else:
            axes_aux = np.copy(self.axes)

        D = dist_from_pts_periodic_boundaries(balls_centres,self.grid,self.M,self.m,axes_aux)

        for i,p in enumerate(balls_centres):
            tmp = list(np.argwhere(D[i,:]<balls_radii[i]))
            if self.simply_connected_top_bottom:
                tmp = [q for q in tmp if self.grid[q[0]][-1]<self.z_1]
                tmp = [q for q in tmp if self.grid[q[0]][-1]>self.z_0]
            aux=aux+[idx[0] for idx in tmp]

        idxs = set(np.arange(self.grid.shape[0]))-set(aux)
        idxs_ = np.sort(list(idxs))

        N=len(idxs_)

        if self.verbose:
            print('The Initial Residual Grid has: ', N,' points.')
        
        if self.swap_res_grid_and_balls:
            if self.simply_connected_top_bottom:
                idxs = np.array([i for i in np.arange(len(self.grid)) if i not in idxs_ 
                                                                      or self.grid[i][-1]<self.z_0 
                                                                      or self.grid[i][-1]>self.z_1])
                idxs_ = np.array([i for i in np.arange(len(self.grid)) if i not in idxs])
            else:
                idxs = np.array([i for i in np.arange(len(self.grid)) if i not in idxs_])
                idxs_ = np.array([i for i in np.arange(len(self.grid)) if i not in idxs])
                
            r = self.d_x*0
            self.balls_centres = self.grid[idxs_]  
            self.balls_radii = np.ones_like(idxs_)*r
        else:
            idxs = np.copy(idxs_)
            self.balls_centres = balls_centres
            self.balls_radii = balls_radii
            
        N=len(idxs)
        res_grid=self.grid[idxs]
        
        """
        Get Rid of Unwanted Connected Components
        """
        D = dist_periodic_boundaries(res_grid,self.M,self.m,self.axes, dim=self.dim)
        self.r_graph = np.sqrt(np.max([self.d_x,self.d_y,self.d_z])**2)*np.sqrt(2)*1.001*self.K_graph
        
        if self.verbose:
            print('Recall: "r_graph" is hard-coded.')
        
        graph = D<self.r_graph

        n_c,C = connected_components(graph)
        n_comp = np.argmax([np.sum(C==i) for i in range(n_c)])
        comp = np.argwhere(C==n_comp)

        comp = np.array([p[0] for p in comp])

        if self.verbose:
            print('We have ',n_c,' connected components, with the following cardinalities: ',
                  [np.sum(C==i) for i in range(n_c)], comp.shape)

        self.res_grid = np.array([p for i,p in enumerate(res_grid) if i in comp])
        self.N=len(self.res_grid)

        if self.verbose:
            print('The Intermediate Residual Grid has: ', self.res_grid.shape[0],' points.')
            
        """
        Make sure that Every Vertex and Every Edges are in a Triangle
        """
        if self.simplify_grid_with_triangles:
            """
            First Round: Remove Vertices which are not in Triangles
            """
            self.D_force_graph = np.ones((self.N,self.N))
            self.make_edges()            
            self.make_triangulation()

            vertices = []
            
            for tr in self.ITRIS:
                vertices += list(tr)

            vertices = np.unique(vertices)  
            self.res_grid = np.array([p for i,p in enumerate(self.res_grid) if i in vertices])                               
            self.N=len(self.res_grid)
            
            """
            Second Round: Remove Edges which are not in Triangles
            """            
            self.D_force_graph = np.ones((self.N,self.N))
            self.make_edges()            
            self.make_triangulation()
           
            self.D_force_graph = np.zeros((self.N,self.N))

            for tr in self.ITRIS:
                self.D_force_graph[tr[0],tr[1]] = 1
                self.D_force_graph[tr[1],tr[0]] = 1
                self.D_force_graph[tr[1],tr[2]] = 1
                self.D_force_graph[tr[2],tr[1]] = 1
                self.D_force_graph[tr[0],tr[2]] = 1
                self.D_force_graph[tr[2],tr[0]] = 1
            """
            Check Path Connected Components
            """            
            n_c,C = connected_components(self.D_force_graph)
            n_comp = np.argmax([np.sum(C==i) for i in range(n_c)])
            comp = np.argwhere(C==n_comp)

            comp = np.array([p[0] for p in comp])

            if self.verbose:
                print('We have ',n_c,' connected components, with the following cardinalities: ',
                      [np.sum(C==i) for i in range(n_c)], comp.shape)            
            
            D_aux = np.copy(self.D_force_graph)
            res_grid_aux = np.copy(self.res_grid)
            
            self.res_grid = np.array([p for i,p in enumerate(self.res_grid) if i in comp])
            self.N=len(self.res_grid)
            self.D_force_graph = np.zeros((self.N,self.N))
            
            idxs_old = {}
            
            for i in range(self.N):
                idxs_old[i] = np.argwhere(np.linalg.norm(res_grid_aux-self.res_grid[i],axis=-1)<0.000001)[0][0]
            
            tmp_ = np.array([C==n_comp])
            tmp = tmp_.T @ tmp_
            
            for i in range(self.N):
                i_old = idxs_old[i]
                for j in range(i):
                    j_old = idxs_old[j]
                    self.D_force_graph[i,j] = D_aux[i_old,j_old]*tmp[i_old,j_old]
                    self.D_force_graph[j,i] = D_aux[i_old,j_old]*tmp[i_old,j_old]
                    
            """
            Final Round
            """
           
            self.make_edges()            
            self.make_triangulation()

            if self.verbose:
                print('Triangulation Check on the Number of Edges: ', np.sum(self.D_force_graph), 2*len(self.idxs_to_edge))
        
        if not self.simplify_grid_with_triangles:
            self.D_force_graph = np.ones((self.N,self.N))
        
        if self.verbose:
            print('The Final Residual Grid has: ', self.N,' points.')

        """
        Setup Indexes to Go Back and Forth from grid to res_grid
        """
        idxs = []

        for p in self.res_grid:
            tmp = np.argmin(np.linalg.norm(self.grid-p,axis=1))
            idxs.append(tmp)

        self.idxs= np.array(idxs)

        self.grid_to_res_grid={}

        for i in range(len(self.grid)):
            tmp = np.argwhere(self.idxs==i)
            if len(tmp)>0:
                self.grid_to_res_grid[i]=tmp[0][0]
            else:
                self.grid_to_res_grid[i]=-1
        """
        Lastly Obtain the Complementary
        """
        tmp = np.array([i for i in np.arange(len(self.grid)) if i not in self.idxs])
        
        if len(tmp)>0:
            self.balls = self.grid[np.array([i for i in np.arange(len(self.grid)) if i not in self.idxs])]  
        else:
            self.balls = tmp

        
    def setup_geodesics_and_boundaries(self,):

        self.D, self.geod = make_geodesics(self.res_grid, self.r_graph, self.M, self.m, self.axes, 
                                           self.dim, 
                                           D_force_graph_ = self.D_force_graph)

        self.idxs_A_0 = np.array([t[0] for t in np.argwhere(self.res_grid[:,-1]<self.z_0)])        
        self.A_0 = self.res_grid[self.idxs_A_0]

        self.A_0_to_res={}
        self.res_to_A_0 = {}

        for i,p in enumerate(self.A_0): 
            tmp = np.argwhere(np.linalg.norm(self.res_grid-p, axis=1)<0.000001)[0][0]
            self.A_0_to_res[i] = tmp
            self.res_to_A_0[tmp] = i

        _, self.geod_A_0 = make_geodesics(self.A_0,self.r_graph, self.M, self.m, self.axes, 
                                          self.dim, 
                                          D_force_graph_ = self.D_force_graph[self.idxs_A_0,:][:,self.idxs_A_0])    

        self.idxs_A_1 = np.array([t[0] for t in np.argwhere(self.res_grid[:,-1]>self.z_1)])
        self.A_1 = self.res_grid[self.idxs_A_1]

        self.A_1_to_res={}
        self.res_to_A_1 = {}

        for i,p in enumerate(self.A_1): 
            tmp = np.argwhere(np.linalg.norm(self.res_grid-p, axis=1)<0.000001)[0][0]
            self.A_1_to_res[i] = tmp
            self.res_to_A_1[tmp] = i

        _, self.geod_A_1 = make_geodesics(self.A_1,self.r_graph, self.M, self.m, self.axes, 
                                          self.dim, 
                                          D_force_graph_ = self.D_force_graph[self.idxs_A_1,:][:,self.idxs_A_1]) 


    def make_edges(self,):

        self.edge_idxs = {}
        idxs_to_edge = []
        self.r_graph_aux = np.sqrt(self.d_x**2 + self.d_y**2+self.d_z**2)*1.001*self.K_graph

        if self.z_axis_boundaries:
            axes_aux = [0,1,2]
        else:
            axes_aux = np.copy(self.axes)

        D_e = dist_periodic_boundaries(self.res_grid,self.M,self.m,self.axes, self.dim)

        aux = 0
        for i in range(D_e.shape[0]):
            for j in range(i):
                if D_e[i,j]<self.r_graph_aux and self.D_force_graph[i,j]>0:
                    self.edge_idxs[(j,i)] = aux
                    idxs_to_edge.append([j,i])
                    aux+=1
                    
        self.idxs_to_edge = np.array(idxs_to_edge)
        self.N_1 = aux
        
    def make_cubic_to_res(self, cubic_grid):
    
        self.cubic_to_res = {}

        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    p = cubic_grid[i,j,k,:]
                    tmp = np.argwhere(np.linalg.norm(self.res_grid-p,axis=1)<0.0000001)

                    if len(tmp)>0:
                        self.cubic_to_res[(i,j,k)] = tmp[0][0]
                    else:
                        self.cubic_to_res[(i,j,k)] = -1

    def make_triangulation(self,):

        cubic_grid = np.transpose(self.grid.reshape((self.nx,self.ny,self.nz,3)),axes=(1,0,2,3))
        self.make_cubic_to_res(cubic_grid)

        self.ITRIS = []
        N_2 = (self.N-self.nx-self.ny+1)*4*4

        B = np.zeros((N_2,self.N_1), dtype=int)

        a_x = int(0 not in self.axes) + self.nx
        a_y = int(1 not in self.axes) + self.ny
        a_z = int(2 not in self.axes) + self.nz

        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):

                    p = cubic_grid[i,j,k,:]

                    for tr in [[(0,0,0),(1,0,0),(1,1,0)],
                                 [(0,0,0),(0,1,0),(1,1,0)],
                                 [(0,0,0),(1,0,0),(0,1,0)],
                                 [(1,0,0),(0,1,0),(1,1,0)],
                                 [(0,0,0),(1,0,0),(1,0,1)],
                                 [(0,0,0),(0,0,1),(1,0,1)],
                                 [(0,0,0),(1,0,0),(0,0,1)],
                                 [(1,0,0),(0,0,1),(1,0,1)], 
                                 [(0,0,0),(0,1,0),(0,1,1)],
                                 [(0,0,0),(0,0,1),(0,1,1)],
                                 [(0,0,0),(0,1,0),(0,0,1)],
                                 [(0,1,0),(0,0,1),(0,1,1)]]:

                        p0 = tuple([(i+tr[0][0])%a_x,(j+tr[0][1])%a_y,(k+tr[0][2])%a_z])
                        p1 = tuple([(i+tr[1][0])%a_x,(j+tr[1][1])%a_y,(k+tr[1][2])%a_z])
                        p2 = tuple([(i+tr[2][0])%a_x,(j+tr[2][1])%a_y,(k+tr[2][2])%a_z])

                        pts = np.vstack([p0,p1,p2]) 

                        if np.max(pts[:,0])<self.nx and np.max(pts[:,1])<self.ny and np.max(pts[:,2])<self.nz: 
                            ITRIS, B = check_tr(p0,p1,p2, self.res_grid, self.grid, self.cubic_to_res, 
                                                self.ITRIS, B, self.edge_idxs)                             

        """
        Setup the 2-Boundary Matrix of the Triangulation and the 1-Boundary Matrix of A_0 and A_1 
        """
        self.ITRIS = np.array(ITRIS)  
        self.N_2=len(ITRIS)
        self.B = B[:self.N_2,:]
        
        axis = -1
        idxs_e = []
        
        for i,e in enumerate(self.idxs_to_edge):
            edge_z = self.res_grid[e,-1]
            tmp0 = np.prod([edge_z<self.z_0])
            tmp1 = np.prod([edge_z>self.z_1])
            
            if np.max([tmp0, tmp1])>0:
                idxs_e.append(i)        
 
        idxs_e = np.array(idxs_e)
        
        A = np.zeros((self.N_1,len(idxs_e)))
        
        for i,idx in enumerate(idxs_e):
            A[idx,i]+=1
            
        self.A = A
        
        if self.verbose:
            print('The Triangulation Has ',self.N_1, ' Edges and ', self.N_2, ' Triangles.')

            
    def make_N(self,):
        N = np.exp(-self.a*self.D/self.eps - (self.b*self.direction+self.c*self.repulsion)/self.eps**2)
        return N

        
    def sample_diffusion(self,):

        idx = np.random.multinomial(1,np.ones_like(self.idxs_A_0)/len(self.idxs_A_0))
        idx = np.argwhere(idx>0)[0][0]
        path = sample_path(self.idxs_A_0[idx], self.N, self.res_grid, self.z_1,)
        tmp = np.argwhere(np.linalg.norm(self.A_1-self.res_grid[path[-1]], axis=1)<0.000001)        
        return path,tmp

    def preprocess_path(self, path):
        return simplify_path(smooth_path(path, self.geod, self.res_grid, self.D))
    
    
    def cut_until_A_0(self,path):
        
        path_ = np.copy(path)
        
        while np.min(np.linalg.norm(self.A_0-self.res_grid[path_[1]],axis=-1))<0.00000001:
            path_ = path_[1:]
            
        while np.min(np.linalg.norm(self.A_1-self.res_grid[path_[-2]],axis=-1))<0.00000001:
            path_ = path_[:-1]
            
        return np.array(path_)
    
    def simplify_path_with_homology(self,path):
        
        if self.cut_paths_until_A_0:
            path_ = self.cut_until_A_0(path)
        else:
            path_ = np.copy(path)

        path_ = homological_simplification_of_path(path_, self.res_grid, self.geod, 
                                                          self.V[self.rank:,:], self.rank, self.GF, 
                                                          self.D, self.N_1, 
                                                          self.edge_idxs, self.r_graph)
        
#        path_ = self.preprocess_path(path_)
        
        if self.cut_paths_until_A_0:
            path_ = self.cut_until_A_0(path_)
        
        return path_

    
    def equivalence(self, v_path):

        b_ = self.GF(v_path.astype(int)%2)
        b = self.V[self.rank:,:]@b_

        aux = np.max(b)
        
        return aux

    def fill_loop(self, idxs_list, BOUNDARY_MAT = None, RED_MAT = None, V_MAT = None,
                  find_chain = False, plot = False):
        """
        Finds the 2-chain which fills the loop in p,q are homologous. Se é "None" prendo A_0 e A_1.
        BOUNDARY_MAT -> la matrice A, che genera il sottospazio di bord
        RED_MAT -> la matrice A_red cioé la versione row-reduced di A[rank(B):,:]
        V_MAT -> la matrice che porta A[rank(B):,:] in forma A_red
        """
              
        if len(idxs_list[0])==1: 
            paths = np.array([self.PATHS[i] for i in idxs_list],dtype=object)
        else:
            paths = idxs_list
                        
        tmp = self.prova(paths, BOUNDARY_MAT = BOUNDARY_MAT, RED_MAT = RED_MAT, V_MAT = V_MAT, verbose = False)
        
        if self.verbose:
            if tmp == 0:
                print('The paths are homologous.')
            else:
                print('The paths are NOT homologous.')
        
        if find_chain or plot:  
            
            loop = np.concatenate(paths,axis=-1,dtype=int)
            
            v = path_to_vec(loop, self.D, self.N_1, self.edge_idxs, self.r_graph, verbose = False)
            v = self.GF(v.astype(int)%2)
            v_ = self.V@v
            
            if BOUNDARY_MAT is None:
                A, A_red, V = self.A, self.A_red, self.V_
                A = self.V@A
            elif np.max(BOUNDARY_MAT)==0:
                A, A_red, V = 0,0,0
            else:
                if RED_MAT is None:
                    A = self.V@BOUNDARY_MAT
                    A_red, V = row_reduce_matrix(A[self.rank:,:],self.GF)
                else:
                    A = BOUNDARY_MAT
#                    A_red = RED_MAT
#                    A_red = V@A[self.rank:,:]                
                    A_red = RED_MAT
                    V = V_MAT              
            
            chain, tmp_ = solve_linear_system_with_boundary(self.B_red, v_, A, A_red, V, self.GF)
                       
            if plot:

                boundary = self.B@chain

                if self.verbose:
                    print('We have a loop with ', sum(v>0), ' edges and a 2-chain with ', sum(chain>0), ' 2-simplexes.')

                fig = plt.figure(figsize=(15,10))

                ax1 = fig.add_subplot(121,projection='3d')
                ax2 = fig.add_subplot(122,projection='3d')

                ax1.set_xlim([self.m_x, self.M_x])
                ax1.set_ylim([self.m_y, self.M_y])
                ax1.set_zlim([self.m_z, self.M_z])

                ax2.set_xlim([self.m_x, self.M_x])
                ax2.set_ylim([self.m_y, self.M_y])
                ax2.set_zlim([self.m_z, self.M_z])


                for tr in self.ITRIS[np.where(chain>0)[0]]:

    #                    vtx = np.vstack([self.res_grid[tr,:],[self.res_grid[tr[0],:]] ])                      
                    vtx = self.res_grid[tr,:]                    
                    tri = Poly3DCollection([vtx], alpha=0.3)
                    tri.set_facecolor('b')
                    tri.set_edgecolor('k')
                    ax1.add_collection3d(tri)

                    vtx = self.traslate_box(vtx)                   
                    tri = Poly3DCollection([vtx], alpha=0.3)
                    tri.set_facecolor('b')
                    tri.set_edgecolor('k')
                    ax2.add_collection3d(tri)

                for p in paths:
                    ax1.plot3D(self.res_grid[p,0], self.res_grid[p,1], self.res_grid[p,2], 
                                   linewidth=1,linestyle='dotted', c='red') 

                    aux = self.traslate_box(self.res_grid[p,:])
                    ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], 
                               linewidth=1,linestyle='dotted', c='red') 

                for idx in np.where(v>0)[0]:
                    e = self.idxs_to_edge[idx,:]
                    ax1.plot3D(self.res_grid[e][:,0], self.res_grid[e][:,1], self.res_grid[e][:,2], 
                               linewidth=8, c='red') 

                    aux = self.traslate_box(self.res_grid[e])
                    ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], linewidth=5, c='red') 

                for idx in np.where(boundary>0)[0]:
                    e = self.idxs_to_edge[idx,:]
                    ax1.plot3D(self.res_grid[e][:,0], self.res_grid[e][:,1], self.res_grid[e][:,2], 
                               linewidth=3, c='green') 

                    aux = self.traslate_box(self.res_grid[e])
                    ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], linewidth=3, c='green') 



                plt.show()         

                if tmp>0:
                    chain = []                
        else:
            chain = []
            
        return chain, tmp
        
    
    def cluster_deadlocks(self, PATHS):

        if len(PATHS)==0:
            return PATHS
        
        if self.verbose:
            print('Currently not working!')
        
        DEADLOCKS=[[PATHS[0]]]

        for p in PATHS[1:]:

            new = 1

            for i,c in enumerate(DEADLOCKS):
                aux = 1
                for q in c:
                    path = close_deadlocks(p,q, self.res_grid, self.geod, self.res_to_A_0,
                                           self.A_0_to_res, self.geod_A_0,self.A_0)

                    if len(path)==0:
                        aux = 0
                        break

                if aux == 1:
                    v_path = path_to_vec(path,self.D,self.N_1,self.edge_idxs, self.r_graph)
                    res = self.equivalence(v_path)
                    if res == 0:
                        new = 0
                        break

            if new == 1:    
                DEADLOCKS.append([p])
            else:
                DEADLOCKS[i].append(p)

                
                
        if self.verbose:
            print([len(c) for c in DEADLOCKS])
                
        return DEADLOCKS
    
    
    def traslate_box(self, coords, v_=None):
        
        axes_aux = np.array([i in self.axes for i in range(self.dim)])
        
        if v_ is None:
            v = np.array([(self.M_x-self.m_x)/2, (self.M_y-self.m_y)/2, (self.M_z-self.m_z)/2])
            v = v*axes_aux
            
        else:
            v = v_
            
        new_coords = self.m + np.remainder(coords-self.m + v, self.M-self.m)
        
        return new_coords
    
    
    def check_cluster(self, n_clus, n = 5, plot=False):
        
        N = len(self.CLUSTERS[n_clus])
        
        for _ in range(n):
            idx0 = np.random.multinomial(1,np.ones((N,))/N)
            idx0 = np.argwhere(idx0>0)[0][0]
            idx1 = np.random.multinomial(1,np.ones((N,))/N)
            idx1 = np.argwhere(idx1>0)[0][0]

            self.fill_loop([self.CLUSTERS[n_clus][idx0],
                           self.CLUSTERS[n_clus][idx1]], plot = plot)
    
    
    def proj_lithium_paths(self,LIST, plot=True):
    
        axes_aux = [0,1,2]
    
        for idx in LIST:
            raw_path = self.Li[:,idx,:]

            p,err = project_path_periodic(raw_path, self.res_grid,self.M,self.m,axes_aux, dim=self.dim)
            p_aux,err_full = project_path_periodic(raw_path, self.grid,self.M,self.m,axes_aux, dim=self.dim)

            if self.verbose:
                print('\nAtom: ',idx, '\nError on Res Grid : ', np.max(err),
                      '\nError on Full Grid : ', np.max(err_full), '\nGrid Diameter : ', (np.sqrt(3)/2)*self.d_x)

            if plot:

                raw = self.traslate_box(raw_path)
                grid = self.traslate_box(self.grid[p_aux])
                res = self.traslate_box(self.res_grid[p])
                balls = self.traslate_box(self.balls)
                
                fig = plt.figure(figsize=(15,15))
                
                ax1 = fig.add_subplot(221,projection='3d')
                ax2 = fig.add_subplot(222,projection='3d')
                ax3 = fig.add_subplot(223,projection='3d')
                ax4 = fig.add_subplot(224,projection='3d')
                
                ax1.plot(raw_path[:,0],raw_path[:,1],raw_path[:,2], c='b', label='Raw Path')
                ax1.plot(self.grid[p_aux,0],self.grid[p_aux,1], self.grid[p_aux,2],c='g',label='Proj on Full Grid')
                ax1.plot(self.res_grid[p,0],self.res_grid[p,1], self.res_grid[p,2],c='r', label='Proj on Res Grid')
                
                ax2.plot(raw[:,0],raw[:,1],raw[:,2], c='b', label='Raw Path')
                ax2.plot(grid[:,0],grid[:,1],grid[:,2],c='g',label='Proj on Full Grid')
                ax2.plot(res[:,0],res[:,1],res[:,2],c='r', label='Proj on Res Grid')
                
                ax3.scatter(self.balls[:,0],self.balls[:,1], self.balls[:,2],c='g',alpha=0.1)
                ax3.plot(raw_path[:,0],raw_path[:,1],raw_path[:,2], c='b', label='Raw Path')
                ax3.plot(self.grid[p_aux,0],self.grid[p_aux,1], self.grid[p_aux,2],c='g',label='Proj on Full Grid')
                ax3.plot(self.res_grid[p,0],self.res_grid[p,1], self.res_grid[p,2],c='r', label='Proj on Res Grid')

                ax4.scatter(balls[:,0],balls[:,1], balls[:,2],c='g',alpha=0.1)
                ax4.plot(raw[:,0],raw[:,1],raw[:,2], c='b', label='Raw Path')
                ax4.plot(grid[:,0],grid[:,1],grid[:,2],c='g',label='Proj on Full Grid')
                ax4.plot(res[:,0],res[:,1],res[:,2],c='r', label='Proj on Res Grid')

                plt.legend()
                plt.tight_layout()
                plt.show()
                
                
    def reeb_graph(self,grid,axis=-1, D_ = None):
        
        g = lambda x: len(x)
        
        f = grid[:,axis]
        values = np.unique(f)
        
        if self.simply_connected_top_bottom:
            values = values[1:]
        
        idxs = np.arange(0,len(f))

        eps = self.r_graph
        
        if D_ is None:
            D_graph = dist_periodic_boundaries(grid,self.M,self.m,self.axes, self.dim)
            D_graph = (D_graph<self.r_graph)
        else:
            D_graph = np.copy(D_)
            
        aux = []
        
        for i,p in enumerate(grid):
            neigh = np.where(D_graph[i,:]>0)
            aux.append(np.max(np.abs(f[neigh]-f[i])))

        delta = 1.001*np.max(aux)
        count = 0
        
        self.REEB_GRAPH = -1*np.ones((grid.shape[0],values.shape[0])).astype(int)
        E = []
        fun = []
        emb = []

        fat_lvl = (f<values[0]+delta)
        lvl_to_grid = np.where(fat_lvl>0)[0]

        n_c,C = connected_components(D_graph[:,fat_lvl][fat_lvl,:])

        for c in range(n_c):
            tmp_ = np.where(C==c)[0]
            tmp = lvl_to_grid[tmp_]
            self.REEB_GRAPH[tmp,0] = count 
            fun.append(g(grid[fat_lvl]))
            emb.append(np.mean(grid[fat_lvl],axis=0))
            count += 1


        for j,v in enumerate(values[1:]):

            i=j+1
#            sub_lvl = (f<=v+delta)
            sub_lvl = (f<=v)
            sup_lvl = (f>v-delta)
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
                    E.append((old_c,count))

                fun.append(g(grid[tmp]))
                emb.append(np.mean(grid[tmp],axis=0))

                count += 1

        fun = np.array(fun)
        emb = np.array(emb)
        self.n = len(emb)-1
        
        return E, fun, emb
    
    def build_graph(self,E,emb,fun):
    
        # create networkx graph
        G=nx.Graph()

        # add vertices
        for p in range(len(emb)):
            G.add_node(p, pos=(emb[p,0],emb[p,1],emb[p,-1]))

        # add edges
        for edge in E:
            G.add_edge(edge[0], edge[1])

        return G
    """
    Plot Functions
    """
    
    def plot_pts(self,pts):
    
        fig = plt.figure(figsize=(15,10))            
        ax1 = fig.add_subplot(121,projection='3d')
        ax2 = fig.add_subplot(122,projection='3d')

        ax1.set_xlim([self.m_x, self.M_x])
        ax1.set_ylim([self.m_y, self.M_y])
        ax1.set_zlim([self.m_z, self.M_z])

        ax2.set_xlim([self.m_x, self.M_x])
        ax2.set_ylim([self.m_y, self.M_y])
        ax2.set_zlim([self.m_z, self.M_z])


        ax1.scatter(pts[:,0],pts[:,1], 
                   pts[:,2],color='g',alpha=0.5)


        prova = self.traslate_box(pts)

        ax2.plot(prova[:,0],prova[:,1], 
                   prova[:,2],color='g',alpha=0.5)
        
        plt.show()
            
    def plot_voids(self,):
        self.plot_pts(self,self.res_grid)
            
    def plot_balls(self,):
        self.plot_pts(self,self.balls)
            
    def plot_triangulation(self, axis=-1):
        
        D_flat = squareform(pdist(self.res_grid))
        
        aux = self.traslate_box(self.res_grid)
        D_flat_aux = squareform(pdist(aux))
        

        for h in [self.x,self.y,self.z][axis]:

            fig = plt.figure(figsize=(15,10))
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            
            ax1.set_xlim([self.m_x, self.M_x])
            ax1.set_ylim([self.m_y, self.M_y])

            ax2.set_xlim([self.m_x, self.M_x])
            ax2.set_ylim([self.m_y, self.M_y])
            
            tmp = np.where(self.res_grid[:,axis] == h)[0]
            
            ax1.scatter(self.res_grid[tmp][:,(axis+1)%3],self.res_grid[tmp][:,(axis+2)%3], c='k')
            ax2.scatter(aux[tmp][:,(axis+1)%3],aux[tmp][:,(axis+2)%3], c='k')

            for edge in self.edge_idxs.keys():
                if np.prod([self.res_grid[p][axis]==h for p in edge])==1:

                    e = np.array(list(edge))
                    if np.max(D_flat[e,:][:,e])>self.r_graph_aux:
                        alpha1=0.1
                    else:
                        alpha1=1
                    ax1.plot(self.res_grid[e][:,(axis+1)%3],self.res_grid[e][:,(axis+2)%3], c='k', alpha=alpha1)

                    if np.max(D_flat_aux[e,:][:,e])>self.r_graph_aux:
                        alpha2=0.1
                    else:
                        alpha2=1
                    ax2.plot(aux[e][:,(axis+1)%3],aux[e][:,(axis+2)%3], c='k', alpha=alpha2)
                    
            for tr in self.ITRIS:
                if np.prod([self.res_grid[p][axis]==h for p in tr])==1:
                    if np.max(D_flat[tr,:][:,tr])>self.r_graph_aux:
                        alpha1=0.1
                        c1='g'
                    else:
                        alpha1=1
                        c1='b'
                    ax1.fill(self.res_grid[tr][:,(axis+1)%3],self.res_grid[tr][:,(axis+2)%3],c=c1,alpha=alpha1)    
                    
                    if np.max(D_flat_aux[tr,:][:,tr])>self.r_graph_aux:
                        alpha2=0.1
                        c2='g'
                    else:
                        alpha2=1
                        c2='b'
                        
                    ax2.fill(aux[tr][:,(axis+1)%3],aux[tr][:,(axis+2)%3],c=c2,alpha=alpha2)                      
            plt.show()

    def plot_path(self,path):
        
     #   path = self.PATHS[path_idx]
                
        fig = plt.figure(figsize=(15,10))

        ax1 = fig.add_subplot(121,projection='3d')
        ax2 = fig.add_subplot(122,projection='3d')

        ax1.set_xlim([self.m_x, self.M_x])
        ax1.set_ylim([self.m_y, self.M_y])
        ax1.set_zlim([self.m_z, self.M_z])

        ax2.set_xlim([self.m_x, self.M_x])
        ax2.set_ylim([self.m_y, self.M_y])
        ax2.set_zlim([self.m_z, self.M_z])

        ax1.plot3D(self.res_grid[path,0], self.res_grid[path,1], self.res_grid[path,2], 
                       linewidth=1,linestyle='dotted', c='red') 

        aux = self.traslate_box(self.res_grid[path,:])
        ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], 
                       linewidth=1,linestyle='dotted', c='red') 

            
    def plot_vector_field(self, axis=1):

        result = self.b*self.v + self.c*self.grad_pen
        
        for h in [self.x,self.y,self.z][axis]:
            aux = []
        
            for i,p in enumerate(self.res_grid):    
                if p[axis]==h:
                    aux.append(i) 

            aux = np.array(aux)

            plt.quiver(self.res_grid[aux][:,(axis+1)%3],self.res_grid[aux][:,(axis+2)%3],
                       result[aux][:,(axis+1)%3],result[aux][:,(axis+2)%3])
            plt.show()
            
    def plot_paths_clustered(self, min_len = 5, k = 0,aggregate = True):
    
        cmap = get_cmap(len(self.CLUSTERS))
        colors = ['r','b','g','k']
        
        aux = sum([1 for c in self.CLUSTERS if len(c)> min_len])
        cnt = -1
        
        if aggregate:
            fig = plt.figure(figsize=(15,10))            
            ax1 = fig.add_subplot(121,projection='3d')
            ax2 = fig.add_subplot(122,projection='3d')
            
            ax1.set_xlim([self.m_x, self.M_x])
            ax1.set_ylim([self.m_y, self.M_y])
            ax1.set_zlim([self.m_z, self.M_z])

            ax2.set_xlim([self.m_x, self.M_x])
            ax2.set_ylim([self.m_y, self.M_y])
            ax2.set_zlim([self.m_z, self.M_z])

        for i,c in enumerate(self.CLUSTERS):
            
            if len(c)>min_len:                
                
                cnt += 1
                
                if aux>4:
                    col = cmap(cnt)
                else:
                    col = colors[cnt]

                if not aggregate:
                    print('Cluster Cardinality: ', len(c))
                    fig = plt.figure(figsize=(15,10))
                    ax1 = fig.add_subplot(121,projection='3d')
                    ax2 = fig.add_subplot(122,projection='3d')
                
                    ax1.set_xlim([self.m_x, self.M_x])
                    ax1.set_ylim([self.m_y, self.M_y])
                    ax1.set_zlim([self.m_z, self.M_z])

                    ax2.set_xlim([self.m_x, self.M_x])
                    ax2.set_ylim([self.m_y, self.M_y])
                    ax2.set_zlim([self.m_z, self.M_z])

                for idx in c:
                    p = self.PATHS[idx]
                    noise = np.random.normal(size=(len(p),3))*k
                    #ax.plot(res_grid[p,0],res_grid[p,1], res_grid[p,2],c=col,alpha=1)
                    
                    ax1.plot(self.res_grid[p,0]+noise[:,0],self.res_grid[p,1]+noise[:,1], 
                               self.res_grid[p,2]+noise[:,2],color=col,alpha=1)
                    
                    
                    prova = self.traslate_box(self.res_grid[p,:])
                    
                    ax2.plot(prova[:,0],prova[:,1], 
                               prova[:,2],color=col,alpha=1)

                if not aggregate:
                    plt.tight_layout()
                    plt.show()
                
        if aggregate:        
            plt.tight_layout()
            plt.show()
        

    def plot_results(self, min_len = 5, k=0, aggregate = False):   
        
        if sum([len(c) for c in self.CLUSTERS])>0:
            if self.verbose:
                print('Clusters')            
            self.plot_paths_clustered(min_len, k = k, aggregate = aggregate)
        if sum([len(c) for c in self.DEADLOCKS])>0:
            if self.verbose:
                print('Deadlocks Not Working!')
#            self.plot_paths_clustered(self.DEADLOCKS, min_len = -1,k=k,aggregate = aggregate)
    
    def network_plot_3D(self, G, fun):

        # Get node positions
        pos = nx.get_node_attributes(G, 'pos')
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        V = np.array([pos[key] for key in pos.keys()])    

        fig = plt.figure(figsize=(10,5))
        ax = fig.add_subplot(projection='3d')
        plt.xlim([self.m_x, self.M_x])
        plt.ylim([self.m_y, self.M_y])
        ax.set_zlim([self.m_z, self.M_z])
        
        ax.scatter(V[:,0],V[:,1],V[:,2], c=fun)

        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        for i,j in enumerate(G.edges()):

            x = np.array((pos[j[0]][0], pos[j[1]][0]))
            y = np.array((pos[j[0]][1], pos[j[1]][1]))
            z = np.array((pos[j[0]][2], pos[j[1]][2]))

        # Plot the connecting lines
            ax.plot(x, y, z, c='black', alpha=0.5)  

        plt.tight_layout()
        plt.show()
        
    def plot_triangulation_and_paths(self,LIST,axis=-1, k=0.1):
        
        D_flat = squareform(pdist(self.res_grid))
        colors = ['r','g','k']

        for h in [self.x,self.y,self.z][axis]:

            for edge in self.edge_idxs.keys():
                if np.prod([self.res_grid[p][axis]==h for p in edge])==1:

                    e = np.array(list(edge))
                    if np.max(D_flat[e,:][:,e])>self.r_graph_aux:
                        alpha=0.1
                    else:
                        alpha=1
                    plt.plot(self.res_grid[e][:,(axis+1)%3],self.res_grid[e][:,(axis+2)%3], c='k', alpha=alpha,zorder=2)


            for tr in self.ITRIS:
                if np.prod([self.res_grid[p][axis]==h for p in tr])==1:
                    if np.max(D_flat[tr,:][:,tr])>self.r_graph_aux:
                        alpha=0.05
                    else:
                        alpha=1

                    plt.fill(self.res_grid[tr][:,(axis+1)%3],self.res_grid[tr][:,(axis+2)%3],c='b',alpha=alpha,zorder=1)    
            
            for i,idx in enumerate(LIST):
                p = self.PATHS[idx]
                coords = self.res_grid[p,:]
                lvl = coords[coords[:,axis]==h]   
                noise = np.random.normal(size=(len(lvl),3))*k
                lvl=lvl+noise
                print('Path ', i,'has ',len(lvl),' points at height ',h,'.')
                
                plt.scatter(lvl[:,0],lvl[:,1], c=colors[i],zorder=3)
                
            plt.show()


    """
    Commands to Run the Diffusion Process and Find Homologous Cycles
    """    

    def check_lithium_ranges(self,plot=True):
        
        axes_aux = [0,1,2]
        
        max_range = np.sqrt(3)*(self.M[0]-self.m[0])

        min_ = []
        max_ = []

        for i in np.arange(self.Li.shape[1]):
            
            subsample = np.arange(0,self.Li.shape[0],1000)
            
            if self.verbose:
                print('Doing atom ',i,' out of ',self.Li.shape[1], end='\r')
                
            p = self.Li[:,i,:]
            p = p[subsample,:]
            D = dist_from_pts_periodic_boundaries(p,p,self.M,self.m,axes_aux)            
            min_.append(0)
            max_.append(np.max(D))

        min_= np.array(min_)
        max_ = np.array(max_)
        range_ = np.abs(max_-min_)/max_range
            
        if self.verbose:
            print('Portion of the box travelled by Lithium atoms: ',np.sort(range_)[::-1][:10])

        if plot:
            plt.boxplot(range_)
            plt.show()
                
    def check_lithium_paths(self,n=5):
        
        n_Li = self.Li.shape[1]

        for _ in range(n):
            idx = np.random.multinomial(1,np.ones((n_Li,))/n_Li)
            idx = np.argwhere(idx>0)[0][0]
            self.proj_lithium_paths([idx])
    
    def reduce_boundary_mat(self,):

        self.GF = galois.GF(2)

        B = self.B.T.astype(float)
        I = np.identity(B.shape[0])

        self.B = self.GF(B.astype(int)%2)
        I = self.GF(I.astype(int)%2)

        BI = np.concatenate((self.B, I), axis=-1)
        Q = BI.row_reduce()
        self.V = Q[:,-B.shape[0]:]
        self.B_red = Q[:,:B.shape[1]] 
        self.pivot_cols = np.where(np.sum(self.B_red==1,axis=0)==1)[0]
        self.rank = np.sum(np.max(self.B_red, axis=-1)>0)
#        self.V = self.V[self.rank:,:]

        self.A = self.GF(self.A.astype(int)%2)        
        A_ = self.V@self.A
        A_ = A_[self.rank:,:]
        
        I = np.identity(A_.shape[0])
        I = self.GF(I.astype(int)%2)

        AI = np.concatenate((A_, I), axis=-1)
        Q = AI.row_reduce()
        self.V_ = Q[:,-A_.shape[0]:]
        self.A_red = Q[:,:A_.shape[1]] 
        self.rank_A = np.sum(np.max(self.A_red, axis=-1)>0)        
        
    def compute_vector_field_penalties(self, v_ = np.array([0,0,1])):
    
        if len(v_.shape)<2:
            self.v = np.ones_like(self.res_grid)*v_
        else:
            self.v = v_

        self.grad_pen = np.array([grad_r2(pt, self.balls_centres, self.balls_radii) for pt in self.res_grid])

        aux = np.isnan(self.grad_pen)
        self.grad_pen[aux]=0*self.grad_pen[aux]
        
        self.direction = np.zeros_like(self.D)
        self.repulsion = np.zeros_like(self.D)

        for i,p in enumerate(self.res_grid):
            if self.verbose:
                print(i,len(self.res_grid), "            ", end="\r")
            for j in range(i):
                path, lengths = find_geod(i,j,self.geod, self.res_grid, self.D)
                lengths = lengths/self.D[i,j]

                line = np.diff(self.res_grid[path], axis=0)   

                self.direction[i,j]= np.sum(np.sum(-self.v[path[:-1]]*(line),axis=1)*lengths)
                self.direction[j,i] = - self.direction[i,j]

                self.repulsion[i,j]= np.sum(np.sum(-self.grad_pen[path[:-1]]*(line),axis=1)*lengths)
                self.repulsion[j,i] = - self.repulsion[i,j]

    def make_markov_chain_matrix(self, eps_ = 1):

        self.eps = eps_                               
        self.N = self.make_N()
        norm = np.sum(self.N,axis=1)

        for i in range(self.N.shape[0]):
            self.N[i,:] = self.N[i,:]/norm[i]


    def equivalence_wrapper(self, loop):
        
        v_path = path_to_vec(loop,self.D,self.N_1,self.edge_idxs,self.r_graph)

        b_ = self.GF(v_path.astype(int)%2)
        b = self.V[self.rank:,:]@b_

        return np.max(b)        
            
    def run_A_0_to_A_1_diffusion(self, n_paths=200, 
                                 save = False,
                                 folder_path = '/Users/ye72al/Library/Mobile Documents/com~apple~CloudDocs/Glass_TDA/'):

        self.DEADLOCKS = []
        self.PATHS = []        
        cnt = 0
        
        while cnt<n_paths:

            if self.verbose:
                print('Paths Sampled:', len(self.PATHS), 
                      '- Deadlocks Found:', len(self.DEADLOCKS), '          ', end = '\r')

            tmp = []
            
            while len(tmp)==0:

                path,tmp = self.sample_diffusion()
                p = self.preprocess_path(path)

                if len(tmp)==0:
                    self.DEADLOCKS.append(p)
                else:
                    self.PATHS.append(p)
                    cnt = len(self.PATHS)

        if save:
            if sum([len(c) for c in self.CLUSTERS])>0:
                for i,c in enumerate(self.CLUSTERS):  
                    LIST = [self.PATHS[i] for i in c]
                    save_paths(LIST,'CLUSTER_'+str(i),folder_path)
            
            if sum([len(c) for c in self.DEADLOCKS])>0:
                print('Deadlocks not Working!')
#                for i,c in enumerate(self.DEADLOCKS):    
#                    LIST = [self.PATHS[i] for i in c]
#                    save_paths(LIST,'DEADLOCKS_'+str(i),folder_path)    

    def random_simplification_wrapper_MP(self,LIST):
        
        path,n,i = LIST
        
        if self.verbose:
            print('Simplifying Path ',i)
        
        p = self.random_simplification_wrapper(path,n)
        
        return p
        


    def random_simplification_wrapper(self,path,n=20,verbose=False):
        
        if self.cut_paths_until_A_0:
            path_ = self.cut_until_A_0(path)
        else:
            path_ = np.copy(path)

#        d = np.sum([self.D[idx,idx+1] for idx in path_[:-1]])
        cnt = 0 
        aux = 0
        
        while cnt < n and aux < 100000:

       #     d_aux = np.sum([self.D[path_[i],path_[i+1]] for i in range(len(path_)-1)])
            
       #     if verbose:
       #         print('Tacchino ',d_aux,cnt)
            
            
            aux+=1
            n_ = len(path_)
            idx0 = np.random.multinomial(1,np.ones((n_,))/n_)
            idx0 = np.argwhere(idx0>0)[0][0]
                                         
            idx1 = np.random.multinomial(1,np.ones((n_,))/n_)
            idx1 = np.argwhere(idx1>0)[0][0]
            
            [p0,p1] = np.sort([idx0,idx1])

            q0 = path_[p0]
            q1 = path_[p1]

            tmp0 = path_[p0:p1+1]
            d_tmp = np.sum([self.D[tmp0[i],tmp0[i+1]] for i in range(len(tmp0)-1)])

            if self.D[q0,q1]+0.0000001<d_tmp:
                tmp1 = np.array([q0])
                loop = close_paths(tmp0,tmp1,self.res_grid,self.geod,self.D)                
                v = path_to_vec(loop, self.D, self.N_1, self.edge_idxs, self.r_graph)
                res = self.prova([loop], BOUNDARY_MAT = 0, RED_MAT = None, verbose = False)
                
                if res == 0:
                    tmp, L = find_geod(q0, q1, self.geod, self.res_grid, self.D) 
                                        
                    path_ = np.array(list(path_[:p0]) + list(tmp) + list(path_[p1+1:]))
                    path_ = simplify_path(path_)                    
                    
                    dd_ = np.sum([self.D[path_[i],path_[i+1]] for i in range(len(path_)-1)])
                    
#                    if dd_ > d_aux + 0.001:
#                        mario
                        
                    cnt = 0
                else:
                    cnt+=1
            else:
                cnt+=1
                
        if self.cut_paths_until_A_0:
            path_ = self.cut_until_A_0(path_)
                                         
        return path_
                                         

    def simplify_paths_with_homology(self,n=100,random=False):
        
        if random:
            if self.MP:
                    pool = mp.Pool(processes=10)
                    aux=pool.map(self.random_simplification_wrapper_MP,((p,n,i) for i,p in enumerate(self.PATHS)))
                    pool.close()
            else:     
                if self.verbose:
                    PATHS_ = []
                    for p in self.PATHS:
                        print('Simplified ',len(PATHS_), ' out of ', len(self.PATHS),'            ',end='\r')
                        PATHS_.append(self.random_simplification_wrapper(p,n))
                    self.PATHS = PATHS_

                else:

                    self.PATHS = [self.random_simplification_wrapper(p,n) for p in self.PATHS]                        
        else:       
            if self.verbose:
                PATHS_ = []
                for p in self.PATHS:
                    print('Simplified ',len(PATHS_), ' out of ', len(self.PATHS),'            ',end='\r')
                    PATHS_.append(self.simplify_path_with_homology(p))
                self.PATHS = PATHS_

            else:
                self.PATHS = [self.simplify_path_with_homology(p) for p in self.PATHS]
                                   
                
    def simplify_deadlocks_with_homology(self,):
        
        if self.verbose:
            DEADLOCKS_ = []
            for p in self.DEADLOCKS:
                print('Simplified ',len(DEADLOCKS_), ' out of ', len(self.DEADLOCKS))
                DEADLOCKS_.append(self.simplify_path_with_homology(p))
            self.DEADLOCKS = DEADLOCKS_
            
        else:
            self.DEADLOCKS = [self.simplify_path_with_homology(p) for p in self.DEADLOCKS_]
                
    def cluster_paths(self, save = False):
        
        self.CLUSTERS = [[0]]
          
        for i_,p in enumerate(self.PATHS[1:]):
            i = i_+1
            
            if self.verbose:
                print('Paths Clustered:', i,'/',len(self.PATHS), 
                      '- Num of Clusters Found:', len(self.CLUSTERS), 
                      '- Clust Cardinalities:',
                      np.sort([len(c) for c in self.CLUSTERS])[::-1][:10],'       ', end='\r')
                
            new = 0
            
            for j,c in enumerate(self.CLUSTERS):
                
                q = self.PATHS[c[0]]
                aux = self.prova([p,q], BOUNDARY_MAT = None, RED_MAT = None, verbose = False)
                
                if aux == 0:
                    new = 1
                    break

            if new == 0:    
                self.CLUSTERS.append([i])
            else:
                self.CLUSTERS[j].append(i)
                
                
        if save:
            self.save_clusters()
            
    def reeb_graph_voids(self,axis=-1, plot = True):
        
        E, fun, emb = self.reeb_graph(self.res_grid,axis=axis, D_ = self.D_force_graph)            
        G = self.build_graph(E,emb,fun)
        
        if plot:
            self.network_plot_3D(G,fun)
               
        
    def reeb_graph_balls(self,axis=-1, plot = True):
        
        E, fun, emb = self.reeb_graph(self.balls,axis=axis)            
        G = self.build_graph(E,emb,fun)
        
        if plot:
            self.network_plot_3D(G,fun)        
         
        
    def prova(self, paths, BOUNDARY_MAT = None, RED_MAT = None, V_MAT = None,verbose = False):
        """
        Check se una lista di paths paths sono riempibili con una 2-chain a meno di un sottospazio generato 
        da una BOUNDARY_MAT. Si sovrappone con fill_loop.
        BOUNDARY_MAT -> la matrice A, che genera il sottospazio di bordo. Se é "None" prendo A_0 e A_1.
        RED_MAT -> la matrice A_red cioé la versione row-reduced di A[rank(B):,:]
        V_MAT -> la matrice che porta A[rank(B):,:] in forma A_red
        """
        path = np.concatenate(paths)
        
        v = path_to_vec(path, self.D, self.N_1, self.edge_idxs, self.r_graph, verbose = verbose)
        v = self.GF(v.astype(int)%2)
        v_ = self.V@v
        b = v_[self.rank:]
            
        if BOUNDARY_MAT is None:
            b_ = self.V_[self.rank_A:,:]@b
        elif np.max(BOUNDARY_MAT)==0:
            b_ = self.V[self.rank:,:]@v
        else:
            if RED_MAT is None:
                A = self.V@BOUNDARY_MAT
                A_red, V = row_reduce_matrix(A[self.rank:,:],self.GF)
            else:
                A = BOUNDARY_MAT
                V = RED_MAT
                A_red = V@A[self.rank:,:]
                
            rank = np.sum(np.max(A_red, axis=-1)>0)
            b_ = V[rank:,:]@b
            
        tmp = np.max(b_)
                    
        return tmp
#        if tmp == 0:
#            print('Il path si puó shrinkare')
#        else:
#            print('Il path NON si puó shrinkare')


    def make_boundary_matrix_from_paths(self, BOUNDARY_PATHS, add_top_bottom = False):

        A = np.zeros((self.N_1,sum([len(p) for p in BOUNDARY_PATHS])))

        for n_col,path_ in enumerate(BOUNDARY_PATHS):
            for i,p_ in enumerate(path_[:-1]):
                q_ = path_[i+1]
                n_row = self.edge_idxs[tuple(np.sort([p_,q_]))]
                A[n_row,n_col] += 1

        A = self.GF(A.astype(int)%2)

        if add_top_bottom:
            A = np.concatenate([self.A,A],axis=-1)

        return A


    def make_overlapping_summary(self, paths, step=1, plot=True, plot_steps = False):

        overl = []
        
        for i,z in enumerate(self.z[:-step]):
            z_l = z
            z_u = self.z[i+step]
            
            axis = -1
            idxs_e = []

            for j,e in enumerate(self.idxs_to_edge):
                edge_z = self.res_grid[e,axis]
#                tmp0 = np.prod([edge_z == z_l])
#                tmp1 = np.prod([edge_z == z_u])
                tmp0 = np.linalg.norm(edge_z-z_l)<0.000001
                tmp1 = np.linalg.norm(edge_z-z_u)<0.000001

                if np.max([tmp0, tmp1])>0:
                    idxs_e.append(j)        

            idxs_e = np.array(idxs_e)

            A = np.zeros((self.N_1,len(idxs_e)))

            for j,idx in enumerate(idxs_e):
                A[idx,j]+=1
                                
            A = self.GF(A.astype(int)%2)

            paths_tmp = [cut_path_from_bottom(p, z_l, self.res_grid, axis=-1) for p in paths]            
            paths_tmp = [cut_path_from_top(p, z_u, self.res_grid, axis=-1) for p in paths_tmp]
            
            if plot_steps:
                _, tmp = self.fill_loop(paths_tmp, BOUNDARY_MAT = A, RED_MAT = None, V_MAT=None, plot = True)
            else:
                tmp = self.prova(paths_tmp, BOUNDARY_MAT = A, RED_MAT = None, V_MAT = None, verbose = False)
            overl.append(tmp)
            
            if self.verbose:
                print('Overlap Summary ',overl)

        overl = np.array(overl)
        
        if self.verbose:
            print('Overlap Summary ',overl)

        if plot:
            
            fig = plt.figure(figsize=(15,10))

            ax1 = fig.add_subplot(121,projection='3d')
            ax2 = fig.add_subplot(122,projection='3d')

           
            idxs = np.where(overl==0)[0]
            
            for idx in idxs:
                z_l = self.z[idx]
                z_u = self.z[idx+step]
                self.fill_loop_on_slice(paths, z_l,z_u,plot=True, FIG = [fig,ax1,ax2])
            plt.show()
            
            plt.plot(self.z[:-step],overl)
            plt.show()
            
        return overl


    def fill_loop_on_slice(self,paths, z_l, z_u, plot=True, FIG = None):
            
        axis = -1
        idxs_e = []

        for j,e in enumerate(self.idxs_to_edge):
            edge_z = self.res_grid[e,axis]
            tmp0 = np.linalg.norm(edge_z-z_l)<0.000001
            tmp1 = np.linalg.norm(edge_z-z_u)<0.000001
            
            if np.max([tmp0, tmp1])>0:
                idxs_e.append(j)        

        idxs_e = np.array(idxs_e)

        A = np.zeros((self.N_1,len(idxs_e)))

        for j,idx in enumerate(idxs_e):
            A[idx,j]+=1

        A = self.GF(A.astype(int)%2)

        paths_tmp = [cut_path_from_bottom(p, z_l, self.res_grid, axis=-1) for p in paths]            
        paths_tmp = [cut_path_from_top(p, z_u, self.res_grid, axis=-1) for p in paths_tmp]

        if not plot:
            tmp = self.prova(paths_tmp, BOUNDARY_MAT = A, RED_MAT = None, verbose = False)
        else:
            
            chain, tmp = self.fill_loop(paths_tmp, BOUNDARY_MAT = A, RED_MAT = None, V_MAT = None,
                                        find_chain = True, plot = False)
            boundary = self.B@chain
            res_grid = self.traslate_box(self.res_grid)
            
            if FIG is None:
                fig = plt.figure(figsize=(15,10))

                ax1 = fig.add_subplot(121,projection='3d')
                ax2 = fig.add_subplot(122,projection='3d')
            else:
                fig, ax1, ax2 = FIG

            ax1.set_xlim([self.m_x, self.M_x])
            ax1.set_ylim([self.m_y, self.M_y])
            ax1.set_zlim([self.m_z, self.M_z])

            ax2.set_xlim([self.m_x, self.M_x])
            ax2.set_ylim([self.m_y, self.M_y])
            ax2.set_zlim([self.m_z, self.M_z])

            ax1.scatter(self.res_grid[:,0],self.res_grid[:,1], self.res_grid[:,2],c='k',alpha=0)
            ax2.scatter(res_grid[:,0],res_grid[:,1], res_grid[:,2],c='k',alpha=0)

            """
            Plot the Paths
            """
            loop = np.concatenate(paths,axis=-1,dtype=int)
            v = path_to_vec(loop, self.D, self.N_1, self.edge_idxs, self.r_graph, verbose = False)

            for idx in np.where(v>0)[0]:
                e = self.idxs_to_edge[idx,:]
                ax1.plot3D(self.res_grid[e][:,0], self.res_grid[e][:,1], self.res_grid[e][:,2], 
                           linewidth=2, c='k') 

                aux = self.traslate_box(self.res_grid[e])
                ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], linewidth=2, c='k') 

            """
            Plot the Paths Slices Which are Homologous
            """                
            loop = np.concatenate(paths_tmp,axis=-1,dtype=int)
            v = path_to_vec(loop, self.D, self.N_1, self.edge_idxs, self.r_graph, verbose = False)

            for tr in self.ITRIS[np.where(chain>0)[0]]:

#                    vtx = np.vstack([self.res_grid[tr,:],[self.res_grid[tr[0],:]] ])                      
                vtx = self.res_grid[tr,:]                    
                tri = Poly3DCollection([vtx], alpha=0.3)
                tri.set_facecolor('b')
                tri.set_edgecolor('k')
                ax1.add_collection3d(tri)

                vtx = self.traslate_box(vtx)                   
                tri = Poly3DCollection([vtx], alpha=0.3)
                tri.set_facecolor('b')
                tri.set_edgecolor('k')
                ax2.add_collection3d(tri)
            
            for p in paths:
                ax1.plot3D(self.res_grid[p,0], self.res_grid[p,1], self.res_grid[p,2], 
                               linewidth=1,linestyle='dotted', c='red') 

                aux = self.traslate_box(self.res_grid[p,:])
                ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], 
                               linewidth=1,linestyle='dotted', c='red') 
            
            for idx in np.where(v>0)[0]:
                e = self.idxs_to_edge[idx,:]
                ax1.plot3D(self.res_grid[e][:,0], self.res_grid[e][:,1], self.res_grid[e][:,2], 
                           linewidth=6, c='red') 

                aux = self.traslate_box(self.res_grid[e])
                ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], linewidth=6, c='red') 
                
            for idx in np.where(boundary>0)[0]:
                e = self.idxs_to_edge[idx,:]
                ax1.plot3D(self.res_grid[e][:,0], self.res_grid[e][:,1], self.res_grid[e][:,2], 
                           linewidth=3, c='green') 

                aux = self.traslate_box(self.res_grid[e])
                ax2.plot3D(aux[:,0], aux[:,1], aux[:,2], linewidth=3, c='green') 

                
            if FIG is None:
                plt.show()
        
        
        return tmp
        

        
        
                  
    