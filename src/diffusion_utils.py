##
# @file diffusion_utils.py
# @brief Utility functions for distance computation, path manipulation, and grid operations with periodic boundaries.
# @version 1.3.0
# @date March 2026

import numpy as np
import galois
import matplotlib.pyplot as plt
import ase.io as io

from scipy.spatial.distance import pdist, squareform, cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra, connected_components

from itertools import chain, combinations, product

import pandas as pd

"""
Utils
"""

def row_reduce_matrix(A,GF):
    """! @brief Row reduce a matrix over a Galois field.

    Computes the row-reduced echelon form of matrix A over the given Galois field,
    along with the transformation matrix V such that V*A = A_red.

    @param A  The matrix to row reduce (shape: (n, m)).
    @param GF  The Galois field to perform arithmetic in.
    @return Tuple (A_red, V) where A_red is the row-reduced matrix and V is the
            transformation matrix.
    """

    I = np.identity(A.shape[0])
    I = GF(I.astype(int)%2)

    AI = np.concatenate((A, I), axis=-1)
    Q = AI.row_reduce()
    V = Q[:,-A.shape[0]:]
    A_red = Q[:,:A.shape[1]] 

    return A_red, V 


def solve_linear_system_with_boundary(B, v, A, A_red, V, GF):
    """! @brief Solve a linear system with boundary constraints.

    Finds a solution, if it exists, of B*x = v + A*mu where A generates
    the boundary subspace.

    @param B      Row-reduced coefficient matrix.
    @param v      Right-hand side vector.
    @param A      Matrix generating the boundary subspace.
    @param A_red  A[rank(B):,:] in row-reduced form via the matrix V.
    @param V      Transformation matrix from row reduction of A.
    @param GF     Galois field for arithmetic operations.
    @return Tuple (x, tmp) where x is the solution vector and tmp indicates
            whether the solution is exact (0) or approximate (1).
    """

    rank = np.sum(np.max(B, axis=-1)>0)
    b = v[rank:]
    b_ = V@b
    mu, tmp0 = solve_linear_system(A_red,b_)
    mu = GF(mu.astype(int)%2)
    v_ = v + A @ mu 
    
    x, tmp1 = solve_linear_system(B,v_)
    x = GF(x.astype(int)%2)
    
    tmp = np.max([tmp0,tmp1])
    
    return x, tmp

def solve_linear_system(B,v):
    """! @brief Solve a linear system B*x = v where B is row-reduced.

    Finds a solution, if it exists, of B*x = v using back-substitution
    on the row-reduced matrix B.

    @param B  Row-reduced coefficient matrix (shape: (n, m)).
    @param v  Right-hand side vector (shape: (n,)).
    @return Tuple (x, tmp) where x is the solution vector (shape: (m,)) and
            tmp indicates if the system is consistent (0) or inconsistent (1).
    """
    
    non_zero = np.where(v==1) 
    pivot_cols = np.where(np.sum(B==1,axis=0)==1)[0]
    rank = np.sum(np.max(B, axis=-1)>0)
    
    x = np.zeros((B.shape[1],), dtype=int)
    x[pivot_cols] = 1
    x[pivot_cols] = np.multiply(x[pivot_cols],v[:rank])
    
    tmp = np.max(v[rank:])
    
    return x, tmp
    
    
def cut_path_from_bottom(path, h, grid, axis=-1):
    """! @brief Trim a path from the bottom at a height threshold.

    Removes leading points from the path whose coordinate along the specified
    axis is at or below the threshold h.

    @param path  Array of grid point indices representing the path.
    @param h     Height threshold value.
    @param grid  Array of grid point coordinates (shape: (npts, dim)).
    @param axis  Coordinate axis to compare against the threshold (default: -1).
    @return Trimmed path as an integer numpy array.
    """

    start_idx = 0
    path_ = np.copy(path)
     
    while grid[path_[start_idx+1]][axis]<=h*1.00001:
        path_ = path_[1:]

    return np.array(path_,dtype=int)

    
def cut_path_from_top(path, h, grid, axis=-1):
    """! @brief Trim a path from the top at a height threshold.

    Removes trailing points from the path whose coordinate along the specified
    axis is at or above the threshold h.

    @param path  Array of grid point indices representing the path.
    @param h     Height threshold value.
    @param grid  Array of grid point coordinates (shape: (npts, dim)).
    @param axis  Coordinate axis to compare against the threshold (default: -1).
    @return Trimmed path as an integer numpy array.
    """

    start_idx = -1
    path_ = np.copy(path)
    
    while grid[path_[start_idx-1]][axis]>=h*0.999999:
        path_ = path_[:-1]

    return np.array(path_,dtype=int)


def get_cmap(n, name='hsv'):
    """! @brief Get a colormap for distinct color assignment.

    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color using a matplotlib colormap.

    @param n     Number of distinct colors needed.
    @param name  Name of the matplotlib colormap to use (default: 'hsv').
    @return A colormap function mapping integer indices to RGBA tuples.
    """
    return plt.cm.get_cmap(name, n)


def get_box_from_cell(cell_file):
    """! @brief Read simulation box dimensions from a cell file.

    Parses a CSV cell file to extract the upper and lower bounds of the
    cubic simulation box.

    @param cell_file  Path to the cell file (CSV format with whitespace delimiter).
    @return Tuple (M, m) where M is the upper bound array [Mx, My, Mz] and
            m is the lower bound array [mx, my, mz].
    """

    cell = pd.read_csv(cell_file, delimiter = '       ')

    M_x = cell['Bx [Angstrom]'].max()
    M_y = M_x
    M_z = M_x
    
    m_x,m_y,m_z = 0,0,0
    
    return np.array([M_x,M_y,M_z]), np.array([m_x,m_y,m_z])

### DEPRECATED?
# def read_data(path, backbone_atoms, flow_atoms):
#     """
#     Read an xyz file and extract the molecule's geometry.

#     The file should be in the format::

#         num_atoms
#         comment
#         ele0 x0 y0 z0
#         ele1 x1 y1 z1
#         ...

#     Parameters
#     ----------
#     path : str
#         A path to a file to read

#     Returns
#     -------
#     val : LazyValues
#         An object storing all the data
#     """
#     N_ATOMS = []
#     ATOMS = []
#     COORDS = []
#     TIMESTEPS = []

#     aux = -1
#     elements = []
#     coords = []
#     comments={}
#     M = []
#     m = []
    
#     tot = 100000000
    
#     with open(path, 'r') as f:
#         for i, line in enumerate(f):  
       
#             if i%tot==1:

#                 time_step = line.strip().split()[0]
#                 TIMESTEPS.append(time_step)
            
#             elif i%tot==3:
                
#                 n_atom = line.strip().split()[0]
                
#                 N_ATOMS.append(float(n_atom))
#                 ATOMS.append(elements)
#                 COORDS.append(coords)
    
#                 elements = []
#                 coords = []
        
#                 aux += 1
#                 tot = float(n_atom)+9
                
#             elif i%tot==5:
    
#                 box = line.strip().split() 
#                 M_x = float(box[1])
#                 m_x = float(box[0])
                                  
#             elif i%tot==6:
                
#                 box = line.strip().split() 
#                 M_y = float(box[1])
#                 m_y = float(box[0])
 
#             elif i%tot==7:
                
#                 box = line.strip().split() 
#                 M_z = float(box[1])
#                 m_z = float(box[0])
                    
#             elif i%tot==0 or i%tot==2 or i%tot==4 or i%tot==8:
            
#                 pass

#             else:
                    
#                 if len(line.strip().split())==5:

#                     _,type_, x, y, z = line.strip().split()

#                     if float(type_) == 8:
#                         ele  = 'Li'
#                     elif float(type_) == 2 or float(type_) == 4 or float(type_) == 5 or float(type_) == 7:
#                         ele = 'S'
#                     else:
#                         ele = 'P'
                
#                 elif len(line.strip().split())==6:
                         
#                     _,_,ele, x, y, z = line.strip().split()

#                 else:
                    
#                     mario

#                 point = (float(x), float(y), float(z))

#                 elements.append(ele)
#                 coords.append(point)                
                    
#         ATOMS.append(elements)
#         COORDS.append(coords)
          
#     M,m = np.array([M_x,M_y,M_z]), np.array([m_x,m_y,m_z])
                               
#     return N_ATOMS, ATOMS[1:], COORDS[1:], TIMESTEPS, M, m


def out_to_atoms(inputfile, backbone_atoms, flow_atoms):
    """! @brief Read an XYZ file and return atom coordinates by type.

    Parses the input XYZ file and separates atom coordinates into backbone
    and flow atom arrays, along with cardinalities and box boundaries.

    @param inputfile      Path to the .xyz input file.
    @param backbone_atoms List of element symbols for backbone atoms.
    @param flow_atoms     List of element symbols for flow (mobile) atoms.
    @return Tuple (Li, P, S, N, TIMESTEPS, M, m) containing atom coordinate
            arrays, atom counts, timestep list, and box boundary arrays.
    """
#    N_ATOMS, ATOMS, COORDS, COMMENTS, M, m = read_xyz_data(inputfile, box)
    N_ATOMS, ATOMS, COORDS, TIMESTEPS, M, m = read_data(inputfile)

    n_backbone = [0 for a in backbone_atoms]
    n_flow = [0 for a in flow_atoms]
    
    for atom in ATOMS[0]:
        if atom in backbone_atoms:
            n_backbone[backbone_atoms.index(atom)] += 1
        elif atom in flow_atoms:
            n_flow[flow_atoms.index(atom)] += 1
            
    backbone = np.zeros((len(N_ATOMS),sum(n_backbone),3))
    flow = np.zeros((len(N_ATOMS),sum(n_flow),3))
    
    for i,n in enumerate(N_ATOMS):        
        
        n_Li, n_P, n_S = -1,-1,-1
        
        for j,atom in enumerate(ATOMS[i]):
            
            if atom=='Li':
                n_Li += 1
                Li[i,n_Li,:]=COORDS[i][j]
            elif atom=='P':
                n_P += 1
                P[i,n_P,:]=COORDS[i][j]
            else:
                n_S += 1
                S[i,n_S,:]=COORDS[i][j]
 
    N = np.array([n_Li,n_P,n_S])+1
    
    return Li,P,S,N,TIMESTEPS,M,m 


def clean_Li(Li,m,M):
    """! @brief Remove atoms outside the simulation box.

    Filters out lithium atoms whose coordinates fall outside the
    simulation box defined by lower bound m and upper bound M.

    @param Li  Lithium atom coordinates array (shape: (timesteps, n_atoms, 3)).
    @param m   Lower bounds of the simulation box (scalar or array of shape (3,)).
    @param M   Upper bounds of the simulation box (scalar or array of shape (3,)).
    @return Filtered lithium coordinates array with out-of-box atoms removed.
    """
    
    idxs = []
    
    if len(M)==1:
        M_x,M_y,M_z = M,M,M
        m_x,m_y,m_z = m,m,m
    else:
        M_x,M_y,M_z = M
        m_x,m_y,m_z = m

    p = np.array([M_x,M_y,M_z])
    P = np.array([m_x,m_y,m_z])
    
    for i in range(Li.shape[1]):
        
        q = Li[:,i,:]
        
        if np.min(P-q)<0 or np.min(q-p)<0:
            idxs.append(i)
            
    if len(idxs)>0:
        idxs= np.array(idxs)    
        new_Li = np.delete(Li, idxs,axis=1)

        return new_Li
    else:
        return Li

def make_geodesics(grid, r, M, m, axes, dim, D_= None, D_force_graph_ = None, verbose = False):
    """! @brief Compute geodesic distances and predecessor matrix on a grid graph.

    Builds a distance graph from the grid points using periodic boundary conditions,
    then computes shortest-path (geodesic) distances and predecessors via Dijkstra's
    algorithm.

    @param grid             Array of grid point coordinates (shape: (npts, dim)).
    @param r                Radius threshold for graph edge connectivity.
    @param M                Upper bounds of the simulation box (shape: (dim,)).
    @param m                Lower bounds of the simulation box (shape: (dim,)).
    @param axes             Array of axes along which periodic boundaries apply.
    @param dim              Dimensionality of the dataset.
    @param D_               Precomputed distance/graph matrix (default: None).
    @param D_force_graph_   Optional mask to force/suppress graph edges (default: None).
    @param verbose          If True, print the number of connected components (default: False).
    @return Tuple (D, geod) where D is the geodesic distance matrix and geod is
            the predecessor matrix for path reconstruction.
    """

    if D_ is None:
        D = dist_periodic_boundaries(grid,M,m,axes, dim)
        
        if D_force_graph_ is None:
            D_force_graph_ = np.ones_like(D_) 
        else:
            D_force_graph = D_force_graph_
            
        graph = D*(D<r)*D_force_graph
    else:
        graph = D_
        
    graph = csr_matrix(graph)
    D, pred = dijkstra(csgraph=graph, directed=False, return_predecessors=True)

    geod = pred.T

    for i in range(geod.shape[0]):
        geod[i,i]=i
    
    if verbose:
        n,_ = connected_components(graph)
        print('Ci sono ',n, ' componenti connesse.')

    
    return D, geod

def find_geod(idx0,idx1,geodesics,pts, D, maxiter = np.inf):
    """! @brief Find the geodesic path between two points.

    Traces the shortest path from idx0 to idx1 using the predecessor matrix,
    collecting coordinates and edge lengths along the way.

    @param idx0       Starting point index.
    @param idx1       Target point index.
    @param geodesics  Predecessor matrix from Dijkstra (shape: (npts, npts)).
    @param pts        Array of grid point coordinates (shape: (npts, dim)).
    @param D          Geodesic distance matrix (shape: (npts, npts)).
    @param maxiter    Maximum number of steps before stopping (default: inf).
    @return Tuple (path, lengths) where path is the array of point indices and
            lengths is the array of edge lengths along the path.
    """

    idx = idx0
    path = [idx]
    coords = [pts[idx]]
    lengths = []
    
    while idx!=idx1 and len(path)<maxiter:
        idx = geodesics[idx,idx1]
        path.append(idx)
        coords.append(pts[idx]) 
        lengths.append(D[idx,path[-2]])
        
    if idx!=idx1 and np.isinf(maxiter):
        path.append(idx1)
        coords.append(pts[idx1]) 
        lengths.append(D[idx,path[-2]])

    path = np.array(path)
    coords = np.array(coords)
    lengths = np.array(lengths)

    return path, lengths

def dist_periodic_boundaries(grid_,M,m,axes, dim=3):
    """! @brief Compute pairwise distances with periodic boundary conditions.

    Computes the full pairwise distance matrix between all grid points,
    accounting for periodic boundary conditions along the specified axes.

    @param grid_  Array of grid point coordinates (shape: (npts, dim)).
    @param M      Upper bounds of the simulation box (shape: (dim,)).
    @param m      Lower bounds of the simulation box (shape: (dim,)).
    @param axes   Array of axes along which periodic boundaries apply.
    @param dim    Dimensionality of the dataset (default: 3).
    @return Symmetric distance matrix of shape (npts, npts).
    """    
#    grid = m + np.remainder(grid_-m, M-m)   

    grid = grid_
    D = cdist(grid,grid)
    DELTAS = M-m
    
    c_idxs = np.arange(dim)
    c_idxs = np.array([i not in axes for i in c_idxs]) 
    DELTAS[c_idxs] = 0
    
    vectors = list([[0,-DELTAS[i],DELTAS[i]] for i in range(len(DELTAS))])
    
    for v in set(product(*vectors)):
        grid_new = grid + v
        D_new = cdist(grid,grid_new)
        D = np.minimum(D_new,D)

    return D

def dist_from_pts_periodic_boundaries(pts_,grid_,M,m,axes, dim=3):
    """! @brief Compute distances from points to a grid with periodic boundaries.

    Computes distances from a set of query points to all grid points,
    accounting for periodic boundary conditions along the specified axes.

    @param pts_   Query points to project (shape: (n, dim) or (dim,)).
    @param grid_  Array of grid point coordinates (shape: (npts, dim)).
    @param M      Upper bounds of the simulation box (shape: (dim,)).
    @param m      Lower bounds of the simulation box (shape: (dim,)).
    @param axes   Array of axes along which periodic boundaries apply.
    @param dim    Dimensionality of the dataset (default: 3).
    @return Distance matrix of shape (n, npts).
    """    
    if len(pts_.shape)<2:
        pts = np.array([pts_])
    else:
        pts = pts_
        
#    pts = m + np.remainder(pts-m, M-m)    
#    grid = m + np.remainder(grid_-m, M-m)    

    pts = pts
    grid = grid_    


    D = cdist(pts,grid)
    DELTAS = M-m
    
    c_idxs = np.arange(dim)
    c_idxs = np.array([i not in axes for i in c_idxs]) 
    DELTAS[c_idxs] = 0
    
    vectors = list([[0,-DELTAS[i],DELTAS[i]] for i in range(len(DELTAS))])
    
    for v in product(*vectors):
        if np.min(v)==np.max(v)==0:
            pass
        else:
            grid_new = grid + v
            D_new = cdist(pts,grid_new)
            D = np.minimum(D_new,D)

    return D


def append_tr(ITRIS, B, edge_idxs, ia,ib,ic):
    """! @brief Append a triangle to the triangle list and update the boundary matrix.

    Adds a triangle defined by vertex indices (ia, ib, ic) to ITRIS and sets
    the corresponding edge entries in the boundary matrix B.

    @param ITRIS      List of triangles (each a list of three vertex indices).
    @param B          Boundary matrix mapping triangles to edges.
    @param edge_idxs  Dictionary mapping sorted vertex pairs to edge indices.
    @param ia         First vertex index (must satisfy ia < ib < ic).
    @param ib         Second vertex index.
    @param ic         Third vertex index.
    @return Tuple (ITRIS, B) with the updated triangle list and boundary matrix.
    """
#    ia,ib,ic = np.sort([ia_,ib_,ic_])
    
    ITRIS.append([ia,ib,ic])
    
    e_ab = edge_idxs[(ia,ib)]
    e_bc = edge_idxs[(ib,ic)]
    e_ac = edge_idxs[(ia,ic)]
    
    B[len(ITRIS)-1,e_ab] = 1
    B[len(ITRIS)-1,e_bc] = 1
    B[len(ITRIS)-1,e_ac] = 1
    
    return ITRIS, B


def check_tr(p0,p1,p2, res_grid, grid, cubic_to_res, ITRIS, B, edge_idxs):
    """! @brief Check and append a triangle from cubic grid indices.

    Converts cubic grid indices to resolution grid indices, sorts them,
    and appends the triangle if all vertices are valid.

    @param p0             First cubic grid index tuple.
    @param p1             Second cubic grid index tuple.
    @param p2             Third cubic grid index tuple.
    @param res_grid       Resolution grid point coordinates.
    @param grid           Full grid point coordinates.
    @param cubic_to_res   Mapping from cubic grid indices to resolution grid indices.
    @param ITRIS          List of triangles to update.
    @param B              Boundary matrix mapping triangles to edges.
    @param edge_idxs      Dictionary mapping sorted vertex pairs to edge indices.
    @return Tuple (ITRIS, B) with the updated triangle list and boundary matrix.
    """
    i0 = cubic_to_res[p0]
    i1 = cubic_to_res[p1]
    i2 = cubic_to_res[p2]
    
    i0, i1, i2 = np.sort([i0,i1,i2]) 

    if i0 > -1:        
        ITRIS, B = append_tr(ITRIS, B, edge_idxs, i0,i1,i2) 
    
    
    return ITRIS, B 

def r_squared(px,py,pz,balls_centres,balls_radii):
    """! @brief Compute a repulsion field value at a point.

    Evaluates a repulsion potential at point (px, py, pz) as the sum of
    inverse power contributions from each ball centre.

    @param px             X-coordinate of the query point.
    @param py             Y-coordinate of the query point.
    @param pz             Z-coordinate of the query point.
    @param balls_centres  Array of ball centre coordinates (shape: (n, 3)).
    @param balls_radii    Array of ball radii (shape: (n,)).
    @return Scalar repulsion field value at the query point.
    """
    p=2
    eps = 0.3
    pt = np.array([px,py,pz])
    out = 0
    
    for i,c in enumerate(balls_centres):
        tmp = (balls_radii[i]-np.linalg.norm(pt-c))/eps       
        aux = np.power(tmp,p)
        out += 1/(eps+aux)

    return out

def aux_fn(px,py,pz,xc,yc,zc,r,eps,p=4):
    """! @brief Compute a single-component gradient contribution for the repulsion field.

    Evaluates the partial derivative contribution of one ball along one
    coordinate axis, used internally by grad_r2.

    @param px   Query point coordinate along the differentiation axis.
    @param py   Query point coordinate along the second axis.
    @param pz   Query point coordinate along the third axis.
    @param xc   Ball centre coordinate along the differentiation axis.
    @param yc   Ball centre coordinate along the second axis.
    @param zc   Ball centre coordinate along the third axis.
    @param r    Radius of the ball.
    @param eps  Smoothing parameter controlling field sharpness.
    @param p    Power exponent for the potential (default: 4).
    @return Scalar gradient contribution along the differentiation axis.
    """
    dx = (px-xc)
    dy = (py-yc)
    dz = (pz-zc)
    a = np.power(dx,2)
    b = np.power(dy,2)
    c = np.power(dz,2)
    
    dist = np.sqrt(a+b+c)
    dist_r = np.abs(dist-r)/eps
    
    num=-p*dx*np.power(dist_r,p-1)
    
    d = eps + np.power(dist_r,p) 
    den = np.power(d,2)*(dist/eps)
        
    out = num/den 
        
    return out
    
def grad_r2(pt,balls_centres,balls_radii,p=2,eps=1):
    """! @brief Compute the negative gradient of the repulsion field at a point.

    Evaluates the gradient of the repulsion potential with respect to
    (x, y, z) at the given point, summing contributions from all balls.

    @param pt             Query point coordinates (shape: (3,)).
    @param balls_centres  Array of ball centre coordinates (shape: (n, 3)).
    @param balls_radii    Array of ball radii (shape: (n,)).
    @param p              Power exponent for the potential (default: 2).
    @param eps            Smoothing parameter controlling field sharpness (default: 1).
    @return Negative gradient vector at the query point (shape: (3,)).
    """
    dx = 0
    dy = 0
    dz = 0

    for i,c in enumerate(balls_centres):

        dx += aux_fn(pt[0],pt[1],pt[2],c[0],c[1],c[2],balls_radii[i],eps,p)
        dy += aux_fn(pt[1],pt[0],pt[2],c[1],c[0],c[2],balls_radii[i],eps,p)
        dz += aux_fn(pt[2],pt[1],pt[0],c[2],c[1],c[0],balls_radii[i],eps,p)
        
    return -np.array([dx,dy,dz])

def sample_path(idx, N, res_grid, z_1, maxiter=100000):
    """! @brief Sample a random walk path on the grid.

    Performs a random walk starting from the given index using the transition
    matrix N, stopping when the path reaches a height above z_1 or after
    maxiter steps.

    @param idx      Starting grid point index.
    @param N        Transition probability matrix (shape: (npts, npts)).
    @param res_grid Array of grid point coordinates (shape: (npts, dim)).
    @param z_1      Height threshold along the last axis to terminate the walk.
    @param maxiter  Maximum number of random walk steps (default: 100000).
    @return Array of grid point indices forming the sampled path.
    """
    PATH = [idx]
    N_idx = N[idx,:]

    for i in range(maxiter):

        res = np.random.multinomial(1,N_idx)
        res = np.argwhere(res>0)[0][0]
        PATH.append(res)

        N_idx = N[res,:]
        
        if res_grid[res,-1]>z_1:
            break

    PATH = np.array(PATH)
    
    return PATH


def smooth_path(path,geod,res_grid,D):
    """! @brief Smooth a path using geodesic interpolation.

    Replaces each consecutive pair of points in the path with the geodesic
    path between them, producing a smoother trajectory on the grid.

    @param path      Array of grid point indices forming the path.
    @param geod      Predecessor matrix for geodesic path reconstruction.
    @param res_grid  Array of grid point coordinates (shape: (npts, dim)).
    @param D         Geodesic distance matrix (shape: (npts, npts)).
    @return Smoothed path as an integer numpy array.
    """
    idx = 0
    path_smooth = []

    while idx < len(path)-1:
        tmp,_= find_geod(path[idx],path[idx+1],geod,res_grid,D)
        path_smooth += list(tmp)
        idx+=1
        
    return np.array(path_smooth)

def close_paths(p_0, p_1, grid, geod,D):
    """! @brief Close two open paths into a loop by connecting their endpoints.

    Connects the start and end points of two paths via geodesics to form
    a closed loop: p_1_start -> p_0_start -> ... -> p_0_end -> p_1_end -> ... -> p_1_start.

    @param p_0   First path as an array of grid point indices.
    @param p_1   Second path as an array of grid point indices.
    @param grid  Array of grid point coordinates (shape: (npts, dim)).
    @param geod  Predecessor matrix for geodesic path reconstruction.
    @param D     Geodesic distance matrix (shape: (npts, npts)).
    @return Closed loop as a numpy array of grid point indices.
    """
    start_point_path,_ = find_geod(p_0[0], p_1[0], geod, grid,D) 
    end_point_path,_ = find_geod(p_0[-1], p_1[-1], geod, grid,D)

    path = list(start_point_path)[::-1] + list (p_0) +\
       list(end_point_path) + list(p_1)[::-1] 
    
    return np.array(path)


def homological_simplification_of_path_pts(idx0, idx1, path, grid, geod, V, rank, GF, D, N, edge_idxs, r_graph):
    """! @brief Simplify a path segment between two indices preserving homology class.

    Attempts to replace the sub-path between idx0 and idx1 with a shorter
    geodesic, provided the resulting loop is homologically trivial.

    @param idx0       First index into the path array.
    @param idx1       Second index into the path array.
    @param path       Array of grid point indices forming the path.
    @param grid       Array of grid point coordinates (shape: (npts, dim)).
    @param geod       Predecessor matrix for geodesic path reconstruction.
    @param V          Transformation matrix from row reduction of the boundary.
    @param rank       Rank of the boundary matrix.
    @param GF         Galois field for arithmetic operations.
    @param D          Geodesic distance matrix (shape: (npts, npts)).
    @param N          Number of 1-simplices (edges) in the complex.
    @param edge_idxs  Dictionary mapping sorted vertex pairs to edge indices.
    @param r_graph    Radius threshold for graph edge connectivity.
    @return Simplified path as a numpy array of grid point indices.
    """
    [p0,p1] = np.sort([idx0,idx1])
    
    q0 = path[p0]
    q1 = path[p1]
    
    tmp0 = path[p0:p1+1]
    d = np.sum([D[idx,idx+1] for idx in tmp0[:-1]])
    
    if D[q0,q1]<d:
        tmp1 = np.array([q0])
        loop = close_paths(tmp0,tmp1,grid,geod,D)
        v = path_to_vec(loop, D, N, edge_idxs, r_graph)
        b_ = GF(v.astype(int)%2)
        b = V@b_

        res = np.max(b)

        if res == 0:
            tmp, _ = find_geod(q0, q1, geod, grid, D) 
            path_ = np.array(list(path[:p0]) + list(tmp) + list(path[p1+1:]))
            path_ = simplify_path(path_)
        else:
            path_ = path
    else:
            path_ = path
        
    return np.array(path_)


def homological_simplification_of_path(path, grid, geod, V, rank, GF, D, N, edge_idxs, r_graph):
    """! @brief Iteratively simplify an entire path preserving its homology class.

    Repeatedly scans the path for segments that can be replaced by shorter
    geodesics without changing the homology class, until no further
    simplification is possible.

    @param path       Array of grid point indices forming the path.
    @param grid       Array of grid point coordinates (shape: (npts, dim)).
    @param geod       Predecessor matrix for geodesic path reconstruction.
    @param V          Transformation matrix from row reduction of the boundary.
    @param rank       Rank of the boundary matrix.
    @param GF         Galois field for arithmetic operations.
    @param D          Geodesic distance matrix (shape: (npts, npts)).
    @param N          Number of 1-simplices (edges) in the complex.
    @param edge_idxs  Dictionary mapping sorted vertex pairs to edge indices.
    @param r_graph    Radius threshold for graph edge connectivity.
    @return Simplified path as a numpy array of grid point indices.
    """
    path_ = path
    path_old = path
    n_ = len(path_)-1
    n = len(path_old)
    
    while n_< n:
        
        n = len(path_old)
        n_ = len(path_)
        path_old = path_
        
        i=0 
        while i < len(path_):
            
            idx0 = path_old[i]
            d = np.sum(np.linalg.norm(np.diff(grid[path_old]),axis=-1))
            
            for j_,idx1 in enumerate(reversed(path_old[i+1:])):
                                
                j = len(path_old)-j_
                if D[idx0,idx1]<d:
                    tmp0 = path_old[i:j]
                    tmp1 = np.array([idx0])
                    loop = close_paths(tmp0,tmp1,grid,geod,D)
                    v = path_to_vec(loop, D, N, edge_idxs, r_graph)
                    b_ = GF(v.astype(int)%2)
                    b = V@b_

                    res = np.max(b)

                    if res == 0:
                        tmp, _ = find_geod(idx0, idx1, geod, grid, D) 
                        path_ = np.array(list(path_old[:i]) + list(tmp) + list(path_old[j:]))
                        path_ = simplify_path(path_)
                        break
 #                       path_ = simplify_path(smooth_path(path_, geod, grid))
                    else:
                        d = D[idx0,idx1]
                else:
                    d = D[idx0,idx1]
                    
            path_old = path_
            n_ = len(path_)
            i += 1
        
    return np.array(path_)


def close_paths_A_0_A_1(p_0,p_1, D, A_0,A_1,res_to_A_0, res_to_A_1, A_0_to_res, A_1_to_res, geod_A_0, geod_A_1):
    """! @brief Close two paths into a loop using boundary subspace geodesics.

    Connects the start points via a geodesic on A_0 and the end points via
    a geodesic on A_1, forming a closed loop through both paths.

    @param p_0           First path as an array of resolution grid indices.
    @param p_1           Second path as an array of resolution grid indices.
    @param D             Geodesic distance matrix.
    @param A_0           Bottom boundary grid point coordinates.
    @param A_1           Top boundary grid point coordinates.
    @param res_to_A_0    Mapping from resolution grid indices to A_0 indices.
    @param res_to_A_1    Mapping from resolution grid indices to A_1 indices.
    @param A_0_to_res    Mapping from A_0 indices to resolution grid indices.
    @param A_1_to_res    Mapping from A_1 indices to resolution grid indices.
    @param geod_A_0      Predecessor matrix for geodesics on A_0.
    @param geod_A_1      Predecessor matrix for geodesics on A_1.
    @return Closed loop as a numpy array of resolution grid point indices.
    """
    start_point_path_,_ = find_geod(res_to_A_0[p_0[0]],res_to_A_0[p_1[0]],geod_A_0,A_0,D) 
    end_point_path_,_ = find_geod(res_to_A_1[p_0[-1]],res_to_A_1[p_1[-1]],geod_A_1,A_1,D)

    start_point_path = [A_0_to_res[p] for p in start_point_path_]
    end_point_path = [A_1_to_res[q] for q in end_point_path_]

    path = list(start_point_path)[::-1] + list (p_0) +\
       list(end_point_path) + list(p_1)[::-1] 
    
    return np.array(path)

def join_points_boundary_conditions(p_, q_, M, m, axes, dim, n=20):
    """! @brief Generate candidate line segments between two points under periodic boundaries.

    Creates straight-line interpolations between p_ and all periodic images
    of q_, wrapping coordinates back into the simulation box.

    @param p_    Start point coordinates (shape: (dim,)).
    @param q_    End point coordinates (shape: (dim,)).
    @param M     Upper bounds of the simulation box (shape: (dim,)).
    @param m     Lower bounds of the simulation box (shape: (dim,)).
    @param axes  Array of axes along which periodic boundaries apply.
    @param dim   Dimensionality of the dataset.
    @param n     Number of interpolation samples per segment (default: 20).
    @return List of interpolated path arrays, one per periodic image.
    """
    p = p_*np.ones((n,dim))
    mu = np.linspace(0,1,n) 
    mu = np.array([mu]).T
    paths = []
    DELTAS = M-m
    
    c_idxs = np.arange(dim)
    c_idxs = np.array([i not in axes for i in c_idxs]) 
    DELTAS[c_idxs] = 0
    
#    vectors = list([[0,-DELTAS[i],DELTAS[i]] for i in range(len(DELTAS))])
    vectors = np.array([np.array([0,DELTAS[i]]) for i in range(len(DELTAS))])
    
    
    for v in set(product(*vectors)):
        
        q = q_*np.ones((n,dim)) + v
        line = p*(1-mu) + q*mu
        line = m + np.remainder(line-m + v, M-m)
        paths.append(line)
        
    return paths

def close_paths_boundary_conditions(p_0, p_1,
                                    res_grid, geod,
                                    M, m, D,
                                    axes, dim):
    """! @brief Close two paths into loops using periodic boundary conditions.

    Generates all candidate closed loops by connecting path endpoints through
    periodic boundary interpolations, projecting onto the grid, and smoothing.

    @param p_0       First path as an array of grid point indices.
    @param p_1       Second path as an array of grid point indices.
    @param res_grid  Array of grid point coordinates (shape: (npts, dim)).
    @param geod      Predecessor matrix for geodesic path reconstruction.
    @param M         Upper bounds of the simulation box (shape: (dim,)).
    @param m         Lower bounds of the simulation box (shape: (dim,)).
    @param D         Geodesic distance matrix (shape: (npts, npts)).
    @param axes      Array of axes along which periodic boundaries apply.
    @param dim       Dimensionality of the dataset.
    @return Array of candidate closed loops (dtype: object).
    """
    start_point_paths_ = join_points_boundary_conditions(res_grid[p_0[0]],
                                                        res_grid[p_1[0]], 
                                                        M, m, axes, dim)
    
    start_point_paths = [simplify_path(smooth_path(project_path(p_, res_grid), geod, res_grid, D)) 
                                             for p_ in start_point_paths_]
    
    end_point_paths_ = join_points_boundary_conditions(res_grid[p_0[-1]],
                                                      res_grid[p_1[-1]], 
                                                      M, m, axes, dim)
    
    end_point_paths = [simplify_path(smooth_path(project_path(p_, res_grid), geod, res_grid, D)) 
                                             for p_ in end_point_paths_]
    
    
    loops = []
    
    for p in start_point_paths:        
        for q in end_point_paths:    
            loop = np.array(list(p)[::-1] + list (p_0) +\
               list(q) + list(p_1)[::-1]) 
    
            loops.append(loop)
    
    return np.array(loops, dtype=object)
    

def close_deadlocks(p_0,p_1, res_grid, geod, res_to_A_0, A_0_to_res, geod_A_0,A_0):
    """! @brief Close two deadlocked paths into a loop if monotonicity is preserved.

    Connects start points via a geodesic on A_0 and end points via a geodesic
    on the resolution grid. Returns an empty list if the connecting path
    is not monotone along the height axis.

    @param p_0          First path as an array of grid point indices.
    @param p_1          Second path as an array of grid point indices.
    @param res_grid     Array of grid point coordinates (shape: (npts, dim)).
    @param geod         Predecessor matrix for geodesics on the resolution grid.
    @param res_to_A_0   Mapping from resolution grid indices to A_0 indices.
    @param A_0_to_res   Mapping from A_0 indices to resolution grid indices.
    @param geod_A_0     Predecessor matrix for geodesics on A_0.
    @param A_0          Bottom boundary grid point coordinates.
    @return Closed loop as a numpy array, or empty list if not monotone.
    """
    start_point_path_,_ = find_geod(res_to_A_0[p_0[0]],res_to_A_0[p_1[0]],geod_A_0,A_0) 
    end_point_path_,_ = find_geod(p_0[-1],p_1[-1],geod,res_grid)

    start_point_path = [A_0_to_res[p] for p in start_point_path_]
    end_point_path = end_point_path_
    
    z_coords = res_grid[start_point_path][:,-1]
    deltas = np.diff(z_coords)
    
    if np.sum(deltas>0)>0 and np.sum(deltas<0)>0:
        return []
    else:
        path = list(start_point_path)[::-1] + list (p_0) +\
           list(end_point_path) + list(p_1)[::-1] 

        return np.array(path)


def simplify_path(path_):
    """! @brief Remove consecutive duplicate points from a path.

    Filters the path so that no two adjacent entries are identical,
    preserving the order of traversal.

    @param path_  Array of grid point indices forming the path.
    @return Simplified path with consecutive duplicates removed.
    """
    if len(path_)==0:
        return path_
    
    path = path_[np.insert(np.diff(path_).astype(np.bool), 0, True)]
    
    
#    path =[path_[i] for i in range(len(path_)-1)  if path_[i]!=path_[i+1]]
    
#    if len(path)==0:
#        return np.array(path)

#    if path[-1]!= path_[-1]:
#        path.append(path_[-1])

#    path = np.array(path)

    return path

def path_to_vec(path,D,N_1,edge_idxs,r_graph,verbose=False):
    """! @brief Convert a path to an edge vector representation.

    Maps each consecutive pair of points in the path to an edge index and
    accumulates edge traversal counts into a vector over all edges.

    @param path       Array of grid point indices forming the path.
    @param D          Geodesic distance matrix (shape: (npts, npts)).
    @param N_1        Total number of edges in the simplicial complex.
    @param edge_idxs  Dictionary mapping sorted vertex pairs to edge indices.
    @param r_graph    Radius threshold; edges longer than this are flagged.
    @param verbose    If True, print warnings for discontinuous edges (default: False).
    @return Edge vector of shape (N_1,) with integer traversal counts.
    """
    v_path = np.zeros(N_1, dtype=int)

    for i,p in enumerate(path[:-1]):
        q = path[i+1]
        if p != q:
            e = np.sort([p,q])

            if D[p,q]>r_graph:
                if verbose:
                    print('Discontinuous Path', p,q,D[p,q],r_graph)
                else:
                    pass
            else:
                edge = edge_idxs[(e[0],e[1])]
                v_path[edge]+=1 
            
    return v_path

def equivalence(A,v_path,rank):
    """
    Broken!
    """
    A=B.T.astype(float)
    b_2 = GF(v_path.astype(int)%2)
    bb = P@b_2

    aux = np.max(bb[rank:])
 
    return aux

def project_path(path_coords, grid):
    """! @brief Project path coordinates onto the nearest grid points.

    For each point in path_coords, finds the closest grid point by
    Euclidean distance and returns the array of grid indices.

    @param path_coords  Array of coordinates to project (shape: (n, dim)).
    @param grid         Array of grid point coordinates (shape: (npts, dim)).
    @return Array of grid point indices (shape: (n,)).
    """
    path = []
     
    for p in path_coords:       
        path.append(np.argmin(np.linalg.norm(grid-p,axis=-1)))
    
    return np.array(path)

def project_path_periodic(path_coords, grid,M,m,axes, dim=3, D_=None):
    """! @brief Project path coordinates onto grid points with periodic boundaries.

    For each point in path_coords, finds the closest grid point accounting
    for periodic boundary conditions, and returns indices and projection errors.

    @param path_coords  Array of coordinates to project (shape: (n, dim)).
    @param grid         Array of grid point coordinates (shape: (npts, dim)).
    @param M            Upper bounds of the simulation box (shape: (dim,)).
    @param m            Lower bounds of the simulation box (shape: (dim,)).
    @param axes         Array of axes along which periodic boundaries apply.
    @param dim          Dimensionality of the dataset (default: 3).
    @param D_           Precomputed distance matrix (default: None).
    @return Tuple (path, err) where path is an array of grid indices and
            err is an array of projection errors.
    """   

    if D_ is None:
        D = dist_from_pts_periodic_boundaries(path_coords,grid,M,m,axes, dim=3)
    else:
        D = D_
        
    path = []
    err = []

    for i,p in enumerate(path_coords):
        proj = np.argmin(D[i,:])
        e = np.min(D[i,:])
        path.append(proj)
        err.append(e)
        
    return np.array(path), np.array(err)


def project_timeseries_array(A,grid,M,m,axes, dim=3):
    """! @brief Project a time-series of coordinate arrays onto grid points.

    Projects each frame in A onto the grid with periodic boundaries and
    accumulates unique grid indices across all frames.

    @param A     Time-series array of coordinates (shape: (T, n, dim)).
    @param grid  Array of grid point coordinates (shape: (npts, dim)).
    @param M     Upper bounds of the simulation box (shape: (dim,)).
    @param m     Lower bounds of the simulation box (shape: (dim,)).
    @param axes  Array of axes along which periodic boundaries apply.
    @param dim   Dimensionality of the dataset (default: 3).
    @return Tuple (p, err) where p is the array of unique grid indices and
            err is the maximum projection error across all frames.
    """
    a = A[0]
    p,e = project_path_periodic(a, grid,M,m,axes, dim)
    
    err = np.max(e)
    
    for a in A[1:]:   
        
        p_,e = project_path_periodic(a, grid,M,m,axes, dim)
        p = np.unique(np.hstack([p,p_]))
    
        err = np.max([err, np.max(e)])
    
    return p,err


def cluster_deadlocks(PATHS,res_grid, geod, res_to_A_0,geod_A_0,A_0):
    """! @brief Cluster deadlocked paths by homological equivalence.

    Groups paths into clusters where paths within each cluster are
    homologically equivalent, using close_deadlocks and path_to_vec
    for equivalence testing.

    @param PATHS        List of paths to cluster.
    @param res_grid     Array of grid point coordinates (shape: (npts, dim)).
    @param geod         Predecessor matrix for geodesic path reconstruction.
    @param res_to_A_0   Mapping from resolution grid indices to A_0 indices.
    @param geod_A_0     Predecessor matrix for geodesics on A_0.
    @param A_0          Bottom boundary grid point coordinates.
    @return List of clusters, each a list of homologically equivalent paths.
    """
    DEADLOCKS=[[PATHS[0]]]
    
    for p in PATHS[1:]:
        
        new=1
        
        for c in DEADLOCKS:
            aux = 1
            for q in c:
                path = close_deadlocks(p,q, res_grid, geod, res_to_A_0,geod_A_0,A_0)
                
                if len(path)==0:
                    aux=0
                    break
            
            if aux==1:
                v_path = path_to_vec(path,D,N_1,edge_idxs)
                res = equivalence(A,v_path,rank)
                if res == 0:
                    c.append(p_s)
                    new = 1
                    break
            
        if new == 1:    
            DEADLOCKS.append([p])

          
        
    return DEADLOCKS

def save_paths(LIST,name, folder_path):
    """! @brief Save a list of paths to a numpy file.

    Pads all paths to equal length with -1 and saves the resulting
    matrix to a .npy file in the specified folder.

    @param LIST         List of paths (arrays of grid point indices).
    @param name         Output file name (without directory prefix).
    @param folder_path  Directory path where the file will be saved.
    """
    n=0
        
    for p in LIST:
        if len(p)>n:
            n=len(p)
            
    paths = np.zeros((len(LIST),n))-1
    
    for i,p in enumerate(LIST):
        paths[i,:len(p)] = p
        
    
    np.save(folder_path+name,paths)

"""
Parallelization Utils
"""

def compute_all(LIST):
    """! @brief Compute directional and penalty integrals along geodesic paths.

    For each pair (i, j) with j < i, finds the geodesic path and integrates
    the direction field v and the gradient penalty field along it.

    @param LIST  Packed argument list: (i, v, grad_pen, multi_idx_to_list, D, geod, res_grid).
    @return List of [directional_integral, penalty_integral] pairs for each j < i.
    """
    i,v,grad_pen,multi_idx_to_list,D,geod,res_grid = LIST
    
    out = []
    
    for j in range(i):
        path, lenghts = find_geod(i,j,geod,res_grid)
        lenghts = lenghts/D[i,j]

        line = np.diff(res_grid[path], axis=0)   

        out0 = np.sum(np.sum(-v*(line),axis=1)*lenghts)
        out1 = np.sum(np.sum(-grad_pen[path[:-1]]*(line),axis=1)*lenghts)
    
        out.append([out0,out1])
    
    return out

def compute_penalties(LIST):
    """! @brief Compute directional and penalty integrals for a single pair using precomputed paths.

    Retrieves the precomputed geodesic path for pair (i, j) and integrates
    the direction and gradient penalty fields along it.

    @param LIST  Packed argument list: (i, j, v, grad_pen, RESULTS, multi_idx_to_list, D, res_grid).
    @return List [directional_integral, penalty_integral].
    """
    i,j,v,grad_pen,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    out0 = np.sum(np.sum(-v*(line),axis=1)*lenghts)
    out1 = np.sum(np.sum(-grad_pen[path[:-1]]*(line),axis=1)*lenghts)
    
    return [out0,out1]

def compute_direction(LIST):
    """! @brief Compute the directional integral for a single pair using precomputed paths.

    Retrieves the precomputed geodesic path for pair (i, j) and integrates
    the direction field v along it.

    @param LIST  Packed argument list: (i, j, v, RESULTS, multi_idx_to_list, D, res_grid).
    @return Scalar directional integral value.
    """
    i,j,v,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    return np.sum(np.sum(-v*(line),axis=1)*lenghts)

def compute_repulsion(LIST):
    """! @brief Compute the repulsion penalty integral for a single pair using precomputed paths.

    Retrieves the precomputed geodesic path for pair (i, j) and integrates
    the gradient penalty field along it.

    @param LIST  Packed argument list: (i, j, grad_pen, RESULTS, multi_idx_to_list, D, res_grid).
    @return Scalar repulsion penalty integral value.
    """
    i,j,grad_pen,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    return np.sum(np.sum(-grad_pen[path[:-1]]*(line),axis=1)*lenghts)

def find_geod_wrap(LIST):
    """! @brief Wrapper for find_geod suitable for multiprocessing pool.map.

    Unpacks arguments from a list and calls find_geod to compute the
    geodesic path between two points.

    @param LIST  Packed argument list: (i, j, geod, res_grid).
    @return List [path, lengths] from find_geod.
    """
    i,j,geod,res_grid = LIST
    path, lenghts = find_geod(i,j,geod,res_grid)
    
    return [path,lenghts]

def rank_mod_p(LIST):
    """! @brief Compute the rank of a matrix modulo a prime p.

    Row-reduces the matrix over the Galois field GF(p) and returns
    the matrix rank.

    @param LIST  Packed argument list: (m, p) where m is the matrix and
                 p is the prime (defaults to 2 if None).
    @return Integer rank of the matrix over GF(p).
    """
    m, p = LIST
    
    if p is None:
        p=2
    
    GF = galois.GF(p)
    m_2 = GF(m.astype(int)%p)
    m_red = m_2.row_reduce()
    
    return np.linalg.matrix_rank(m_red)


"""
Wrappers
"""