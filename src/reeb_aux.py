##
# @file reeb_aux.py
# @brief Auxiliary functions for Reeb graph construction including Numba-optimised grid operations.
# @version 1.3.0
# @date March 2026

import numpy as np
from numba import jit, int32, prange
from diffusion_utils import dist_from_pts_periodic_boundaries, out_to_atoms
from itertools import chain, combinations, product
import multiprocessing as mp
from scipy.sparse.csgraph import connected_components
import os

# Configure Numba threading to prevent concurrent access issues
# Use default threading layer
# Limit threads to prevent oversubscription
os.environ['NUMBA_NUM_THREADS'] = str(min(4, mp.cpu_count()))


import networkx as nx
import pyomo.environ as pyo
import logging
logging.getLogger('pyomo.core').setLevel(logging.ERROR)

# Thread synchronization for Numba functions
import threading
_numba_lock = threading.Lock()

def thread_safe_numba_call(func, *args, **kwargs):
    """! @brief Wrapper to ensure thread-safe execution of Numba functions.

    Acquires a global lock before calling the Numba-compiled function to
    prevent concurrent access issues.

    @param func    The Numba-compiled function to call.
    @param args    Positional arguments forwarded to func.
    @param kwargs  Keyword arguments forwarded to func.
    @return The return value of func(*args, **kwargs).
    """
    with _numba_lock:
        return func(*args, **kwargs)

# Thread-safe wrappers for commonly used Numba functions
def safe_dist_from_pts_periodic_boundaries_numba(A_, B_, M, m, axes, dim=3):
    """! @brief Thread-safe wrapper for dist_from_pts_periodic_boundaries_numba.

    @param A_    First set of points (shape: (n, dim)).
    @param B_    Second set of points (shape: (m, dim)).
    @param M     Upper bounds of the simulation box (shape: (dim,)).
    @param m     Lower bounds of the simulation box (shape: (dim,)).
    @param axes  Array of axes along which periodic boundaries apply.
    @param dim   Dimensionality of the dataset (default: 3).
    @return Distance matrix of shape (n, m).
    """
    return thread_safe_numba_call(dist_from_pts_periodic_boundaries_numba, A_, B_, M, m, axes, dim)

def safe_point_cloud_frechet_mean_numba(D, M, m, subsample=10, tol=0.001, maxiter=50):
    """! @brief Thread-safe wrapper for point_cloud_frechet_mean_numba.

    @param D          Point cloud data array (shape: (T, L, dim)).
    @param M          Upper bounds of the simulation box (shape: (dim,)).
    @param m          Lower bounds of the simulation box (shape: (dim,)).
    @param subsample  Subsampling stride for initial point selection (default: 10).
    @param tol        Convergence tolerance (default: 0.001).
    @param maxiter    Maximum number of iterations (default: 50).
    @return Frechet mean coordinates (shape: (L, dim)).
    """
    return thread_safe_numba_call(point_cloud_frechet_mean_numba, D, M, m, subsample, tol, maxiter)

def safe_flow_to_fmean_dist(flow, backbone, M, m):
    """! @brief Thread-safe wrapper for flow_to_fmean_dist.

    @param flow      Flow atom coordinates (shape: (T, n_flow, 3)).
    @param backbone  Backbone Frechet mean coordinates (shape: (n_backbone, 3)).
    @param M         Upper bounds of the simulation box (shape: (3,)).
    @param m         Lower bounds of the simulation box (shape: (3,)).
    @return Distance matrix between backbone atoms and flow atoms across time.
    """
    return thread_safe_numba_call(flow_to_fmean_dist, flow, backbone, M, m)


def circular_max_flow(reeb):
    """! @brief Compute the maximum flow through a circular Reeb graph.

    Constructs a flow network from the Reeb graph by adding source and
    sink nodes connected to the bottom and top level-set vertices, then
    solves the max-flow problem using Gurobi via Pyomo.

    @param reeb  Reeb graph object with attributes E (edges), REEB_GRAPH,
                 G (NetworkX graph), and D_reeb_w (edge weights).
    @return Maximum flow value through the network.
    """
    edges = reeb.E+1
    
    pts_bottom = np.unique(reeb.REEB_GRAPH[:,0])
    pts_bottom = pts_bottom[pts_bottom>-1]
    pts_top = np.unique(reeb.REEB_GRAPH[:,-1])
    pts_top = pts_top[pts_top>-1]
    
    add_edges_bottom = np.vstack([np.zeros_like(pts_bottom),pts_bottom + 1]).T
    add_edges_top = np.vstack([pts_top + 1, np.zeros_like(pts_top) + np.max(reeb.G.nodes) + 2]).T
    
    constr_edges = np.zeros((len(add_edges_bottom),4))
#    constr_edges = np.hstack([add_edges_bottom, add_edges_top])
    constr_edges[:,:2] = add_edges_bottom
    constr_edges[:,2:] = add_edges_top
    
    edges = np.vstack([add_edges_bottom,edges,add_edges_top])
    edges_to_idx = {}

    for e in edges:
        edges_to_idx[tuple(e)] = np.argmin(np.linalg.norm(edges-e,axis=-1)) 
        
    N = np.sum(reeb.D_reeb_w)
    cap = {}
    
    for e in reeb.E:
        cap[tuple(e+1)] = reeb.G.edges[tuple(e)]['capacity']
        
    for e in add_edges_bottom:
        cap[tuple(e)] = N
    
    for e in add_edges_top:
        cap[tuple(e)] = N
    
    cost = pyo.ConcreteModel()
    cost.n = len(edges) 

    cost.f = pyo.Var(np.arange(cost.n),
                     domain=pyo.NonNegativeReals,initialize=0)
    
    flow = sum([cost.f[edges_to_idx[(0,v+1)]] for v in pts_bottom]) 

    cost.obj=pyo.Objective(rule=flow, sense=pyo.maximize)
    cost.costr=pyo.ConstraintList()

    for e in edges:
        cost.costr.add(cost.f[edges_to_idx[tuple(e)]] <= cap[tuple(e)])

    for e in constr_edges:
        cost.costr.add(cost.f[edges_to_idx[tuple(e[:2])]] == cost.f[edges_to_idx[tuple(e[2:])]])


    for v in list(reeb.G.nodes):
        
        neigh = list(reeb.G[v])
        
        if v in pts_bottom:
            neigh = neigh+[-1]  
        elif v in pts_top:
            neigh = neigh+[np.max(reeb.G.nodes) + 1]  
            
        if len(neigh)>0:
            flx_in = sum([cost.f[edges_to_idx[(u+1,v+1)]] for u in neigh if u<v])
            flx_out = sum([cost.f[edges_to_idx[(v+1,u+1)]] for u in neigh if u>v])
            cost.costr.add(flx_in-flx_out==0)

    solver = pyo.SolverFactory('gurobi')
    S=solver.solve(cost)

    flx = cost.obj()

    f_eval = np.zeros_like(edges[:,0])

    for key,val in cost.f.extract_values().items():
                f_eval[key]=val 

    flow_eval = np.sum([f_eval[edges_to_idx[(0,v+1)]] for v in pts_bottom]) 
    
    return flow_eval


def sphere_sample(n=100, seed=0):
    """! @brief Sample points uniformly on the unit sphere.

    Generates n random points on the 3D unit sphere by normalising
    standard Gaussian samples.

    @param n     Number of points to sample (default: 100).
    @param seed  Random seed for reproducibility (default: 0).
    @return Array of unit-sphere points (shape: (n, 3)).
    """
    np.random.seed(seed)
    pts = np.random.normal(size=(n,3))    
    norms = np.linalg.norm(pts,axis=-1)
    
    for i in range(3):
        pts[:,i] = pts[:,i]/norms
    return pts


def sample_matrix_grid(n=100,seed=0):
    """! @brief Generate a grid of rotation matrices from sphere samples.

    Constructs n rotation matrices by applying Gram-Schmidt orthogonalisation
    to frames built from sphere-sampled direction vectors.

    @param n     Number of rotation matrices to generate (default: 100).
    @param seed  Random seed for sphere sampling (default: 0).
    @return Array of rotation matrices (shape: (n, 3, 3)).
    """
    pts = sphere_sample(n, seed)

    MAT = []
    b,c = [0,1,0],[0,0,1]
    
    for a in pts:
        A = np.array([a,b,c]).T
        m,_ = gramschmidt(A)
        m = m[:,[2,1,0]]
        if np.linalg.det(m)<0:
            m = m[:,[0,2,1]]
        MAT.append(m)
        
    return np.array(MAT)

def canonical_matrix_grid():
    """! @brief Generate the three canonical axis-permutation rotation matrices.

    Returns the identity matrix with its rows cyclically permuted to produce
    rotations aligning each axis to the first position.

    @return Array of three permutation matrices (shape: (3, 3, 3)).
    """
    MAT = []
    m = np.identity(3)
    
    for p in [(0,1,2),(1,2,0),(2,0,1)]:
        MAT.append(m[np.array(p)])
    
    return np.array(MAT)

def gramschmidt(A):
    """! @brief Apply the Gram-Schmidt orthogonalisation process to a matrix.

    Decomposes matrix A into an orthogonal matrix Q and an upper triangular
    matrix R such that A = Q * R.

    @param A  Input matrix to decompose (shape: (n, k)).
    @return Tuple (Q, R) where Q is orthogonal (shape: (n, k)) and R is
            upper triangular (shape: (k, k)).
    """
    R = np.zeros((A.shape[1], A.shape[1]))
    Q = np.zeros(A.shape)
    
    for k in range(0, A.shape[1]):
        R[k, k] = np.sqrt(np.dot(A[:, k], A[:, k]))
        Q[:, k] = A[:, k]/R[k, k]
        for j in range(k+1, A.shape[1]):
            R[k, j] = np.dot(Q[:, k], A[:, j])
            A[:, j] = A[:, j] - R[k, j]*Q[:, k]
    return Q, R

@jit(nopython=False, fastmath=True)
def make_cubic_to_res(cubic_grid,res_grid,nx,ny,nz):
    """! @brief Build mappings between cubic grid and resolution grid indices.

    For each cell in the cubic grid, finds the matching point in the resolution
    grid and records bidirectional index mappings.

    @param cubic_grid  Cubic grid coordinates (shape: (nx, ny, nz, 3)).
    @param res_grid    Resolution grid coordinates (shape: (npts, 3)).
    @param nx          Number of grid cells along x.
    @param ny          Number of grid cells along y.
    @param nz          Number of grid cells along z.
    @return Tuple (cubic_to_res, res_to_cubic) where cubic_to_res maps
            (i,j,k) to resolution index (-1 if absent) and res_to_cubic
            maps resolution index to (i,j,k).
    """
    cubic_to_res = np.zeros((nx,ny,nz), dtype=int32)-1
    res_to_cubic = np.zeros((len(res_grid),3), dtype=int32)
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                p = cubic_grid[i,j,k,:]
                tmp = np.argwhere(np.sum(np.abs(res_grid-p),axis=1)<0.00001)

                if len(tmp)>0:                    
                    cubic_to_res[i,j,k] = tmp[0][0]
                    res_to_cubic[tmp[0][0],:] = [i,j,k]
                    
    return cubic_to_res, res_to_cubic


@jit(nopython=False, fastmath=True)
def make_graph(grid,res_grid,nx,ny,nz,dim,M,m):
    """! @brief Build a neighbourhood graph on the resolution grid using a cubic lattice.

    Reshapes the flat grid into a cubic lattice, maps to resolution indices,
    and creates edges between each resolution point and its 26-connected
    cubic neighbours (with periodic wrapping in x and y).

    @param grid      Flat grid coordinates (shape: (nx*ny*nz, 3)).
    @param res_grid  Resolution grid coordinates (shape: (npts, 3)).
    @param nx        Number of grid cells along x.
    @param ny        Number of grid cells along y.
    @param nz        Number of grid cells along z.
    @param dim       Dimensionality of the dataset.
    @param M         Upper bounds of the simulation box (shape: (dim,)).
    @param m         Lower bounds of the simulation box (shape: (dim,)).
    @return Tuple (graph_try, cnt) where graph_try is an edge array
            (shape: (npts*27, 2)) and cnt is the number of valid edges.
    """
    cubic_grid = np.transpose(grid.reshape((nx,ny,nz,3)),axes=(1,0,2,3))
    cubic_to_res, res_to_cubic = make_cubic_to_res(cubic_grid,res_grid,nx,ny,nz)

    graph_try = np.zeros((len(res_grid)*27,2),dtype=int32)
    cnt = 0
        
    for p_res in range(len(res_grid)): 
        
        i,j,k = res_to_cubic[p_res,:]
        p = cubic_grid[i,j,k,:]

        for v in [(1,0,0),(0,1,0),(0,0,1),
                   (-1,0,0),(0,-1,0),(0,0,-1),
                   (1,1,0),(0,1,1),(1,0,1),
                   (-1,-1,0),(0,-1,-1),(-1,0,-1),
                   (1,-1,0),(0,1,-1),(1,0,-1),
                   (-1,1,0),(0,-1,1),(-1,0,1),
                   (1,1,-1),(1,-1,1),(-1,1,1),
                   (-1,-1,1),(-1,1,-1),(1,-1,-1),
                   (1,1,1),(-1,-1,-1),(0,0,0)]:

            q_cub = ((i+v[0])%nx,(j+v[1])%ny,(k+v[2])%(nz+1))
            
            if q_cub[-1]<nz and q_cub[-1]>-1:
               
                q_res = cubic_to_res[q_cub]
                if q_res>-1:
                    graph_try[cnt,:]=[p_res,q_res]
                    cnt+=1
    
    return graph_try, cnt

@jit(nopython=False, parallel=False, fastmath=True)
def make_graph_fast(grid,res_idxs,
                   res_to_grid,grid_to_res,
                   grid_to_cubic, cubic_to_grid,
                   nx,ny,nz,dim,M,m):
    """! @brief Build a neighbourhood graph using precomputed grid-to-cubic mappings.

    A faster variant of make_graph that uses precomputed index mappings
    instead of recomputing them, iterating only over resolution grid points.

    @param grid            Full grid coordinates (shape: (N, 3)).
    @param res_idxs        Boolean/indicator array marking resolution grid membership.
    @param res_to_grid     Mapping from resolution indices to full grid indices.
    @param grid_to_res     Mapping from full grid indices to resolution indices.
    @param grid_to_cubic   Mapping from full grid indices to (i, j, k) tuples.
    @param cubic_to_grid   Mapping from (i, j, k) to full grid indices (shape: (nx, ny, nz)).
    @param nx              Number of grid cells along x.
    @param ny              Number of grid cells along y.
    @param nz              Number of grid cells along z.
    @param dim             Dimensionality of the dataset.
    @param M               Upper bounds of the simulation box (shape: (dim,)).
    @param m               Lower bounds of the simulation box (shape: (dim,)).
    @return Tuple (graph_try, cnt) where graph_try is an edge array and
            cnt is the number of valid edges.
    """
    graph_try = np.zeros((len(res_to_grid)*27,2),dtype=int32)
    cnt = 0
        
    for i_res,idx in enumerate(res_to_grid):
        
        p = grid[idx]
        i,j,k = grid_to_cubic[idx]


        for v in [(1,0,0),(0,1,0),(0,0,1),
                   (-1,0,0),(0,-1,0),(0,0,-1),
                   (1,1,0),(0,1,1),(1,0,1),
                   (-1,-1,0),(0,-1,-1),(-1,0,-1),
                   (1,-1,0),(0,1,-1),(1,0,-1),
                   (-1,1,0),(0,-1,1),(-1,0,1),
                   (1,1,-1),(1,-1,1),(-1,1,1),
                   (-1,-1,1),(-1,1,-1),(1,-1,-1),
                   (1,1,1),(-1,-1,-1),(0,0,0)]:

            q_cub = ((i+v[0])%nx,(j+v[1])%ny,(k+v[2])%(nz+1))

            if q_cub[-1]<nz and q_cub[-1]>-1:

                q_grid = cubic_to_grid[q_cub[0],q_cub[1],q_cub[2]]

                if res_idxs[q_grid]>0:
                    q_res = grid_to_res[q_grid]
                    graph_try[cnt,:]=[i_res,q_res]
                    cnt+=1
                            
    return graph_try, cnt


@jit(nopython=False, fastmath=True)
def match_grids(new_grid, old_grid):
    """! @brief Find index correspondences between two grids.

    For each point in old_grid, finds its matching point in new_grid
    (within tolerance 1e-5) and records the index mapping.

    @param new_grid  Target grid coordinates (shape: (n, dim)).
    @param old_grid  Source grid coordinates (shape: (m, dim)).
    @return Index mapping array (shape: (m,)) where entry i gives the
            new_grid index for old_grid[i], or -1 if no match found.
    """
    idxs = np.zeros_like(old_grid, dtype=int32)-1
    
    for i,p in enumerate(old_grid):
        tmp = np.argwhere(np.sum(np.abs(new_grid-p),axis=1)<0.00001)

        if len(tmp)>0:                    
            idxs[i] = tmp[0][0]
    return idxs

@jit(nopython=False, fastmath=True)
def match_graphs(old_to_new, old_graph):
    """! @brief Remap a graph's vertex indices using an old-to-new index mapping.

    Translates edges from old_graph into the new index space, discarding
    edges where either endpoint has no valid mapping.

    @param old_to_new  Index mapping from old to new vertex indices (shape: (n,)).
    @param old_graph   Edge list in old indices (shape: (E, 2)).
    @return Tuple (graph, cnt) where graph is the remapped edge array and
            cnt is the number of valid edges.
    """
    graph = np.zeros((len(old_to_new),2),dtype=int32)
    cnt = 0

    for (p,q) in old_graph:
        if old_to_new[p]>-1 and old_to_new[q]>-1:
            graph[cnt,:] = [p,q]
            
    return graph, cnt


@jit(nopython=False, parallel=True)
def dist_from_pts_periodic_boundaries_pool(LIST):
    """! @brief Pool-compatible wrapper for dist_from_pts_periodic_boundaries_numba.

    Unpacks arguments from a list and delegates to the Numba-optimised
    periodic boundary distance function.

    @param LIST  Packed argument list: (A, B, M, m, axes, dim).
    @return Distance matrix from dist_from_pts_periodic_boundaries_numba.
    """
    A,B,M,m,axes, dim = LIST

    d = dist_from_pts_periodic_boundaries_numba(A,B,M,m,axes, dim)
    
    return d


@jit(nopython=True, parallel=True, fastmath=True)
def numba_dist(a, b):
    """! @brief Compute pairwise Euclidean distances using Numba parallelisation.

    Calculates the distance d(a_i, b_j) for every pair of points
    a_i in a and b_j in b using parallel loops.

    @param a  First array of points (shape: (n, dim)).
    @param b  Second array of points (shape: (m, dim)).
    @return Distance matrix of shape (n, m).
    """
    dist = np.zeros((a.shape[0],b.shape[0]))
    for a_i in prange(a.shape[0]):
        for b_j in prange(b.shape[0]):
            dist_ai_bj = 0.0
            for i in range(a.shape[1]):
                diff = b[b_j, i] - a[a_i, i]
                dist_ai_bj = dist_ai_bj + diff**2
            dist[a_i,b_j] = dist_ai_bj**0.5
    return dist

@jit(nopython=True, parallel=False, fastmath=True)
def dist_from_pts_periodic_boundaries_numba(A_,B_,M,m,axes, dim=3):
    """! @brief Compute pairwise distances with periodic boundaries using Numba.

    Numba-optimised version of dist_from_pts_periodic_boundaries. Wraps
    coordinates into the simulation box and evaluates distances across
    all periodic images.

    @param A_    First set of points (shape: (n, dim)).
    @param B_    Second set of points (shape: (m, dim)).
    @param M     Upper bounds of the simulation box (shape: (dim,)).
    @param m     Lower bounds of the simulation box (shape: (dim,)).
    @param axes  Array of axes along which periodic boundaries apply.
    @param dim   Dimensionality of the dataset (default: 3).
    @return Distance matrix of shape (n, m).
    """    

    A = m + np.remainder(A_-m, M-m)
    B = m + np.remainder(B_-m, M-m)

    D = numba_dist(A, B)
    DELTAS = M-m
    
    c_idxs = np.arange(dim)
    c_idxs = np.array([i not in axes for i in c_idxs]) 
    DELTAS[c_idxs] = 0.0
    
    for v_0 in [0.0,-DELTAS[0],DELTAS[0]]:
        for v_1 in [0.0,-DELTAS[1],DELTAS[1]]:
            for v_2 in [0.0,-DELTAS[2],DELTAS[2]]:

                v = np.array([v_0,v_1,v_2])
                
                if np.min(v)==np.max(v)==0:
                    pass
                else:
                    B_new = B + v
                    D_new = numba_dist(A, B_new)
                    D = np.minimum(D_new,D)

    return D


# @jit(nopython=False, parallel=True, fastmath=True)
def numba_dist_print(a, b):
    """! @brief Compute pairwise Euclidean distances with debug printing.

    Same as numba_dist but with print statements for debugging. Not
    JIT-compiled to allow print output.

    @param a  First array of points (shape: (n, dim)).
    @param b  Second array of points (shape: (m, dim)).
    @return Distance matrix of shape (n, m).
    """
    print("computing distances between ", a, " and ", b)
    # Ensure contiguous arrays and float64 dtype for scalar math
    # a_c = np.ascontiguousarray(a)
    # b_c = np.ascontiguousarray(b)
    dist = np.zeros((a.shape[0], b.shape[0]), dtype=np.float64)
    assert a.shape[1] == b.shape[1]
    for a_i in range(a.shape[0]):
        for b_j in range(b.shape[0]):
            dist_ai_bj = 0.0
            for i in range(a.shape[1]):
                diff = b[b_j, i] - a[a_i, i]
                dist_ai_bj = dist_ai_bj + diff**2
            dist[a_i, b_j] = dist_ai_bj ** 0.5
            # dist[a_i,b_j] = dist[a_i,b_j]**0.5
    return dist


@jit(nopython=False, parallel=True, fastmath=True)
def dist_from_pts_periodic_boundaries_numba_print(A_,B_,M,m,axes, dim=3):
    """! @brief Compute periodic boundary distances with debug printing.

    Debug variant of dist_from_pts_periodic_boundaries_numba that prints
    array shapes for diagnostic purposes.

    @param A_    First set of points (shape: (n, dim)).
    @param B_    Second set of points (shape: (m, dim)).
    @param M     Upper bounds of the simulation box (shape: (dim,)).
    @param m     Lower bounds of the simulation box (shape: (dim,)).
    @param axes  Array of axes along which periodic boundaries apply.
    @param dim   Dimensionality of the dataset (default: 3).
    @return Distance matrix of shape (n, m).
    """    
    
    print("A_ has shape: ", A_[0].shape)
    print("B_ has shape: ", B_.shape)
    A = m + np.remainder(A_-m, M-m)
    B = m + np.remainder(B_-m, M-m)
    A = np.ascontiguousarray(A)
    B = np.ascontiguousarray(B)
    D = numba_dist_print(A, B)
    DELTAS = M-m
    
    c_idxs = np.arange(dim)
    c_idxs = np.array([i not in axes for i in c_idxs]) 
    DELTAS[c_idxs] = 0.0
    for v_0 in [0.0,-DELTAS[0],DELTAS[0]]:
        for v_1 in [0.0,-DELTAS[1],DELTAS[1]]:
            for v_2 in [0.0,-DELTAS[2],DELTAS[2]]:

                v = np.array([v_0,v_1,v_2])
                
                if np.min(v)==np.max(v)==0:
                    pass
                else:
                    B_new = B + v
                    D_new = numba_dist_print(A, B_new)
                    D = np.minimum(D_new,D)

    return D

@jit(nopython=False, parallel=True, fastmath=True)
def point_cloud_frechet_mean_numba(D,M,m,subsample=10,tol=0.001, maxiter = 50):
    """! @brief Compute the Frechet mean of a point cloud with periodic boundaries (Numba).

    For each landmark point, computes the Frechet mean across time steps
    using the Procrustes mean algorithm with periodic boundary wrapping.

    @param D          Point cloud data (shape: (T, L, dim)).
    @param M          Upper bounds of the simulation box (shape: (dim,)).
    @param m          Lower bounds of the simulation box (shape: (dim,)).
    @param subsample  Subsampling stride for initial points (default: 10).
    @param tol        Convergence tolerance (default: 0.001).
    @param maxiter    Maximum number of iterations per mean computation (default: 50).
    @return Frechet mean coordinates (shape: (L, dim)).
    """
    frechet = np.zeros_like(D[0,:,:])
    N = D.shape[0]
    L = D.shape[1]
    
    for i in range(L):
        
        pt = D[:,i,:]
    #    idxs = np.random.choice(np.arange(0,N),size=(subsample,), replace=False) 
        idxs = np.arange(0,N,subsample)
        initial_pts = pt[idxs,:]
        frechet[i,:] = frechet_mean_numba(pt, initial_pts, M, m, tol = tol, maxiter = maxiter)
    
    return frechet

@jit(nopython=False, parallel=False, fastmath=True)
def procustes_mean(p, D, M, m, tol = 0.001, maxiter = 50):
    """! @brief Compute the Procrustes (Frechet) mean from a single initial point.

    Iteratively refines the mean estimate by centring the data around
    the current estimate, averaging, and wrapping back with periodic
    boundary conditions until convergence.

    @param p        Initial point estimate (shape: (dim,)).
    @param D        Data points to average (shape: (N, dim)).
    @param M        Upper bounds of the simulation box (shape: (dim,)).
    @param m        Lower bounds of the simulation box (shape: (dim,)).
    @param tol      Convergence tolerance on cost change (default: 0.001).
    @param maxiter  Maximum number of iterations (default: 50).
    @return Tuple (mean_coords, cost) where mean_coords is the converged
            mean (shape: (dim,)) and cost is the final mean squared distance.
    """
    cnt = 0
    err = 100
    cost = 0
    N = D.shape[0]
    mean_coords = p
    axes = np.array([0,1,2])
    center = (M+m)/2

    while err > tol and cnt < maxiter:
        mean_tmp = mean_coords
        v = center - mean_tmp
        cost_tmp = cost

        new_coords = m + np.remainder(D+v-m, M-m)    
        mean = np.sum(new_coords,axis=0)/N           
        mean_coords = m + np.remainder(mean-v-m, M-m)
        a = np.ascontiguousarray(mean_coords).reshape(1,-1)
        b = np.ascontiguousarray(mean_tmp).reshape(1,-1)

        cost = dist_from_pts_periodic_boundaries_numba(D,a,M,m,axes)    
        cost = np.mean(np.power(cost,2))

        err = np.abs(cost-cost_tmp)

        cnt+=1

    return mean_coords, cost


@jit(nopython=False, parallel=False, fastmath=True)
def frechet_mean_numba(D,initial_points,M,m,tol=0.001, maxiter = 50):
    """! @brief Compute the Frechet mean by selecting the best among multiple initial points (Numba).

    Runs procustes_mean from each initial point and returns the result
    with the lowest cost.

    @param D               Data points to average (shape: (N, dim)).
    @param initial_points  Candidate initial points (shape: (L, dim)).
    @param M               Upper bounds of the simulation box (shape: (dim,)).
    @param m               Lower bounds of the simulation box (shape: (dim,)).
    @param tol             Convergence tolerance (default: 0.001).
    @param maxiter         Maximum iterations per initial point (default: 50).
    @return Best Frechet mean coordinates (shape: (dim,)).
    """
    L = initial_points.shape[0]
    
    MEANS = np.zeros_like(initial_points)
    COSTS = np.zeros_like(initial_points[:,0])
    
    for i in range(L):
        
        mean, cost = procustes_mean(initial_points[i,:], D, M, m, tol=tol, maxiter = maxiter)
    
        MEANS[i,:] = mean
        COSTS[i] = cost
    
    idx = np.argmin(COSTS)
    
    return MEANS[idx,:]

def point_cloud_frechet_mean(D,M,m,subsample=10,tol=0.0001, maxiter = 50):
    """! @brief Compute the Frechet mean of a point cloud with periodic boundaries.

    Non-Numba version of point_cloud_frechet_mean_numba. Uses random
    subsampling for initial point selection.

    @param D          Point cloud data (shape: (T, L, dim)).
    @param M          Upper bounds of the simulation box (shape: (dim,)).
    @param m          Lower bounds of the simulation box (shape: (dim,)).
    @param subsample  Number of random initial points to try (default: 10).
    @param tol        Convergence tolerance (default: 0.0001).
    @param maxiter    Maximum iterations per mean computation (default: 50).
    @return Frechet mean coordinates (shape: (L, dim)).
    """
    frechet = np.zeros_like(D[0,:,:])
    N = D.shape[0]
    L = D.shape[1]
    
    for i in range(L):
        
        pt = D[:,i,:]
        idxs = np.random.choice(np.arange(0,N),size=(subsample,), replace=False) 
        initial_pts = pt[idxs,:]
        
        frechet[i,:] = frechet_mean(pt,initial_pts,M,m,tol=0.0001, maxiter = 50)
    
    return frechet


def frechet_mean(D,initial_points,M,m,tol=0.0001, maxiter = 50):
    """! @brief Compute the Frechet mean by selecting the best among multiple initial points.

    Non-Numba version of frechet_mean_numba. Iterates the Procrustes mean
    algorithm from each initial point and returns the lowest-cost result.

    @param D               Data points to average (shape: (N, dim)).
    @param initial_points  Candidate initial points (shape: (L, dim)).
    @param M               Upper bounds of the simulation box (shape: (dim,)).
    @param m               Lower bounds of the simulation box (shape: (dim,)).
    @param tol             Convergence tolerance (default: 0.0001).
    @param maxiter         Maximum iterations per initial point (default: 50).
    @return Best Frechet mean coordinates (shape: (dim,)).
    """
    axes = np.array([0,1,2])
    center = (M+m)/2
    N = D.shape[0]
    
    MEANS = np.zeros_like(initial_points)
    COSTS = np.zeros_like(initial_points[:,0])*1000000
    

    for i,p in enumerate(initial_points):
        
        mean_coords = p
        
        cnt = 0
        err = 100
        cost = 0

        while err > tol and cnt < maxiter:
            mean_tmp = mean_coords
            v = center - mean_tmp
            cost_tmp = cost
            
            new_coords = m + np.remainder(D+v-m, M-m)    
            mean = np.mean(new_coords,axis=0)
            mean_coords = m + np.remainder(mean-v-m, M-m)

            cost = dist_from_pts_periodic_boundaries_numba(D,mean_coords.reshape(1,-1),M,m,axes)    
            cost = np.mean(np.power(cost,2))

            err = np.abs(cost-cost_tmp)            
            cnt+=1
        
        MEANS[i,:] = mean_coords
        COSTS[i] = cost
    
    idx = np.argmin(COSTS)    
    return MEANS[idx,:]


def  estimate_radius(inputfile, backbone_atoms, flow_atoms, n_tmp_grid):
    """! @brief Estimate exclusion radii for backbone atoms from simulation data.

    Reads atom trajectories, computes Frechet means for backbone and flow
    atoms, and derives per-atom exclusion radii based on minimum distances
    to flow atoms.

    @param inputfile       Path to the input XYZ trajectory file.
    @param backbone_atoms  List of element symbols for backbone atoms.
    @param flow_atoms      List of element symbols for flow atoms.
    @param n_tmp_grid      Temporal subsampling stride for the trajectory.
    @return Tuple (fmean, M_Li, radii) where fmean is the backbone Frechet
            mean, M_Li is the preprocessed flow coordinates, and radii is
            the array of per-atom exclusion radii.
    """
    backbone, flow, N, timesteps, M, m = out_to_atoms(inputfile)
    # Li,P,S,N,timesteps,M,m = out_to_atoms(inputfile)

    tmp = np.arange(0,len(Li),n_tmp_grid)
    
    Li = Li[tmp,:,:]
    P = P[tmp,:,:]
    S = S[tmp,:,:]

    Li = m + np.remainder(Li-m, M-m)    
    P = m + np.remainder(P-m, M-m)    
    S = m + np.remainder(S-m, M-m)                
    BACKBONE = np.hstack([P,S])
   
    n_Li = Li.shape[1]
    n_P = P.shape[1]
    n_S = S.shape[1]

    M_P = safe_point_cloud_frechet_mean_numba(P, M, m, subsample=min([20,len(tmp)]), tol=0.001, maxiter = 20)            
    M_S = safe_point_cloud_frechet_mean_numba(S, M, m, subsample=min([20,len(tmp)]), tol=0.001, maxiter = 20)

    fmean = np.concatenate([M_P, M_S])
    M_Li = preprocess_PATHS(fmean, BACKBONE, Li, M, m)
            
#    D = Li_to_backbone_dist(Li,P,S,M,m)
    D = Li_to_backbone_dist(M_Li,M_P,M_S,M,m)

    r_P = np.min(D[:n_P,:],axis=-1)
    r_S = np.min(D[n_P:,:],axis=-1)
    
    return fmean, M_Li, np.concatenate([r_P,r_S])


@jit(nopython=False, parallel=False, fastmath=True)
def Li_to_backbone_dist(Li,P,S,M,m):
    """! @brief Compute distances from backbone atoms to lithium ions across time steps.

    For each time step, concatenates P and S backbone coordinates and
    computes distances to all Li positions with periodic boundaries.

    @param Li  Lithium ion coordinates (shape: (T, n_Li, 3)).
    @param P   Phosphorus atom coordinates (shape: (T, n_P, 3)).
    @param S   Sulphur atom coordinates (shape: (T, n_S, 3)).
    @param M   Upper bounds of the simulation box (shape: (3,)).
    @param m   Lower bounds of the simulation box (shape: (3,)).
    @return Distance matrix (shape: (n_P + n_S, T * n_Li)).
    """
    tmp = len(Li)
    axes_aux = np.array([0,1,2])
    RES = np.zeros((len(P[0,:,0])+len(S[0,:,0]),tmp*len(Li[0,:,0])))
        
    for t in range(tmp):
        L = Li[t,:,:] 
        P_ = P[t,:,:]
        S_ = S[t,:,:] 

        B = np.concatenate((P_,S_))
        RES[:,t*len(Li[0,:,0]):(t+1)*len(Li[0,:,0])] = dist_from_pts_periodic_boundaries_numba(B,L,M,m,axes_aux)
    
    return RES

@jit(nopython=False, parallel=False, fastmath=True)
def flow_to_fmean_dist(flow,backbone,M,m):
    """! @brief Compute distances from backbone Frechet mean to flow atoms across time steps.

    For each time step, computes distances from the backbone mean
    positions to all flow atom positions with periodic boundaries.

    @param flow      Flow atom coordinates (shape: (T, n_flow, 3)).
    @param backbone  Backbone Frechet mean coordinates (shape: (n_backbone, 3)).
    @param M         Upper bounds of the simulation box (shape: (3,)).
    @param m         Lower bounds of the simulation box (shape: (3,)).
    @return Distance matrix (shape: (n_backbone, T * n_flow)).
    """
    tmp = len(flow)
    axes_aux = np.array([0,1,2])
    RES = np.zeros((len(backbone[:,0]),tmp*len(flow[0,:,0])))
        
    for t in range(tmp):
        L = flow[t,:,:] 
        RES[:,t*len(flow[0,:,0]):(t+1)*len(flow[0,:,0])] = dist_from_pts_periodic_boundaries_numba(backbone,L,M,m,axes_aux)
    
    return RES


@jit(nopython=False, parallel=False, fastmath=True)
def preprocess_ions(fmean, backbone_t, ions, M, m):
    """! @brief Shift ion coordinates towards the Frechet mean frame.

    Translates each ion position using a distance-weighted average of
    backbone-to-Frechet-mean displacement vectors, with periodic wrapping.

    @param fmean       Backbone Frechet mean coordinates (shape: (n_backbone, 3)).
    @param backbone_t  Backbone coordinates at the current time step (shape: (n_backbone, 3)).
    @param ions        Ion coordinates to preprocess (shape: (n_ions, 3)).
    @param M           Upper bounds of the simulation box (shape: (3,)).
    @param m           Lower bounds of the simulation box (shape: (3,)).
    @return Preprocessed ion coordinates (shape: (n_ions, 3)).
    """
    axes = np.array([0,1,2])
    matrix = dist_from_pts_periodic_boundaries_numba(ions,backbone_t,M,m,axes)
    new_ions = np.zeros_like(ions)
    center = (M+m)/2
    
    for i in range(len(ions)):
        p = ions[i,:]
        w = np.exp(-np.power(matrix[i,:],2))
        w = w/np.sum(w)
                
        v = center-backbone_t
        centr_fmean = m + np.remainder(fmean+v-m, M-m)
        v_tang = centr_fmean-center
    
        v_mean = np.sum(v_tang.T*w,axis=1).T
                
        new_coords = p+v_mean
        tmp = m + np.remainder(new_coords-m, M-m)
        new_ions[i,:] = tmp
                
    return new_ions 

@jit(nopython=False, parallel=False, fastmath=True)
def preprocess_PATHS(fmean, BACKBONE, PATHS, M, m):
    """! @brief Preprocess ion trajectories by shifting towards the Frechet mean frame.

    Applies preprocess_ions to each time step of the trajectory to align
    ion paths with the Frechet mean backbone configuration.

    @param fmean     Backbone Frechet mean coordinates (shape: (n_backbone, 3)).
    @param BACKBONE  Backbone coordinates across time (shape: (T, n_backbone, 3)).
    @param PATHS     Ion trajectories to preprocess (shape: (T, n_ions, 3)).
    @param M         Upper bounds of the simulation box (shape: (3,)).
    @param m         Lower bounds of the simulation box (shape: (3,)).
    @return Preprocessed trajectories (shape: (T, n_ions, 3)).
    """
    new_PATHS = np.zeros_like(PATHS)
    time_steps = PATHS.shape[0]
    
    for t in range(time_steps):
        backbone_t = BACKBONE[t,:,:]
        ions = PATHS[t,:,:]
        new_PATHS[t,:,:]=preprocess_ions(fmean, backbone_t, ions, M, m)
        
    return new_PATHS


@jit(nopython=False, parallel=False, fastmath=True)
def make_neigh(grid,p_graph,i_graph,f):
    """! @brief Compute the maximum neighbourhood function variation on the grid.

    For each grid point, finds the maximum absolute difference between its
    function value and those of its neighbours (including itself), then
    returns the global maximum across all points.

    @param grid     Grid point coordinates (shape: (npts, dim)).
    @param p_graph  CSR row pointer array for the adjacency graph.
    @param i_graph  CSR column index array for the adjacency graph.
    @param f        Function values at grid points (shape: (npts,)).
    @return Global maximum neighbourhood variation (scalar).
    """
    aux = np.zeros_like(grid[:,0])
    
    for i in range(len(grid)):

        neigh = np.zeros((p_graph[i+1]-p_graph[i]+1),dtype=int32)
        neigh[-1]=i
        
        for j in range(p_graph[i],p_graph[i+1]):
            neigh[j-p_graph[i]]=i_graph[j]
            
        aux[i] = np.max(np.abs(f[neigh]-f[i]))

    return np.max(aux)

@jit(nopython=False, parallel=False, fastmath=True)
def make_neigh_(grid,p_graph,i_graph,f):
    """! @brief Compute the maximum neighbourhood function variation (scalar variant).

    Same as make_neigh but tracks the maximum as a running scalar instead
    of allocating an auxiliary array, for improved memory efficiency.

    @param grid     Grid point coordinates (shape: (npts, dim)).
    @param p_graph  CSR row pointer array for the adjacency graph.
    @param i_graph  CSR column index array for the adjacency graph.
    @param f        Function values at grid points (shape: (npts,)).
    @return Global maximum neighbourhood variation (scalar).
    """
    aux = 0
    
    for i in range(len(grid)):
        neigh = np.zeros((p_graph[i+1]-p_graph[i]+1),dtype=int32)
        neigh[-1]=i
        
        for j in range(p_graph[i],p_graph[i+1]):
            neigh[j-p_graph[i]]=i_graph[j]
            
        a = np.max(np.abs(f[neigh]-f[i]))
        
        if a>aux:
            aux=a
        
    return aux


@jit(nopython=False, parallel=False, fastmath=True)
def connect_lvl_sets_aux(i, n_c, C, grid, lvl_to_grid, REEB_GRAPH, norm, count):
    """! @brief Connect components of adjacent level sets for Reeb graph construction.

    Assigns Reeb graph node labels to grid points at level i based on their
    connected component membership C, then computes edge weights between
    level i and level i-1 nodes by counting shared grid points.

    @param i              Current level-set index.
    @param n_c            Number of connected components at this level.
    @param C              Component labels for points at this level (shape: (n_level,)).
    @param grid           Grid point coordinates (shape: (npts, 3)).
    @param lvl_to_grid    Mapping from level-set indices to grid indices.
    @param REEB_GRAPH     Reeb graph label matrix (shape: (npts, n_levels)), modified in place.
    @param norm           Normalisation factor for edge weights and node sizes.
    @param count          Running counter for Reeb graph node identifiers.
    @return Tuple (REEB_GRAPH_col, AUX_w, AUX_f, AUX_e, count) where
            REEB_GRAPH_col is the updated column i of REEB_GRAPH, AUX_w
            contains weighted edges, AUX_f contains node sizes, AUX_e
            contains node mean coordinates, and count is the updated counter.
    """
    AUX_w = np.zeros((n_c*len(np.unique(REEB_GRAPH[:,i-1])),3))
    AUX_f = np.zeros((n_c,))
    AUX_e = np.zeros((n_c,3)) 

    cnt = 0
    
    for c in range(n_c):
        tmp_ = np.where(C==c)[0]
        tmp = lvl_to_grid[tmp_]  

        REEB_GRAPH[tmp,i] = count 
        idxs_aux = np.unique(REEB_GRAPH[tmp,i-1])        
        idxs_aux = idxs_aux[idxs_aux>-1]
        AUX_f[c] = np.max(np.array([np.rint(len(tmp)/norm),1]))
        AUX_e[c,:] =  np.array([np.mean(grid[tmp][:,ax]) for ax in range(3)])

        for j in prange(len(idxs_aux)):
            old_c = idxs_aux[j]
            now = (REEB_GRAPH[:,i] == count)
            old = (REEB_GRAPH[:,i-1] == old_c)
            inters = np.sum(np.multiply(now,old))

            if inters>0:
                AUX_w[cnt+j,:] = np.array([old_c,count,np.max(np.array([inters//norm,1]))])

            
        cnt+= len(idxs_aux)        
        count += 1

    aux = AUX_w[:,-1]>0

    return REEB_GRAPH[:,i], AUX_w[aux,:], AUX_f, AUX_e, count
