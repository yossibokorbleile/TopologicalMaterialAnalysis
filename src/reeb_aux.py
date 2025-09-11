import numpy as np
from numba import jit, int32, prange
from diffusion_utils import dist_from_pts_periodic_boundaries, out_to_atoms
from itertools import chain, combinations, product
import multiprocessing as mp
from scipy.sparse.csgraph import connected_components


import networkx as nx
import pyomo.environ as pyo
import logging
logging.getLogger('pyomo.core').setLevel(logging.ERROR)


def circular_max_flow(reeb):
    
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
    
    np.random.seed(seed)
    pts = np.random.normal(size=(n,3))    
    norms = np.linalg.norm(pts,axis=-1)
    
    for i in range(3):
        pts[:,i] = pts[:,i]/norms
    return pts


def sample_matrix_grid(n=100,seed=0):
    
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
    
    MAT = []
    m = np.identity(3)
    
    for p in [(0,1,2),(1,2,0),(2,0,1)]:
        MAT.append(m[np.array(p)])
    
    return np.array(MAT)


def gramschmidt(A):
    """
    Applies the Gram-Schmidt method to A
    and returns Q and R, so Q*R = A.
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

@jit(nopython=True, fastmath=True)
def make_cubic_to_res(cubic_grid,res_grid,nx,ny,nz):

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


@jit(nopython=True, fastmath=True)    
def make_graph(grid,res_grid,nx,ny,nz,dim,M,m):

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

@jit(nopython=True, parallel=False, fastmath=True)    
def make_graph_fast(grid,res_idxs,
                   res_to_grid,grid_to_res,
                   grid_to_cubic, cubic_to_grid,
                   nx,ny,nz,dim,M,m):

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


@jit(nopython=True, fastmath=True)
def match_grids(new_grid, old_grid):

    idxs = np.zeros_like(old_grid, dtype=int32)-1
    
    for i,p in enumerate(old_grid):
        tmp = np.argwhere(np.sum(np.abs(new_grid-p),axis=1)<0.00001)

        if len(tmp)>0:                    
            idxs[i] = tmp[0][0]
    return idxs

@jit(nopython=True, fastmath=True)
def match_graphs(old_to_new, old_graph):

    graph = np.zeros((len(old_to_new),2),dtype=int32)
    cnt = 0

    for (p,q) in old_graph:
        if old_to_new[p]>-1 and old_to_new[q]>-1:
            graph[cnt,:] = [p,q]
            
    return graph, cnt


#@jit(nopython=True, parallel=True)
def dist_from_pts_periodic_boundaries_pool(LIST):
    
    A,B,M,m,axes, dim = LIST

    d = dist_from_pts_periodic_boundaries_numba(A,B,M,m,axes, dim)
    
    return d


@jit(nopython=True, parallel=True, fastmath=True)
def numba_dist(a, b):
    """
    A,B -> array of pts (npts,dim), calcolo distanza d(a_i,b_j) per ogni a_i in a e b_j in b
    """
    dist = np.zeros((a.shape[0],b.shape[0]))
    for a_i in prange(a.shape[0]):
        for b_j in prange(b.shape[0]):
            for i in range(a.shape[1]):
                dist[a_i,b_j] += (b[b_j,i] - a[a_i, i])**2
            dist[a_i,b_j] = dist[a_i,b_j]**0.5
    return dist

@jit(nopython=True, parallel=True, fastmath=True)
def dist_from_pts_periodic_boundaries_numba(A_,B_,M,m,axes, dim=3):
    """
    A,B -> array of pts (npts,dim), calcolo distanza d(a,b) per ogni a in A e b in B
    M -> array with upper bounds of the box (dim,)
    m -> array with lower bound of the box (dim,)
    axes -> array with the axes in which the boundary condition is applied
    dim -> dimension of the dataset
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

@jit(nopython=True, parallel=True, fastmath=True)
def point_cloud_frechet_mean_numba(D,M,m,subsample=10,tol=0.001, maxiter = 50):

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

@jit(nopython=True, parallel=False, fastmath=True)
def procustes_mean(p, D, M, m, tol = 0.001, maxiter = 50):
    
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


@jit(nopython=True, parallel=False, fastmath=True)
def frechet_mean_numba(D,initial_points,M,m,tol=0.001, maxiter = 50):
    
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


def estimate_radius(inputfile, backbone_atoms, flow_atoms, n_tmp_grid):
    
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

    M_P = point_cloud_frechet_mean_numba(P, M, m, subsample=min([20,len(tmp)]), tol=0.001, maxiter = 20)            
    M_S = point_cloud_frechet_mean_numba(S, M, m, subsample=min([20,len(tmp)]), tol=0.001, maxiter = 20)

    fmean = np.concatenate([M_P, M_S])
    M_Li = preprocess_PATHS(fmean, BACKBONE, Li, M, m)
            
#    D = Li_to_backbone_dist(Li,P,S,M,m)
    D = Li_to_fmean_dist(M_Li,M_P,M_S,M,m)

    r_P = np.min(D[:n_P,:],axis=-1)
    r_S = np.min(D[n_P:,:],axis=-1)
    
    return fmean, M_Li, np.concatenate([r_P,r_S])


@jit(nopython=True, parallel=True, fastmath=True)
def Li_to_backbone_dist(Li,P,S,M,m):
    
    tmp = len(Li)
    axes_aux = np.array([0,1,2])
    RES = np.zeros((len(P[0,:,0])+len(S[0,:,0]),tmp*len(Li[0,:,0])))
        
    for t in prange(tmp):
        L = Li[t,:,:] 
        P_ = P[t,:,:]
        S_ = S[t,:,:] 

        B = np.concatenate((P_,S_))
        RES[:,t*len(Li[0,:,0]):(t+1)*len(Li[0,:,0])] = dist_from_pts_periodic_boundaries_numba(B,L,M,m,axes_aux)
    
    return RES

@jit(nopython=True, parallel=True, fastmath=True)
def Li_to_fmean_dist(Li,P,S,M,m):
    
    tmp = len(Li)
    axes_aux = np.array([0,1,2])
    RES = np.zeros((len(P[:,0])+len(S[:,0]),tmp*len(Li[0,:,0])))
        
    for t in prange(tmp):
        L = Li[t,:,:] 
        B = np.concatenate((P,S))
        RES[:,t*len(Li[0,:,0]):(t+1)*len(Li[0,:,0])] = dist_from_pts_periodic_boundaries_numba(B,L,M,m,axes_aux)
    
    return RES


@jit(nopython=True, parallel=False, fastmath=True)            
def preprocess_ions(fmean, backbone_t, ions, M, m):
    
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

@jit(nopython=True, parallel=False, fastmath=True)            
def preprocess_PATHS(fmean, BACKBONE, PATHS, M, m):
    
    new_PATHS = np.zeros_like(PATHS)
    time_steps = PATHS.shape[0]
    
    for t in range(time_steps):
        backbone_t = BACKBONE[t,:,:]
        ions = PATHS[t,:,:]
        new_PATHS[t,:,:]=preprocess_ions(fmean, backbone_t, ions, M, m)
        
    return new_PATHS


@jit(nopython=True, parallel=True, fastmath=True)
def make_neigh(grid,p_graph,i_graph,f):

    aux = np.zeros_like(grid[:,0])
    
    for i in prange(len(grid)):

        neigh = np.zeros((p_graph[i+1]-p_graph[i]+1),dtype=int32)
        neigh[-1]=i
        
        for j in range(p_graph[i],p_graph[i+1]):
            neigh[j-p_graph[i]]=i_graph[j]
            
        aux[i] = np.max(np.abs(f[neigh]-f[i]))

    return np.max(aux)

@jit(nopython=True, parallel=False, fastmath=True)
def make_neigh_(grid,p_graph,i_graph,f):
    aux = 0
    
    for i in prange(len(grid)):
        neigh = np.zeros((p_graph[i+1]-p_graph[i]+1),dtype=int32)
        neigh[-1]=i
        
        for j in range(p_graph[i],p_graph[i+1]):
            neigh[j-p_graph[i]]=i_graph[j]
            
        a = np.max(np.abs(f[neigh]-f[i]))
        
        if a>aux:
            aux=a
        
    return aux


@jit(nopython=True, parallel=False, fastmath=True)    
def connect_lvl_sets_aux(i, n_c, C, grid, lvl_to_grid, REEB_GRAPH, norm, count):
    
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
