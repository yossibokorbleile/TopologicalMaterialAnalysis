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
    
    I = np.identity(A.shape[0])
    I = GF(I.astype(int)%2)

    AI = np.concatenate((A, I), axis=-1)
    Q = AI.row_reduce()
    V = Q[:,-A.shape[0]:]
    A_red = Q[:,:A.shape[1]] 

    return A_red, V 


def solve_linear_system_with_boundary(B, v, A, A_red, V, GF):
    """
    Finds a solution - if it exists - of B*x = v + A*mu
    
    B is row-reduced form
    A generates the boundary subspace
    A_red is A[rank(B):,:] in row reduced form via the matrix V
    tmp tells if the solution is exact
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
    """
    Finds a solution - if it exists - of B*x = v
    
    B is row-reduced form
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

    start_idx = 0
    path_ = np.copy(path)
    
    while grid[path_[start_idx+1]][axis]<=h*1.00001:
        path_ = path_[1:]

    return np.array(path_,dtype=int)

    
def cut_path_from_top(path, h, grid, axis=-1):

    start_idx = -1
    path_ = np.copy(path)
    
    while grid[path_[start_idx-1]][axis]>=h*0.999999:
        path_ = path_[:-1]

    return np.array(path_,dtype=int)


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def get_box_from_cell(cell_file):
    
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
    """
    Read the .xyz file and returns atoms coordinates by type, cardinalities and box boundaries
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
    """
    Remove Lithium Atoms Outside the Box
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
    """
    Makes a 
    
    grid -> array of pts (npts,dim)
    M -> array with upper bounds of the box (dim,)
    m -> array with lower bound of the box (dim,)
    axes -> array with the axes in which the boundary condition is applied
    dim -> dimension of the dataset
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
    """
    pts -> punti da proiettare
    grid -> array of pts (npts,dim)
    M -> array with upper bounds of the box (dim,)
    m -> array with lower bound of the box (dim,)
    axes -> array with the axes in which the boundary condition is applied
    dim -> dimension of the dataset
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
    """
    N.B. ia<ib<ic
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
    
    i0 = cubic_to_res[p0]
    i1 = cubic_to_res[p1]
    i2 = cubic_to_res[p2]
    
    i0, i1, i2 = np.sort([i0,i1,i2]) 

    if i0 > -1:        
        ITRIS, B = append_tr(ITRIS, B, edge_idxs, i0,i1,i2) 
    
    
    return ITRIS, B 

def r_squared(px,py,pz,balls_centres,balls_radii):
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
      
    dx = 0
    dy = 0
    dz = 0

    for i,c in enumerate(balls_centres):

        dx += aux_fn(pt[0],pt[1],pt[2],c[0],c[1],c[2],balls_radii[i],eps,p)
        dy += aux_fn(pt[1],pt[0],pt[2],c[1],c[0],c[2],balls_radii[i],eps,p)
        dz += aux_fn(pt[2],pt[1],pt[0],c[2],c[1],c[0],balls_radii[i],eps,p)
        
    return -np.array([dx,dy,dz])

def sample_path(idx, N, res_grid, z_1, maxiter=100000):
    
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
    
    idx = 0
    path_smooth = []

    while idx < len(path)-1:
        tmp,_= find_geod(path[idx],path[idx+1],geod,res_grid,D)
        path_smooth += list(tmp)
        idx+=1
        
    return np.array(path_smooth)

def close_paths(p_0, p_1, grid, geod,D):
    
    start_point_path,_ = find_geod(p_0[0], p_1[0], geod, grid,D) 
    end_point_path,_ = find_geod(p_0[-1], p_1[-1], geod, grid,D)

    path = list(start_point_path)[::-1] + list (p_0) +\
       list(end_point_path) + list(p_1)[::-1] 
    
    return np.array(path)


def homological_simplification_of_path_pts(idx0, idx1, path, grid, geod, V, rank, GF, D, N, edge_idxs, r_graph):
    """
    idx0,idx1 (p0,p1) -> are indexes to sample from the path
    q0,q1 -> are the indexes sampled from the paths and refer to pts in the grid
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
    
    start_point_path_,_ = find_geod(res_to_A_0[p_0[0]],res_to_A_0[p_1[0]],geod_A_0,A_0,D) 
    end_point_path_,_ = find_geod(res_to_A_1[p_0[-1]],res_to_A_1[p_1[-1]],geod_A_1,A_1,D)

    start_point_path = [A_0_to_res[p] for p in start_point_path_]
    end_point_path = [A_1_to_res[q] for q in end_point_path_]

    path = list(start_point_path)[::-1] + list (p_0) +\
       list(end_point_path) + list(p_1)[::-1] 
    
    return np.array(path)

def join_points_boundary_conditions(p_, q_, M, m, axes, dim, n=20):
    
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
    
    path = []
     
    for p in path_coords:       
        path.append(np.argmin(np.linalg.norm(grid-p,axis=-1)))
    
    return np.array(path)

def project_path_periodic(path_coords, grid,M,m,axes, dim=3, D_=None):
    """
    path_coords -> coordinates of the points of the path
    grid -> array of pts to project the path on (npts,dim)
    M -> array with upper bounds of the box (dim,)
    m -> array with lower bound of the box (dim,)
    axes -> array with the axes in which the boundary condition is applied
    dim -> dimension of the dataset
    D_ -> distance matrix if already computed
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
    
    a = A[0]
    p,e = project_path_periodic(a, grid,M,m,axes, dim)
    
    err = np.max(e)
    
    for a in A[1:]:   
        
        p_,e = project_path_periodic(a, grid,M,m,axes, dim)
        p = np.unique(np.hstack([p,p_]))
    
        err = np.max([err, np.max(e)])
    
    return p,err


def cluster_deadlocks(PATHS,res_grid, geod, res_to_A_0,geod_A_0,A_0):
    
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
    
    i,j,v,grad_pen,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    out0 = np.sum(np.sum(-v*(line),axis=1)*lenghts)
    out1 = np.sum(np.sum(-grad_pen[path[:-1]]*(line),axis=1)*lenghts)
    
    return [out0,out1]

def compute_direction(LIST):
    
    i,j,v,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    return np.sum(np.sum(-v*(line),axis=1)*lenghts)

def compute_repulsion(LIST):
    
    i,j,grad_pen,RESULTS,multi_idx_to_list,D,res_grid = LIST
    
    path, lengths = RESULTS[multi_idx_to_list[(i,j)]]
    lenghts = lengths/D[i,j]
               
    line = np.diff(res_grid[path], axis=0)   

    return np.sum(np.sum(-grad_pen[path[:-1]]*(line),axis=1)*lenghts)

def find_geod_wrap(LIST):
    i,j,geod,res_grid = LIST 
    path, lenghts = find_geod(i,j,geod,res_grid)
    
    return [path,lenghts]

def rank_mod_p(LIST):
    
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