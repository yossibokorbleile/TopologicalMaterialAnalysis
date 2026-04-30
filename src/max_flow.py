import pickle

import numpy as np

from reeb_aux import *
from reeb_graph import Reeb_Graph


def compute_backbone(
    ase_atoms,
    backbone_atoms,
    flow_atoms,
    reeb_stride,
    common_backbone=True,
    apply_frechet_mean=False,
    print_mean_stucture=False,
    repeat=(1, 1, 1),
):
    """Compute backbone radii and flow coordinates"""
    print("processing backbone")
    n_steps = len(ase_atoms)

    if repeat != (1, 1, 1):
        unwrapped_atoms = unwrap_trajectory(
            ase_atoms, np.diagonal(ase_atoms[0].get_cell())
        )
        repeated_frames = []
        for atoms in unwrapped_atoms:
            new_atoms = atoms.repeat(repeat)
            new_atoms.wrap()
            repeated_frames.append(new_atoms)
    else:
        repeated_frames = ase_atoms

    cell = np.diagonal(repeated_frames[0].get_cell())
    atoms = backbone_atoms + flow_atoms
    n_backbone_types = len(backbone_atoms)

    atom_locations = []
    for i in range(0, n_steps, reeb_stride):
        atom_locations.append(repeated_frames[i].get_positions())

    atom_types = repeated_frames[0].get_chemical_symbols()
    atom_coords = np.ascontiguousarray(atom_locations)

    m = np.array([0, 0, 0])
    backbone_coords = []
    flow_coords = []
    backbone_idxs = []
    flow_idxs = []
    atom_types = np.array(atom_types)

    atom_idxs = []
    for i in range(len(atoms)):
        atom_idxs.append(np.where(np.isin(atom_types, [atoms[i]]))[0].tolist())
        if i < n_backbone_types:
            backbone_idxs.extend(atom_idxs[i])
        else:
            flow_idxs.extend(atom_idxs[i])
    backbone_atom_types = atom_types[backbone_idxs]
    backbone_coords = atom_coords[:, backbone_idxs, :]
    flow_coords = atom_coords[:, flow_idxs, :]

    if apply_frechet_mean:
        print("Using Fréchet mean for backbone computation. This may take a while...")
        backbone_mean = safe_point_cloud_frechet_mean_numba(
            backbone_coords,
            cell,
            m,
            subsample=min([20, len(atom_coords)]),
            tol=0.001,
            maxiter=20,
        )

    else:
        unwrapped_atoms = unwrap_trajectory(
            repeated_frames, np.diagonal(repeated_frames[0].get_cell())
        )
        mean_positions = np.mean([at.get_positions() for at in unwrapped_atoms], axis=0)
        mean_positions = np.remainder(mean_positions, cell)
        backbone_mean = mean_positions[backbone_idxs]

    M_flow = preprocess_PATHS(backbone_mean, backbone_coords, flow_coords, cell, m)

    if print_mean_stucture:
        from ase import Atoms

        back_structure = Atoms(
            symbols=backbone_atom_types,
            positions=backbone_mean,
            cell=cell,
            pbc=True,
        )

        flow_structure = Atoms(
            symbols=atom_types[flow_idxs],
            positions=M_flow[0, :, :],
            cell=cell,
            pbc=True,
        )

        mean_structure = back_structure + flow_structure
        mean_structure.write("mean_structure.xyz")

    if common_backbone is True:
        radii = []
        for a in backbone_atoms:
            a_mean = backbone_mean[backbone_atom_types == a]
            dist_matrix = flow_to_fmean_dist(M_flow, a_mean, cell, m)
            r = rdf_onset_radius(dist_matrix.flatten(), np.prod(cell))
            # r = np.min(dist_matrix)
            radii.append(r)
            print(r, np.min(dist_matrix))
        backbone_radii = []
        for a in backbone_atom_types:
            backbone_radii.append(float(radii[backbone_atoms.index(a)]))
    else:
        backbone_radii = np.min(
            flow_to_fmean_dist(M_flow, backbone_mean, cell, m), axis=-1
        )
    print("backbone_radii", radii)
    return backbone_mean, backbone_radii, cell, m


def compute_max_flow(
    backbone_mean,
    backbone_radii,
    M,
    m,
    grid_spacing,
    volume,
    fat,
    reeb_stride,
    stride,
    periodic,
    save_rb=False,
    verbose=False,
):
    """Compute reeb_graph and max flow"""

    print("computing circular max flow")
    grid_size = int((volume ** (1 / 3)) / grid_spacing)

    # calculate flow in each cardinal direction and take the average
    # rotate the system and repeat the calculation to get a more isotropic measure of flow
    flows = []
    for i in range(3):
        reeb = Reeb_Graph(
            inputfile=None,
            backbone=backbone_mean,
            radii=backbone_radii,
            M=M,
            m=m,
            grid_size=int(grid_size),
            periodic=periodic,
            fat_radius=float(fat),
            covering=np.array([-1, 1]),
            reeb_stride=int(reeb_stride),
            transform_points=None,
            swap_res_grid_and_balls=False,  # if this is set to false, you use the complementary of the thickening
            relax_z_axis=[int(grid_size) // 2, -(int(grid_size) // 2 + 1)],
            verbose=verbose,
            save_RAM=True,
            stride=int(stride),
            MP=False,
        )
        reeb.make_reeb_graph(plot=False)
        if save_rb:
            with open(f"{save_rb}_dir_{i}.pkl", "wb") as f:
                pickle.dump(reeb, f)
        flow = circular_max_flow(reeb) * (reeb.unit_2d)
        flows.append(flow)
        print(f"Flow in direction {i}: {flow}")

        # Rotate the system for the next iteration
        backbone_mean = np.roll(backbone_mean, shift=1, axis=-1)
        M = np.roll(M, shift=1, axis=0)
        m = np.roll(m, shift=1, axis=0)
    return np.mean(flows)


import matplotlib.pyplot as plt


def rdf_onset_radius(distances_flat, box_volume, n_bins=200, rrange=6, threshold=0.01):
    """
    Find the onset of the first peak in the RDF
    """
    edges = np.linspace(0, rrange, n_bins + 1)
    xval = (edges[1:] + edges[:-1]) / 2
    volbin = (4 / 3) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    h, _ = np.histogram(distances_flat, bins=n_bins, range=(0, rrange))
    h[0] = 0
    gr = (h / volbin) / (distances_flat.size / box_volume)

    peak = np.argmax(gr)
    peak_height = gr[peak]
    for i in range(peak, -1, -1):
        if gr[i] < threshold * peak_height:
            onset_idx = i
            break

    plt.plot(xval, gr)
    plt.axvline(xval[onset_idx], color="red", linestyle="--")
    plt.xlabel("Distance")
    plt.ylabel("g(r)")
    plt.show()

    return xval[onset_idx]


def unwrap_trajectory(atoms, cell):
    unwrapped_atoms_list = [at.copy() for at in atoms]
    crossings = np.zeros((len(atoms[0]), 3))
    previous_positions = unwrapped_atoms_list[0].get_positions()
    for at in unwrapped_atoms_list:
        current_positions = at.get_positions()
        difference = current_positions - previous_positions
        crossings = crossings + np.floor_divide(difference + 0.5 * cell, cell).astype(
            int
        )
        previous_positions = current_positions
        new_positions = current_positions - cell * crossings
        at.set_positions(new_positions)
    return unwrapped_atoms_list
