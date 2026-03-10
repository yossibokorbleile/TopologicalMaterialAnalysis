##
# @file example_basic_persistence.py
# @brief Example: Computing basic persistent homology of a ZIF structure.
#
# This example demonstrates the core ToMA workflow:
# 1. Loading a structure configuration (atom types and radii)
# 2. Sampling an atomic structure from a trajectory file
# 3. Computing persistent homology using weighted alpha complexes
# 4. Calculating the Accumulated Persistence Function (APF)
#
# @section ex1_run Running this example
# @code{.sh}
# cd TopologicalMaterialAnalysis
# python examples/example_basic_persistence.py
# @endcode

import sys
sys.path.insert(0, "src")

import oineus
from toma_functions import read_configuration, sample_at, oineus_process, calculate_APF

def main():
    """! Main function demonstrating basic persistent homology computation."""

    # Step 1: Read the configuration (atom types, radii, repetitions)
    print("=== Step 1: Reading configuration ===")
    atoms, radii, sample_start, sample_end, sample_step, repeat_x, repeat_y, repeat_z = \
        read_configuration("examples/settings.ini", "ZIF")

    print(f"Atom types: {atoms}")
    print(f"Radii: {radii}")
    print(f"Cell repetitions: ({repeat_x}, {repeat_y}, {repeat_z})")

    # Step 2: Sample the structure at time step 0
    print("\n=== Step 2: Loading structure ===")
    points = sample_at(
        file_path="examples/ZIF_test.xyz",
        format="xyz",
        sample_index=0,
        repeat_x=repeat_x,
        repeat_y=repeat_y,
        repeat_z=repeat_z,
        atom_list=atoms,
        radius_list=radii
    )
    print(f"Loaded {len(points)} atoms")
    print(f"Atom types found: {list(points['Atom'].unique())}")
    print(f"Coordinate ranges:")
    print(f"  x: [{points['x'].min():.2f}, {points['x'].max():.2f}]")
    print(f"  y: [{points['y'].min():.2f}, {points['y'].max():.2f}]")
    print(f"  z: [{points['z'].min():.2f}, {points['z'].max():.2f}]")

    # Step 3: Compute persistent homology
    print("\n=== Step 3: Computing persistent homology ===")
    params = oineus.ReductionParams()
    params.n_threads = 4

    dcmp, filt, dgm_0, dgm_1, dgm_2 = oineus_process(points, params)

    print(f"Dimension 0: {len(dgm_0)} features (connected components)")
    print(f"Dimension 1: {len(dgm_1)} features (loops/rings)")
    print(f"Dimension 2: {len(dgm_2)} features (voids/cavities)")

    # Step 4: Compute APFs
    print("\n=== Step 4: Computing Accumulated Persistence Functions ===")
    apf_0 = calculate_APF(dgm_0)
    apf_1 = calculate_APF(dgm_1)
    apf_2 = calculate_APF(dgm_2)

    print(f"APF dim 0 - Total accumulated persistence: {apf_0['lifetime'].iloc[-1]:.4f}")
    print(f"APF dim 1 - Total accumulated persistence: {apf_1['lifetime'].iloc[-1]:.4f}")
    print(f"APF dim 2 - Total accumulated persistence: {apf_2['lifetime'].iloc[-1]:.4f}")

    print("\n=== Done! ===")

if __name__ == "__main__":
    main()
