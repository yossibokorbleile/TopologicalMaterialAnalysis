##
# @file example_multi_timestep.py
# @brief Example: Analysing a material across multiple time steps.
#
# This example demonstrates how to iterate over multiple time steps
# from a molecular dynamics trajectory, computing persistent homology
# at each step and comparing the results via the APF.
#
# @section ex3_run Running this example
# @code{.sh}
# cd TopologicalMaterialAnalysis
# python examples/example_multi_timestep.py
# @endcode

import sys
sys.path.insert(0, "src")

import oineus
from toma_functions import (
    read_configuration, sample_at,
    oineus_process, calculate_APF, write_dgm_csv
)

def main():
    """! Main function demonstrating multi-timestep analysis."""

    # Load configuration
    print("=== Loading configuration ===")
    atoms, radii, _, _, _, repeat_x, repeat_y, repeat_z = \
        read_configuration("examples/settings.ini", "ZIF")

    params = oineus.ReductionParams()
    params.n_threads = 4

    # Define the time step range
    sample_start = 0
    sample_end = 2
    sample_step = 1

    all_dgms_1 = []
    all_apfs_1 = []

    print(f"\n=== Processing time steps {sample_start} to {sample_end-1} ===")

    for t in range(sample_start, sample_end, sample_step):
        print(f"\n--- Time step {t} ---")

        # Load the structure at this time step
        points = sample_at(
            file_path="examples/ZIF_test.xyz",
            format="xyz",
            sample_index=t,
            repeat_x=repeat_x,
            repeat_y=repeat_y,
            repeat_z=repeat_z,
            atom_list=atoms,
            radius_list=radii
        )
        print(f"  Loaded {len(points)} atoms")

        # Compute persistent homology
        dcmp, filt, dgm_0, dgm_1, dgm_2 = oineus_process(points, params)

        # Calculate APF
        apf_1 = calculate_APF(dgm_1)

        all_dgms_1.append(dgm_1)
        all_apfs_1.append(apf_1)

        print(f"  Dim 0: {len(dgm_0)} features")
        print(f"  Dim 1: {len(dgm_1)} features")
        print(f"  Dim 2: {len(dgm_2)} features")
        print(f"  APF dim 1 total: {apf_1['lifetime'].iloc[-1]:.4f}")

    # Summary comparison
    print(f"\n=== Summary across {len(all_dgms_1)} time steps ===")
    for i, (dgm, apf) in enumerate(zip(all_dgms_1, all_apfs_1)):
        print(f"  Step {sample_start + i * sample_step}: "
              f"{len(dgm)} dim-1 features, "
              f"APF total = {apf['lifetime'].iloc[-1]:.4f}")

    print("\n=== Done! ===")

if __name__ == "__main__":
    main()
