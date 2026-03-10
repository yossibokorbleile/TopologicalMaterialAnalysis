##
# @file example_kernel_image_cokernel.py
# @brief Example: Computing kernel, image, and cokernel persistence.
#
# This example demonstrates how to compute relative homology using
# kernel/image/cokernel persistence. The subcomplex is defined by
# the top and bottom layers of the structure (by z-coordinate),
# allowing analysis of how surface topology relates to bulk topology.
#
# @section ex2_background Background
# Given a subcomplex L (surface layers) included in a complex K (full structure),
# the inclusion map induces maps on homology. The kernel, image, and cokernel
# of these maps capture different aspects of the topological relationship:
# - **Kernel**: Features born in L that die when mapped to K
# - **Image**: Features from L that persist in K
# - **Cokernel**: Features in K not coming from L
#
# @section ex2_run Running this example
# @code{.sh}
# cd TopologicalMaterialAnalysis
# python examples/example_kernel_image_cokernel.py
# @endcode

import sys
sys.path.insert(0, "src")

import oineus
import pandas
import numpy
from toma_functions import (
    read_configuration, sample_at,
    oineus_kernel_image_cokernel, calculate_APF
)

def main():
    """! Main function demonstrating kernel/image/cokernel persistence."""

    # Load configuration
    print("=== Loading configuration ===")
    atoms, radii, _, _, _, repeat_x, repeat_y, repeat_z = \
        read_configuration("examples/settings.ini", "ZIF")

    # Sample the structure
    print("\n=== Loading structure ===")
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

    # Define the subcomplex thresholds (top and bottom 10%)
    top_pt = max(points["z"])
    bot_pt = min(points["z"])
    height = abs(top_pt - bot_pt)
    thickness = 0.1  # 10% from top and bottom
    upper_threshold = top_pt - thickness * height
    lower_threshold = bot_pt + thickness * height

    print(f"\n=== Subcomplex definition ===")
    print(f"Structure z-range: [{bot_pt:.2f}, {top_pt:.2f}] (height: {height:.2f})")
    print(f"Thickness: {thickness*100:.0f}%")
    print(f"Upper threshold: {upper_threshold:.2f} (points above are in subcomplex)")
    print(f"Lower threshold: {lower_threshold:.2f} (points below are in subcomplex)")

    # Set up parameters
    params = oineus.ReductionParams()
    params.n_threads = 4
    params.kernel = True
    params.image = True
    params.cokernel = True

    # Compute kernel/image/cokernel persistence
    print("\n=== Computing kernel/image/cokernel persistence ===")
    kicr, dgm_0, dgm_1, dgm_2 = oineus_kernel_image_cokernel(
        points, params, upper_threshold, lower_threshold
    )

    # Print codomain (full complex) results
    print(f"\nCodomain (full complex) persistence:")
    print(f"  Dim 0: {len(dgm_0)} features")
    print(f"  Dim 1: {len(dgm_1)} features")
    print(f"  Dim 2: {len(dgm_2)} features")

    # Extract and print kernel persistence
    for dim in range(3):
        kernel_dgm = pandas.DataFrame(
            kicr.kernel_diagrams().in_dimension(dim),
            columns=["birth", "death"]
        )
        image_dgm = pandas.DataFrame(
            kicr.image_diagrams().in_dimension(dim),
            columns=["birth", "death"]
        )
        cokernel_dgm = pandas.DataFrame(
            kicr.cokernel_diagrams().in_dimension(dim),
            columns=["birth", "death"]
        )
        print(f"\nDimension {dim}:")
        print(f"  Kernel:   {len(kernel_dgm)} features")
        print(f"  Image:    {len(image_dgm)} features")
        print(f"  Cokernel: {len(cokernel_dgm)} features")

    print("\n=== Done! ===")

if __name__ == "__main__":
    main()
