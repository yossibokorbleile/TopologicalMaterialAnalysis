import oineus as oin
from process import *

atoms, radii, repeat_x, repeat_y, repeat_z = read_configuration("structure-types.ini", "ZIF-TEST")
xyz = load_atom_file("/home/ubuntu/ZIF-Data/ZIF-76/3/ZIF_300_wrapped.xyz")

points = sample_at(xyz, 4005, repeat_x, repeat_y, repeat_z, atoms, radii)
sub = sub_complex(points, 23.3, 0.03)
K,L,L_to_K = oineus_pair(points, sub)
params = oin.ReductionParams()
params.n_threads=32
params.verbose=True
params.kernel=True
kicr = oin.compute_kernel_image_cokernel_reduction(K, L, L_to_K, params)
