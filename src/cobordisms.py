import oineus
from functools import cmp_to_key
import math
from collections import Counter

def generate_permutation(n_verts, sub_verts):
	"""! Generate a permutation of the vertices so that the subvertices are first
	@param n_verts	number of vertices
	@param sub_verts	subvertices

	@return permute_verts	permutation of the vertices
	@return inverse_perm	inverse permutation of the vertices
	"""
	permute_verts = [-1 for i in range(n_verts)]
	inverse_perm = [-1 for i in range(n_verts)]
	n_sub_verts = len(sub_verts)
	n_L = 0
	n_K = 0
	for i in range(n_verts):
		# print(i)
		if i in sub_verts:
			permute_verts[i] = n_L
			inverse_perm[n_L] = i
			n_L += 1
			# if i < len(sub_verts):
			# 	permute_verts[i] = i
			# 	inverse_perm[i] = i
			# else:
			# 	permute_verts[i] = i-len(sub_verts)
			# 	inverse_perm[i-len(sub_verts)] = i
				# permute_verts.append(i-n_verts+len(sub_verts))
		else:
			permute_verts[i] = n_sub_verts + n_K
			inverse_perm[n_sub_verts + n_K] = i
			n_K += 1
			# if i >= len(sub_verts):
			# 	permute_verts[i] = i
			# 	inverse_perm[i] = i
			# 	# permute_verts.append(i)
			# else:
			# 	permute_verts[i] = i+len(sub_verts)
			# 	inverse_perm[i+len(sub_verts)] = i
			# 	# permute_verts.append(i-len(sub_verts)+n_verts)
	return permute_verts, inverse_perm

def generate_permutation_AB(n_verts, sub_verts_A, sub_verts_B):
	"""! Generate a permutation of the vertices so that AuB is first
	@param n_verts	number of vertices
	@param sub_verts_A	subvertices of A
	@param sub_verts_B	subvertices of B

	@return permute_verts	permutation of the vertices
	"""
	permute_verts = []
	for i in range(n_verts):
		# print(i)
		if i in sub_verts_A:
			if i < len(sub_verts_A):
				permute_verts.append(i)
			else:
				permute_verts.append(i-n_verts+len(sub_verts_A))
		elif i in sub_verts_B:
			if i < len(sub_verts_A):
				permute.append(i+len(sub_verts_A))
			elif i < len(sub_verts_B)+len(sub_verts_A):
				permute_verts.append(i-n_verts+len(sub_verts_A)+len(sub_verts_B))
			else:
				permute_verts.append(i-len(sub_verts_A)-len(sub_verts_B)+n_verts)
		else:
			if i < len(sub_verts_A)+len(sub_verts_B):
				permute_verts.append(i-len(sub_verts_A)-len(sub_verts_B)+n_verts)
			else:
				permute_verts.append(i)
	return permute_verts


def dim_compare(x, y):
	"""! Comparison to compare list of simplicies to get them in the order for oineus
	@param x	simplex to compare
	@param y	simplex to compare

	@return -1 if x<y, 1 otherwise
	"""
	if len(x[0]) == len(y[0]):
		if x[0][0] <= y[0][0]:
			return -1
		else:
			return 1
	elif len(x[0]) < len(y[0]):
		return -1
	else:
		return 1

def translate_column(col, perm):
	"""! Translate a column of indices to the new order
	@param col	column of indices
	@param perm	permutation

	@return col	column of indices in the new order
	"""
	for i in range(len(col)):
		col[i] = perm[col[i]]
	col.sort()
	return col

def construct_pair(simplices, n_verts, sub_verts, permute):
	"""! Construct a pair of simplices
	@param simplices	simplices
	@param n_verts	number of vertices
	@param sub_verts	subvertices
	@param permute	permutation

	@return simplices_all	simplices in the new order
	@return simplices_sub	subsimplices in the new order
	"""
	simplices_not_sub = []
	simplices_sub = []
	# permute = generate_permutation(n_verts, sub_verts)
	simplices = sorted(simplices, key=cmp_to_key(dim_compare))
	sub_verts_set = set(sub_verts)
	for s in simplices:
		permuted_verts = [permute[i] for i in s[0]]
		if set(s[0]).issubset(sub_verts_set):
			simplices_sub.append([permuted_verts, s[1]])
		else:
			simplices_not_sub.append([permuted_verts, s[1]])
	set_sub = {tuple(s[0]) for s in simplices_sub}
	set_not_sub = {tuple(s[0]) for s in simplices_not_sub}
	assert set_sub.intersection(set_not_sub) == set()
	simplices_all = simplices_sub + simplices_not_sub
	simplices_sub = [[i,s[0],s[1]] for i, s in enumerate(simplices_sub)]
	simplices_all = [[i,s[0],s[1]] for i, s in enumerate(simplices_all)]
	return simplices_all, simplices_sub

def convert_vertices(vertices, permute):
	"""! Convert the vertices to the new order
	@param vertices	vertices
	@param permute	permutation

	@return new	vertices in the new order
	"""
	new = [permute[v] for v in vertices]
	new.sort()
	return new

def rref_d_phi(d_phi):
	"""! Convert the d_phi to a matrix
	@param d_phi	d_phi

	@return matrix	matrix of the d_phi
	"""
	columns = len(d_phi)
	rows = max([max(col) for col in d_phi])+1
	# print("rows: ", rows, " columns: ", columns)
	matrix = [[0 for j in range(columns)] for i in range(rows)]
	# print("matrix length: ", len(matrix))
	for i, col in enumerate(d_phi):
		for j in col:
			# print("i: ", i, " j: ", j)
			matrix[j][i] = 1
	matrix = sympy.Matrix(matrix)
	# print("matrix: ")
	# print(matrix)
	return matrix


def find_births(dgm_A, dgm_B, dgm_AB):
	"""! Find the points in Cok(d_phi)
	@param dgm_A	time and index diagram of A
	@param dgm_B	time and index diagram of B
	@param dgm_AB	time and index diagram of AB
	@return birth_ids	indices of the births in the order they appear in the cokernel
	"""

	A_birth_ids = dgm_A['birth idx'].tolist()
	B_birth_ids = dgm_B['birth idx'].tolist()
	AB_birth_ids = dgm_AB['birth idx'].tolist()
	all_births = set(A_birth_ids) | set(B_birth_ids) | set(AB_birth_ids)
	# print(all_births)
	# print(A_birth_ids)
	# print(B_birth_ids)
	# print(AB_birth_ids)
	birth_ids = []
	birth_times = []
	for idx in all_births:
		# print("idx is ", idx)
		if idx in A_birth_ids and idx in B_birth_ids and idx in AB_birth_ids:
			# death_ids.append(idx)
			continue
		elif idx not in A_birth_ids and idx not in B_birth_ids and idx in AB_birth_ids:
			birth_ids.append(idx)
			# birth_times.append(dgm_AB['birth'][idx])
		elif idx in A_birth_ids and idx not in B_birth_ids and idx not in AB_birth_ids:
			# birth_ids.append(idx)
			# birth_times.append(dgm_A['birth'][idx])
			continue
		elif idx in A_birth_ids and idx in B_birth_ids and idx not in AB_birth_ids:
			# death_ids.append(idx)
			continue
		else:
			continue
	return birth_ids#, birth_times

def match_points_index(sorted_ids, cok_phi, birth_ids):
	"""! Match the points in Cok(d_phi)
	@param sorted_ids	indices of the simplices in the order they appear in the cokernel
	@param cok_phi		cokernel of d_phi
	@param birth_ids	indices of the births in the order they appear in the cokernel
	@param birth_times	times of the births in the order they appear in the cokernel
	"""	
	# print("birth ids are: ", birth_ids)
	duplicates = [item for item, count in Counter(sorted_ids).items() if count > 1]
	# print("Values that appear twice in sorted_ids:", duplicates)
	print(duplicates)

	# Get the indices for each duplicate value
	duplicate_indices = {}
	for i, value in enumerate(sorted_ids):
		if value in duplicates:
			if value in duplicate_indices:
				duplicate_indices[value].append(i)
			else:
				duplicate_indices[value] = [i]

	# print("Indices for duplicate values:", duplicate_indices)
	
	matched_births = {}
	for i in duplicates:
		print("i is ", i,  " and ", duplicate_indices[i])
		print(cok_phi[duplicate_indices[i][1]])
		if cok_phi[duplicate_indices[i][1]] != [] and cok_phi[duplicate_indices[i][1]][-1] in birth_ids:
			matched_births[cok_phi[duplicate_indices[i][1]][-1]] = i
	# print("matched births are: ", matched_births)
	inf_births = []
	for i in birth_ids:
		if i not in matched_births.keys():
			inf_births.append(i)
	# print("inf births are: ", inf_births)
	index_dgm = []
	for i in matched_births.keys():
		index_dgm.append([i, matched_births[i]])
	for i in inf_births:
		index_dgm.append([i, -1])

	return index_dgm

def convert_index_diagram(index_dgm, kicr):
	"""! Convert the index diagram to the original diagram
	@param index_dgm	index diagram
	@param kicr		KICR object
	"""
	
	persistence_diagram = []
	for pt in index_dgm:
		print("looking at ", pt)
		print(kicr.fil_K.simplex(int(pt[0])).value)
		if pt[1] != -1:
			persistence_diagram.append([kicr.fil_K.simplex(int(pt[0])).value, kicr.fil_K.simplex(int(pt[1])).value, pt[0], pt[1], len(kicr.fil_K.simplex(int(pt[0])).vertices)])
		else:
			persistence_diagram.append([kicr.fil_K.simplex(int(pt[0])).value, math.inf, pt[0], -1, len(kicr.fil_K.simplex(int(pt[0])).vertices)])
	return persistence_diagram
	
def d_phi_generator(kicr_AB, kicr_A, kicr_B, sorted_ids_AB, sorted_ids_A, sorted_ids_B, A_verts, B_verts, perm_AB, perm_A, perm_B, inv_AB, inv_A, inv_B):
	"""! Generate the d_phi matrix
	@param kicr_AB	KICR object of AB
	@param kicr_A	KICR object of A
	@param kicr_B	KICR object of B
	@param sorted_ids_AB	indices of the simplices in the order they appear in the cokernel
	@param sorted_ids_A	indices of the simplices in the order they appear in the cokernel
	@param sorted_ids_B	indices of the simplices in the order they appear in the cokernel
	@param A_verts	vertices of A
	@param B_verts	vertices of B
	@param perm_AB	permutation of the vertices
	@param perm_A	permutation of the vertices
	@param perm_B	permutation of the vertices
	@param inv_AB	inverse permutation of the vertices
	@param inv_A	inverse permutation of the vertices
	@param inv_B	inverse permutation of the vertices

	@return d_phi	matrix of the d_phi
	@return sorted_ids	indices of the simplices in the order they appear in the d_phi
	"""
	d_phi = []
	sorted_ids = []
	simplices = kicr_AB.fil_K.simplices()
	for i in range(len(simplices)):
		simplex = simplices[i]
		og_verts = sorted(convert_vertices(simplex.vertices, inv_AB))
		sorted_id_A = sorted_ids_A[tuple(og_verts)]
		sorted_id_B = sorted_ids_B[tuple(og_verts)]
		print(i, " looking at simplex ", simplex, " with sorted_A ", sorted_id_A, " and sorted_B ", sorted_id_B)
		in_A = True
		in_B = True
		new_col_A = []
		new_col_B = []
		for v in og_verts:
			if v not in A_verts:
				in_A = False
				break
			if v not in B_verts:
				in_B = False
				break
		zero_col_A = False
		zero_col_B = False
		lowest_one_in_A = False
		lowest_one_in_B = False
		if kicr_A.decomposition_im.r_data[sorted_id_A] == []:
			zero_col_A = True
			# print("zero column in R_im^A")
			cur_column = kicr_A.decomposition_im.v_data[sorted_id_A]
			# print(cur_column)
			new_col_A = []
			for s in cur_column:
				new_col_A.append(sorted_ids_AB[tuple(convert_vertices(kicr_A.fil_K.simplices()[kicr_A.new_order_to_old()[s]].vertices, inv_A))])
			new_col_A.sort()
			# print(new_col_A)
		if not zero_col_A:
			# print("non-zero column in R_im^A")
			lowest_one = kicr_A.new_order_to_old()[kicr_A.decomposition_im.r_data[sorted_id_A][-1]]
			# print("lowest one in R_im^A is ", lowest_one)
			lowest_one_vertices = convert_vertices(kicr_A.fil_K.simplices()[lowest_one].vertices, inv_A)
			lowest_one_in_A = True
			for v in lowest_one_vertices:
				if v not in A_verts:
					lowest_one_in_A = False
					break
			if lowest_one_in_A:
				# print("lowest one in R_im^A is in A")
				cur_column = kicr_A.decomposition_im.v_data[sorted_id_A]
				# print(cur_column)
				new_col_A = []
				for s in cur_column:
					new_col_A.append(sorted_ids_AB[tuple(convert_vertices(kicr_A.fil_K.simplices()[kicr_A.new_order_to_old()[s]].vertices, inv_A))])
				new_col_A.sort()
			else: 
				continue
				# print("lowest one in R_im^A is not in A")	
		# print("sorted_id_B is ", sorted_id_B)
		if kicr_B.decomposition_im.r_data[sorted_id_B] != []:
			# print("non-zero column in R_im^B")
			lowest_one = kicr_A.new_order_to_old()[kicr_B.decomposition_im.r_data[sorted_id_B][-1]]
			# print("lowest one in R_im^B is ", lowest_one)
			lowest_one_vertices = convert_vertices(kicr_B.fil_K.simplices()[lowest_one].vertices, inv_B)
			lowest_one_in_B = True
			for v in lowest_one_vertices:
				if v not in B_verts:
					lowest_one_in_B = False
					break
			if lowest_one_in_B:
				# print("lowest one in R_im^B is in B")
				cur_column = kicr_B.decomposition_im.v_data[sorted_id_B]
				# print(cur_column)
				new_col_B = []
				for s in cur_column:
					new_col_B.append(sorted_ids_AB[tuple(convert_vertices(kicr_B.fil_K.simplices()[kicr_B.new_order_to_old()[s]].vertices, inv_B))])
				new_col_B.sort()
			else:
				# print("lowest one in R_im^B is not in B")	
				continue
		# print("zero_col_A is ", zero_col_A)
		# print("zero_col_B is ", zero_col_B)
		if new_col_A != [] and new_col_B != []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			if in_A and not in_B:
				# print("in_A and not in_B")
				d_phi.append(new_col_A)
				sorted_ids.append(i)
				d_phi.append(new_col_B)
				sorted_ids.append(i)
			elif not in_A and in_B:
				# print("not in_A and in_B")
				d_phi.append(new_col_B)
				sorted_ids.append(i)
				d_phi.append(new_col_A)
				sorted_ids.append(i)
		elif new_col_A != [] and new_col_B == []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			d_phi.append(new_col_A)
			sorted_ids.append(i)
		elif new_col_A == [] and new_col_B != []:
			# print("new_col_A is ", new_col_A)
			# print("new_col_B is ", new_col_B)
			d_phi.append(new_col_B)
			sorted_ids.append(i)
	return d_phi, sorted_ids

