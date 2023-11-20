import oineus
import diode
import numpy
import pandas
from functools import cmp_to_key
import math

def oineus_pair(points : pandas.DataFrame, sub : list):
	"""! Given a set of points, and the points that are in the subset L, construct the complexes and map between them. The subcomplex L will consists of all simplices whose vertex are in the subset.

	@param points		pandas.DataFrame containing the points and their weights
	@param sub			a list containing the indices of the points on which we construct the subcomplex

	@return K			list of simplices for the entire complex, as needed by oineus
	@return L			list of simplices for the subcomplex, as needed by oineus
	@return L_to_K		list which tells you how to map the simplices in L to the simplices in K
	"""
	
	points["sub"]=sub
	points.sort_values(by="sub", ascending=False)
	simplices = diode.fill_weighted_alpha_shapes(points[["x","y","z","w"]].to_numpy())
	for i in range(len(simplices)):
		simplices[i] = [sorted(simplices[i][0]), simplices[i][1]]
	simplices = sorted(simplices, key=cmp_to_key(oineus_compare))
	L = []
	not_L = []
	#L_to_K = []
	#id_L = 0
	#K_to_L = [-1 for i in range(len(simplices))]
	for i,s in enumerate(simplices):
		if len(s[0])==1:
			if sub[s[0][0]]==True:
				#K_to_L[i]= id_L            
				L.append([s[0], s[1]])
				#L_to_K.append(i)
				#id_L +=1
			else:
				not_L.append([s[0],s[1]])
		else:
			sub_complex = True
			for v in s[0]:
				if sub[v] == False:
					sub_complex=False
					break
			if sub_complex == True:
				#verts = [K_to_L[v] for v in s[0]]
				L.append([s[0], s[1]])
				#L_to_K.append(i)
			  	#K_to_L[i] = id_L
				#id_L +=1
			else:
				not_L.append([s[0],s[1]])
	K = []
	for s in L:
		K.append(s)
	for s in not_L:
		K.append(s)
	L = [[i,s[0],s[1]] for i, s in enumerate(L)]
	K = [[i,s[0],s[1]] for i, s in enumerate(K)]
	return K, L#, L_to_K

def oineus_compare(x, y):
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

def sub_complex(points : pandas.DataFrame, z_upper : float, z_lower : float):
	"""! Given the points, and the upper and lower thresholds in the 'z'-component. 

	@param points		pandas.DataFrame containing of the points.
	@param z_upper		float giving the upper threshold, any point above this is in the subcomplex
	@param z_lower		float giving the lower threshold, any point below this is in the subcomplex

	@return sub_comp	list containing the indices of the points on which we build the subcomplex
	"""
	sub_comp = []
	for i in range(points.shape[0]):
		if (points["z"][i] >= z_upper) or (points["z"][i]    <= z_lower):
			sub_comp.append(True)
		else:
			sub_comp.append(False)
	return sub_comp    

def kernel_image_cokernel(points : pandas.DataFrame, kernel : bool, image : bool, cokernel : bool, n_threads : int, upper_threshold : float, lower_threshold : float):
	"""! Given points, and parameters for oineus, calculate the kernel/image/cokernel persistence as desired.

	@param points			pandas.DataFrame of the points, with columns 'x','y','z','w' corresponding to the coordinates and weights respectively
	@param kernel 			boolean parameter to set if kernel persistence is calculated
	@param image 			boolean parameter to set if image persistence is calculated
	@param cokernel 		boolean parameter to set if cokernel persistence is calculated
	@param n_threads		number of threads to use in oineus
	@param upper_threshold	float, z-coordinate above which points are in the subcomplex 
	@param lower_threshold	float z-coordinate below which points are in the subcomplex


	@return kicr			oineus object which contains the kernel, image, cokernel persistence diagrams as required, can also calculate ones that weren't initially specificed
	"""
	params = oineus.ReductionParams()
	params.n_threads = n_threads
	params.kernel = kernel
	params.image = image
	params.cokernel = cokernel
	sub = sub_complex(points, upper_threshold, lower_threshold)
	K, L = oineus_pair(points, sub)
	L = oineus.list_to_filtration_float(L, params)
	K = oineus.list_to_filtration_float(K, params)
	kicr = oineus.KerImCokReduced_float(K,L,params,False)
	return kicr


points = numpy.random.random((10,4))
points = pandas.DataFrame(points, columns=["x","y","z","w"])

kicr = kernel_image_cokernel(points, False, True, False, 4, math.floor(max(points["z"])), math.ceil(min(points["z"])))

image_pd = kicr.image_diagrams().in_dimension(1)
for pt in image_pd:
	print(pt)

