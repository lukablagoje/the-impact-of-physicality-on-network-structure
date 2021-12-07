# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:01:30 2021

@author: lukab
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
import timeit
start = timeit.default_timer()
object_1 = [0.2, 4.5, 198, 0.003]
object_2 = [0.3, 2.0, 999, 0.001]
object_3 = [0.1, 9.2, 321, 0.023]
list_of_objects = [object_1, object_2, object_3]

# make a 4x3 matrix from list of objects
X = np.array([(0,0,2.0001,1),(0,0,0,0),(0,0,10.1,1),(0,0,14.1,1),(0,0,4,1),(0,0,10.1,0),(0,0,14.1,0),(0,0,2,0)])

#calculate pairwise distances
distances = pdist(X)

#make a square matrix from result
distances_as_2d_matrix = squareform(distances)

print(distances)
print(distances_as_2d_matrix)

end = timeit.default_timer()
print("Time for finding distance with the tree:" + str(end-start))
