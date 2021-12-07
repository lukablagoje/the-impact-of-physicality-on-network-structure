# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:57:53 2021

@author: lukab
"""

import numpy as np
from sklearn.neighbors import KDTree
rng = np.random.RandomState(0)
X = np.array([(0,0,0),(0,0,1),(0,0,2)]) # 10 points in 3 dimensions
mask = [1,1,0]
neuron_1_index_start = 5
neuron_1_index_end = 10
tree = KDTree(X, leaf_size=2)  
distance = 0.1     
from sklearn.neighbors import KDTree
import numpy as np
tree = KDTree(X, leaf_size=2) 
count = tree.query_radius(X[:1], r=distance, count_only=True)   
print(X)
print("Count")
print(count)


ind = tree.query_radius(X, r=distance)  
print("Indices")
print(ind) 
#print(np.argmax(ind > neuron_1_number_of_points))
if np.max(ind[0]) > neuron_1_index_start and neuron_1_index_end > np.max(ind[0]) :
    print(np.max(ind[0]))
    print('Match')
        # indices of neighbors within distance 0.3

dist = tree.query_radius(X, r=distance,return_distance=True)  
print("Distances")
print(dist)