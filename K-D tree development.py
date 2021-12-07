    # -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:57:53 2021

@author: lukab
"""


#https://medium.com/sicara/fast-custom-knn-sklearn-cython-de92e5a325c
import pyximport
pyximport.install()
import numpy as np
from sklearn.neighbors import BallTree
from pykdtree.kdtree import KDTree
import timeit
from cython_metric import *
print('Balltree - Native Python')
threshold_distance = 100
X = np.array([(0,0,2.0001,1),(0,0,0,0),(0,0,10.1,1),(0,0,14.1,1),(0,0,4,1),(0,0,10.1,0),(0,0,14.1,0),(0,0,2,0)])
# Redefine function with if then loops?
def euclidean_distance_different_labels(x, y):
    return (1/(np.sqrt(np.sum((x[:-1]-y[:-1])**2)))) * ((x[-1:] != y[-1:]) * (-1) + (x[-1:] == y[-1:])*(1))

def euclidean_distance_different_labels_new(x, y):
    if x[-1:] != y[-1:]:
        return np.sqrt(np.sum((x[:-1]-y[:-1])**2))
    else:
        return 50000
    
tree = BallTree(X, leaf_size=2,metric='pyfunc', func=mydist)
# Find smallest positive number - that is the closest neighbor between the point sets
start = timeit.default_timer()
dist, ind = tree.query(X, k=2)  
dist[np.isneginf(dist)] = 0
minimum_distance_index = np.where(dist ==np.amin(dist))
print(ind)
print(dist)
list_of_coordinates = list(zip(minimum_distance_index[0], minimum_distance_index[1]))
print(list_of_coordinates)
true_minimum_distance =-1 * (1/ dist[list_of_coordinates[0][0],list_of_coordinates[0][1]])
print("The minimum distance between two different neurons is: " + str(np.round(true_minimum_distance,4)))

point_1 = X[ind[list_of_coordinates[0][0],list_of_coordinates[0][1]]]
point_2 = X[ind[list_of_coordinates[1][0],list_of_coordinates[1][1]]]
end = timeit.default_timer()
print('And it is between points ' + str(list(point_1[:-1]))+ ' and ' + str(list(point_2[:-1])))
print("Time for finding distance with the tree:" + str(end-start))


# print('KD Tree')
# #np.seterr(divide='ignore')
# from sklearn.neighbors import KDTree
# threshold_distance = 50
# X = np.array([(0,0,2.01),(0,0,0),(0,0,10.1),(0,0,14.1),(0,0,4),(0,0,10.1),(0,0,14.1),(0,0,2)])
# X_mask = np.array([1,0,0,0,0,0,0,0])
# print(X.shape)
# print(X_mask.shape)
# tree = KDTree(X, leaf_size=2)
# # Find smallest positive number - that is the closest neighbor between the point sets
# start = timeit.default_timer()
# dist, ind = tree.query(X, k=2)
# print(ind)
# print(dist)

# def apply_mask_indices(ind):
#     filtered_ind = []
#     for x in np.nditer(ind.T.copy(), flags=['external_loop'], order='F'):
#        if np.not_equal(X_mask[x[0]],X_mask[x[1]]):
#            #print(x)
#            filtered_ind.append(x[1])
#     return np.array(filtered_ind)
# #def apply_filt_ind_distances(dist,filtered):
    
# filtered_indices = apply_mask_indices(ind)
# print(filtered_indices)
# for x in np.nditer(filtered_indices):
#     print(x)
# print("Minimum distance is:")
# print(np.max(dist[filtered_indices][:,1]))


# # point_1 = X[ind[list_of_coordinates[0][0],list_of_coordinates[0][1]]]
# # point_2 = X[ind[list_of_coordinates[1][0],list_of_coordinates[1][1]]]
# end = timeit.default_timer()
# # print('And it is between points ' + str(list(point_1[:-1]))+ ' and ' + str(list(point_2[:-1])))
# print("Time for finding distance with the tree:" + str(end-start))