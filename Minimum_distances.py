# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 14:23:21 2021

@author: lukab
"""


from Important_functions_optimized import *
from vedo import *
import numpy as np
import pickle
import pandas as pd
import timeit
from os import listdir
from os.path import isfile, join
region  = 'mALT(L)_'
edge_list = pd.read_csv('mALT(L)_dataset_small/edge_list_mALT(L)')
mypath = 'mALT(L)_dataset/Meshes/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
neuron_minimum_distances = {}
neuron_maximum_distances = {}
neuron_mean_distances = {}
neuron_median_distances = {}
areas = {}
times_distances = {}
for index,row in edge_list.iterrows():
    print("Index: " + str(index) + ' / ' + str(len(edge_list)))
    start = timeit.default_timer()
    name_1 = 'mALT(L)_' + str(row['bodyId_pre']) + '.ply'
    name_2 = 'mALT(L)_' +  str(row['bodyId_post']) + '.ply'
    if (name_1 in onlyfiles) and (name_2 in onlyfiles):
        mesh_1 = load(mypath + 'mALT(L)_' + str(row['bodyId_pre']) + '.ply')
        mesh_2 = load(mypath + 'mALT(L)_' +  str(row['bodyId_post']) + '.ply')
        mesh_1.distanceToMesh(mesh_2, signed=False, negate=False)
        min_distance_1 = min(mesh_1.pointdata["Distance"])
        max_distance_1 = max(mesh_1.pointdata["Distance"])
        mean_distance_1 = np.mean((mesh_1.pointdata["Distance"]))
        median_distance_1 = np.median((mesh_1.pointdata["Distance"]))
        print('Minimum distance:' + str(min_distance_1))
        print('Maximum distance:' + str(max_distance_1))
        print('Mean distance:' + str(mean_distance_1))   
        print('Median distance:' + str(median_distance_1))
        neuron_minimum_distances[(row['bodyId_pre'],row['bodyId_post'])] = min_distance_1
        neuron_maximum_distances[(row['bodyId_pre'],row['bodyId_post'])] = max_distance_1
        neuron_mean_distances[(row['bodyId_pre'],row['bodyId_post'])] = mean_distance_1
        neuron_median_distances[(row['bodyId_pre'],row['bodyId_post'])] = median_distance_1
        with open(region + "neuron_minimum_distances.pkl", "wb") as h:
            pickle.dump(neuron_minimum_distances , h)
            
        with open(region + "neuron_maximum_distances.pkl", "wb") as h:
            pickle.dump(neuron_maximum_distances , h)
            
        with open(region + "neuron_mean_distances.pkl", "wb") as h:
            pickle.dump(neuron_mean_distances , h)
            
        with open(region + "neuron_median_distances.pkl", "wb") as h:
            pickle.dump(neuron_median_distances , h)
            
        with open(region+ "areas.pkl", "wb") as h:
            pickle.dump(areas, h)
        with open(region+ "times_distances.pkl", "wb") as h:
            pickle.dump(times_distances, h)
        end = timeit.default_timer()
        times_distances[(row['bodyId_pre'],row['bodyId_post'])] = end - start
        print('Elapsed time: ' + str(end - start))
        
with open(region + "neuron_minimum_distances.pkl", "wb") as h:
    pickle.dump(neuron_minimum_distances , h)
    
with open(region + "neuron_maximum_distances.pkl", "wb") as h:
    pickle.dump(neuron_maximum_distances , h)
    
with open(region + "neuron_mean_distances.pkl", "wb") as h:
    pickle.dump(neuron_mean_distances , h)
    
with open(region + "neuron_median_distances.pkl", "wb") as h:
    pickle.dump(neuron_median_distances , h)
    
with open(region+ "areas.pkl", "wb") as h:
    pickle.dump(areas, h)
with open(region+ "times_distances.pkl", "wb") as h:
    pickle.dump(times_distances, h)