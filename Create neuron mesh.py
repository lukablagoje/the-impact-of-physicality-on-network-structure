# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 09:33:34 2021

@author: lukab
"""

from Important_functions_optimized import *
from vedo import Mesh, merge, show
import numpy as np
import pickle
import pandas as pd
import timeit

from os import listdir
from os.path import isfile, join
mypath = 'ME(R)/CSV/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
volumes = {}
areas = {}
bounding_boxes = {}
times_loading = {}
segments = {}
segmentwise_lengths = {}
shortest_lengths = {}
count = 0
segmentwise_areas = {}
segmentwise_volumes = {}

for bodyid in onlyfiles:
    start = timeit.default_timer()
    print('Count: ' + str(count) + ' of ' + str(len(onlyfiles)))
    neuron_df = load_neuron_data(mypath+bodyid)
    if neuron_df == 0 :
        count += 1
    else:
        segmentwise_areas[bodyid[:-4]] = segmentwise_area(neuron_df)
        segmentwise_volumes[bodyid[:-4]] = segmentwise_volume(neuron_df)
        number_of_segments = len(neuron_df[0].keys()) + len(neuron_df[2].keys())
        segmentwise_lengths[bodyid[:-4]] = segmentwise_length(neuron_df)
        shortest_lengths[bodyid[:-4]] = longest_length(neuron_df)
        print('The number of segments is ' + str(number_of_segments))
       
        print('Neuron loaded')
        scaling = 1
        neuron_mesh = create_neuron_mesh(neuron_df,30,scaling)
        # neuron_mesh = neuron_mesh_1.clone().decimate(fraction=0.25).color('blue')
        # segments[bodyid] = number_of_segments
        print("Mesh created")
        bounding_boxes[bodyid[:-4]] =  neuron_mesh.bounds()
        # areas[bodyid] = neuron_mesh.area()
        # print('Volume of the neuron:' +str(volumes[bodyid]))
        # print('Area of the neuron:' +str(areas[bodyid]))
        # print('Bounds of the neuron:' +str(bounding_boxes[bodyid]))
        print('Saving neuron')
        # show(neuron_mesh_1).close()
        # neuron_mesh.write(bodyid    + ".ply")
        end = timeit.default_timer()
        print('Elapsed time: ' + str(end - start))
        times_loading[bodyid[:-4]] = end-start
        count += 1


    
with open("ME(R)_segmentwise_areas.pkl", "wb") as h:
    pickle.dump(segmentwise_areas, h)
    
with open("ME(R)_segmentwise_volumes.pkl", "wb") as h:
    pickle.dump(segmentwise_volumes, h)
    
with open("ME(R)_number_of_segments.pkl", "wb") as h:
    pickle.dump(number_of_segments, h)

with open("ME(R)_segmentwise_lengths.pkl", "wb") as h:
    pickle.dump(segmentwise_lengths, h)
    
with open("ME(R)_bounding_boxes.pkl", "wb") as h:
    pickle.dump(bounding_boxes, h)
    
with open("ME(R)_times_loading.pkl", "wb") as h:
    pickle.dump(times_loading, h)
    
    # with open("mALT(L)_dec_areas.pkl", "wb") as h:
    #     pickle.dump(areas, h)
    
    # with open("mALT(L)_dec_load_times.pkl", "wb") as h:
    #     pickle.dump(times_loading, h)
    
    # with open("mALT(L)_dec_number_segments.pkl", "wb") as h:
    #     pickle.dump(segments, h)
    
    #with open("mALT(L)_dec_segmentwise_lengths.pkl", "wb") as h:
    #    pickle.dump(segmentwise_lengths, h)
    
    #with open("mALT(L)_dec_longest_lengths.pkl", "wb") as h:
    #    pickle.dump(longest_lengths, h)

    # with open("mALT(L)_segmentwise_volumes.pkl", "wb") as h:
    #    pickle.dump(segmentwise_volumes, h)
        
    # with open("mALT(L)_segmentwise_areas.pkl", "wb") as h:
    #    pickle.dump(segmentwise_areas, h)