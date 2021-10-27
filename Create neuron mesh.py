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
mypath = 'mALT(L)_dataset_small/CSV/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
volumes = {}
areas = {}
bounding_boxes = {}
times_loading = {}
segments = {}
count = 0
for bodyid in onlyfiles:
    start = timeit.default_timer()
    print('Count: ' + str(count) + ' of ' + str(len(onlyfiles)))
    neuron_df = load_neuron_data(mypath+bodyid)
    if neuron_df == 0 :
        count += 1
    else:
        number_of_segments = len(neuron_df[0].keys()) + len(neuron_df[2].keys())
        print('The number of segments is ' + str(number_of_segments))
        print('Neuron loaded')
        scaling = 1
        neuron_mesh = create_neuron_mesh(neuron_df,30,scaling)
        segments[bodyid] = number_of_segments
        print("Mesh created")
        volumes[bodyid] = neuron_mesh.volume()
        bounding_boxes[bodyid] =  neuron_mesh.bounds()
        areas[bodyid] = neuron_mesh.area()
        print('Volume of the neuron:' +str(volumes[bodyid]))
        print('Area of the neuron:' +str(areas[bodyid]))
        print('Bounds of the neuron:' +str(bounding_boxes[bodyid]))
        print('Saving neuron')

        neuron_mesh.write(bodyid    + ".ply")
        end = timeit.default_timer()
        print('Elapsed time: ' + str(end - start))
        times_loading[bodyid] = end-start
        count += 1


with open("mALT(L)_volumes.pkl", "wb") as h:
    pickle.dump(volumes, h)

with open("mALT(L)_bounding_boxes.pkl", "wb") as h:
    pickle.dump(bounding_boxes, h)

with open("mALT(L)_areas.pkl", "wb") as h:
    pickle.dump(areas, h)

with open("mALT(L)_load_times.pkl", "wb") as h:
    pickle.dump(times_loading, h)

with open("mALT(L)_number_segments.pkl", "wb") as h:
    pickle.dump(segments, h)
