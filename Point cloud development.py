# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:35:54 2021

@author: lukab
"""

from vedo import *
import numpy as np
import pandas as pd

points_df = pd.read_csv('ME(R)_neurons_sample.points.csv')


point = np.array((0,0,0))
r = 1
def spherical_sample(point,r):
#radius

# point_1 = point + np.array([0,0,r])
# point_2 =  point + np.array([0,0,-r])
# point_3 = point + np.array([0,r,0])
# point_4 =  point + np.array([0,-r,0])
# point_5 = point + np.array([r,0,0])
# point_6 =  point + np.array([-r,0,0])

# point_7 = point + np.array([r * 1/2,r * 1/2,r* np.sqrt(2)/2])
# point_8 = point + np.array([r * 1/2,-r* 1/2,r* np.sqrt(2)/2])
# point_9 = point + np.array([-r* 1/2,r* 1/2,r* np.sqrt(2)/2])
# point_10 = point + np.array([-r* 1/2,-r* 1/2,r* np.sqrt(2)/2])

# point_11 = point + np.array([r* 1/2,r *1/2,-r* np.sqrt(2)/2])
# point_12 = point + np.array([r* 1/2,-r* 1/2,-r* np.sqrt(2)/2])
# point_13 = point + np.array([-r* 1/2,r* 1/2,-r* np.sqrt(2)/2])
# point_14 = point + np.array([-r* 1/2,-r* 1/2,-r* np.sqrt(2)/2])

    return np.array([point + np.array([0,0,r]),point + np.array([0,0,-r]),point + np.array([0,r,0]),point + np.array([0,-r,0]),
                        point + np.array([r,0,0]),point + np.array([-r,0,0]),point + np.array([r * 1/2,r * 1/2,r* np.sqrt(2)/2]),point + np.array([r * 1/2,-r* 1/2,r* np.sqrt(2)/2]),
                        point + np.array([-r* 1/2,r* 1/2,r* np.sqrt(2)/2]),point + np.array([-r* 1/2,-r* 1/2,r* np.sqrt(2)/2]),np.array([r* 1/2,r *1/2,-r* np.sqrt(2)/2]),
                        point + np.array([r* 1/2,-r* 1/2,-r* np.sqrt(2)/2]), point + np.array([-r* 1/2,r* 1/2,-r* np.sqrt(2)/2]),point + np.array([-r* 1/2,-r* 1/2,-r* np.sqrt(2)/2])])
point_cloud = Points(spherical_sample(point,r),r=55)
sphere = Sphere(pos=(0,0,0),r=1).alpha(0.5)

show(point_cloud,sphere,axes=2).close()