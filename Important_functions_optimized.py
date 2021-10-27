# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 12:24:42 2021

@author: lukab
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 16:59:10 2021

@author: lukab
"""

from vedo import *
from vedo import Mesh, merge, show
import numpy as np
import pickle
import pandas as pd
def make_cylinders(radius,coordinates,number_vertices,a,end_point):
    tube_list = []
    
    n = 2*number_vertices
    
    # Defining faces used in creation of the mesh
    faces = [[i, i+1, i+3, i+2] for i in range(0, n-2, 2)]
    faces.append([n-2, n-1, 1, 0])
    for i in range(0,len(coordinates)-1):
        r_start = a *radius[0]
        if end_point == True:
            r_end = 0
        else:
            r_end =a*radius[1]
        angle = 0
        number_of_vertices = number_vertices
        start_x = coordinates[i][0]
        start_y = coordinates[i][1]
        start_z = coordinates[i][2]
        end_x = coordinates[i+1][0]
        end_y = coordinates[i+1][1]
        end_z = coordinates[i+1][2]
        mesh_coordinates = []
        
        for j in range(0, number_of_vertices):
            angle = (j / number_of_vertices) * 2 * np.pi
            circle_x_start = start_x + r_start * np.cos(angle)
            circle_y_start = start_y + r_start * np.sin(angle)
            circle_x_end = end_x + r_end * np.cos(angle)
            circle_y_end = end_y + r_end * np.sin(angle)
            mesh_coordinates.append((circle_x_start, circle_y_start, start_z))
            mesh_coordinates.append((circle_x_end, circle_y_end, end_z))
        tube = Mesh([mesh_coordinates, faces]).triangulate()
        tube_list.append(tube)
 
    
    tube_first = merge(tube_list)
    tube_first = tube_first.clean()
    return tube_first

def make_branch_cylinders(radius,coordinates,number_vertices,a,end_points):
    tube_list = []
    n = 2*number_vertices
    faces = [[i, i+1, i+3, i+2] for i in range(0, n-2, 2)]
    faces.append([n-2, n-1, 1, 0])

    for i in range(0,len(coordinates)-1):
        r_start = a* radius[0]
        if i+1 in end_points:
            r_end = 0
        else:
            r_end = a * radius[i]

        angle = 0
        number_of_vertices = number_vertices
        start_x = float(coordinates[0][0])
        start_y = float(coordinates[0][1])
        start_z = float(coordinates[0][2])
        end_x = float(coordinates[i+1][0])
        end_y = float(coordinates[i+1][1])
        end_z = float(coordinates[i+1][2])
        mesh_coordinates = []
    
        for j in range(0, number_of_vertices):
            angle = (j / number_of_vertices) * 2 * np.pi
            circle_x_start = start_x + r_start * np.cos(angle)
            circle_y_start = start_y + r_start * np.sin(angle)
            circle_x_end = end_x + r_end * np.cos(angle)
            circle_y_end = end_y + r_end * np.sin(angle)
            mesh_coordinates.append((circle_x_start, circle_y_start, start_z))
            mesh_coordinates.append((circle_x_end, circle_y_end, end_z))
        tube = Mesh([mesh_coordinates, faces]).triangulate()
        tube_list.append(tube)
 
    

    tube_first = merge(tube_list)
    tube_first = tube_first.clean()
    return tube_first

def create_neuron_mesh(neuron_properties,number_vertices,a):
    tube_points_coordinates = neuron_properties[0]
    tube_points_radii= neuron_properties[1]
    branching_points_coordinates= neuron_properties[2]
    branching_points_radii = neuron_properties[3]
    tube_points_ends = neuron_properties[4]
    branching_points_ends = neuron_properties[5]
    tube_mesh_list = []
    for key in tube_points_coordinates.keys():
        coordinates = tube_points_coordinates[key]
        radii = tube_points_radii[key]
        if key in tube_points_ends:
            end_point = True
            tube = make_cylinders(radii , coordinates,number_vertices,a,end_point)
        else:
            end_point = False
            tube = make_cylinders(radii , coordinates,number_vertices,a,end_point)
        tube_mesh_list.append(tube)
    
    final_tube = merge(tube_mesh_list)

    
    branch_mesh_list = []
    for key in branching_points_coordinates.keys():
        coordinates = branching_points_coordinates[key]
        radius = branching_points_radii[key]
        if key in branching_points_ends.keys():
            end_points = branching_points_ends[key]
            branch =  make_branch_cylinders(radius, coordinates,number_vertices,a,end_points)
        else:
            branch = make_branch_cylinders(radius, coordinates,number_vertices,a,end_points)
        branch_mesh_list.append(branch)
    

    final_branch =  merge(branch_mesh_list)
    
    
    final_mesh = merge(final_branch, final_tube)
    #cap = final_mesh.cap(returnCap=True).triangulate()
    #final_tube = merge(final_mesh, cap).clean()
        #tube_cut = final_tube.clone().cutWithPlane().cap()
    final_tube = final_mesh.fillHoles(size=10*np.max(list(tube_points_radii.values()))).clean()
    #final_tube = final_mesh.clean()
    return final_tube

def load_neuron_data(neuron_path):
    body_id_skeleton = pd.read_csv(neuron_path)
    print(len(body_id_skeleton))
    if len(body_id_skeleton) > 0:
        # Finding branching points - used to create "branching cylinders"
        dictionary = dict(body_id_skeleton['rowId_1'].value_counts())
        branching_row_id_1 = {}
        
        for row_id_1 in dictionary.keys():
            branching_row_id_1[row_id_1] = 0
        for row_id_1 in dictionary.keys():
            if dictionary[row_id_1] > 1:
                branching_row_id_1[row_id_1] += dictionary[row_id_1]
            else:
                branching_row_id_1.pop(row_id_1)
                
        ending_rows = []
        for row_id_2 in body_id_skeleton['rowId_2']:
            if len(body_id_skeleton[body_id_skeleton['rowId_1'] == row_id_2]) == 0:
                ending_rows.append(row_id_2)
        
        branching_points_coordinates = {}
        branching_points_radii = {}
        branching_point_ends = {}
        for row_id_1 in branching_row_id_1.keys():
            branching_points_coordinates[row_id_1] = []
            branching_points_radii[row_id_1] = []
            branch_rows = body_id_skeleton[body_id_skeleton['rowId_1'] == row_id_1]
            branching_points_coordinates[row_id_1].append((branch_rows.iloc[0]['x_1'],branch_rows.iloc[0]['y_1'],branch_rows.iloc[0]['z_1']))
            count = 0
            branching_point_ends[row_id_1]  = []
            for index,row in branch_rows.iterrows():
                branching_points_coordinates[row_id_1].append((branch_rows.iloc[count]['x_2'],branch_rows.iloc[count]['y_2'],branch_rows.iloc[count]['z_2']))
                if count == 0:
                    branching_points_radii[row_id_1].append(branch_rows.iloc[count]['radius_1'])
                else:
                    branching_points_radii[row_id_1].append(branch_rows.iloc[count]['radius_2'])
                    if row['rowId_2'] in ending_rows:
                        branching_point_ends[row_id_1].append(count)
                count += 1
    
                    
                    
        # Completing the rest of the neuron (tubes)
        non_branch_tube_id_list = body_id_skeleton['rowId_1'].to_list()
        for row_id in body_id_skeleton['rowId_1']:
            if row_id in branching_points_coordinates.keys():
                non_branch_tube_id_list.remove(row_id)
                
                    
        tube_points_coordinates = {}
        tube_points_radii = {}
        tube_ends = [0]
        for row_id_1 in non_branch_tube_id_list:
            row = body_id_skeleton[body_id_skeleton['rowId_1'] == row_id_1]
            tube_points_coordinates[row_id_1] = []
            tube_points_radii[row_id_1] = []
            tube_points_coordinates[row_id_1].append((row['x_1'].values[0],row['y_1'].values[0],row['z_1'].values[0]))
            tube_points_coordinates[row_id_1].append((row['x_2'].values[0],row['y_2'].values[0],row['z_2'].values[0]))
            tube_points_radii[row_id_1].append(row['radius_1'].values[0])
            tube_points_radii[row_id_1].append(row['radius_2'].values[0])
            if row['rowId_2'].values[0] in ending_rows:
                tube_ends.append(row_id_1)
        starting_row =  np.min(body_id_skeleton['rowId_1'])
        row = body_id_skeleton[body_id_skeleton['rowId_1'] == starting_row]
        tube_points_coordinates[0] = []
        tube_points_radii[0] = []
        tube_points_coordinates[0].append((row['x_1'].values[0],row['y_1'].values[0],row['z_1'].values[0]))
        tube_points_coordinates[0].append((row['x_1'].values[0],row['y_1'].values[0],row['z_1'].values[0]-1))
        tube_points_radii[0].append(row['radius_1'].values[0])
        tube_points_radii[0].append(row['radius_1'].values[0])        
        
        result_list = [tube_points_coordinates,tube_points_radii,branching_points_coordinates,branching_points_radii,tube_ends,branching_point_ends]
        return result_list
    else:
        return 0