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
import math
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
            end_point = True
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

def make_cylinders_points(radius,coordinates,number_vertices,a,end_point):
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
        tube_list.append(mesh_coordinates)
 
    
    return tube_list

def make_branch_cylinders_points(radius,coordinates,number_vertices,a,end_points):
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
        tube_list.append(mesh_coordinates)
    return tube_list

def create_neuron_point_cloud(neuron_properties,number_vertices,a):
    tube_points_coordinates = neuron_properties[0]
    tube_points_radii= neuron_properties[1]
    branching_points_coordinates= neuron_properties[2]
    branching_points_radii = neuron_properties[3]
    tube_points_ends = neuron_properties[4]
    branching_points_ends = neuron_properties[5]
    number_vertices = 5
    a = 1
    tube_mesh_list = []
    for key in tube_points_coordinates.keys():
        end_point = True
        coordinates = tube_points_coordinates[key]
        radii = tube_points_radii[key]
        cylinders_points = make_cylinders_points(radii , coordinates,number_vertices,a,end_point)
        tube_mesh_list.append(cylinders_points)
    
    point_list = []
    for i in range(0,len(tube_mesh_list)):
        for j in range(0,len(tube_mesh_list[i])):
            for k in range(0,len(tube_mesh_list[i][j])):
                point_list.append(tube_mesh_list[i][j][k])
    branch_mesh_list = []
    for key in branching_points_coordinates.keys():
        coordinates = branching_points_coordinates[key]
        radius = branching_points_radii[key]
        if key in branching_points_ends.keys():
            end_points = branching_points_ends[key]
            branch =  make_branch_cylinders_points(radius, coordinates,number_vertices,a,end_points)
        else:
            branch = make_branch_cylinders_points(radius, coordinates,number_vertices,a,end_points)
        branch_mesh_list.append(branch)

    for i in range(0,len(branch_mesh_list)):
        for j in range(0,len( branch_mesh_list[i])):
            for k in range(0,len( branch_mesh_list[i][j])):
                point_list.append(branch_mesh_list[i][j][k])    
    return point_list


def load_neuron_data(neuron_path):
    body_id_skeleton = pd.read_csv(neuron_path)
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
            #CHECK RADIUS
            branching_points_radii[row_id_1].append(branch_rows.iloc[0]['radius_1'])
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
    
    
def segmentwise_length(neuron_df):
    total_length = 0
    for tube_key in neuron_df[0]:
        coordinate_start =neuron_df[0][tube_key]
        coordinate_end = neuron_df[0][tube_key]
        length = math.sqrt((coordinate_start[0][0]-coordinate_end[1][0])**2 + (coordinate_start[0][1] - coordinate_end[1][1])**2 + (coordinate_start[0][2]-coordinate_end[1][2])**2)
        total_length += length
    for branch_key in neuron_df[2]:
        coordinate_start =neuron_df[2][branch_key]
        coordinate_end = neuron_df[2][branch_key]
        length = math.sqrt((coordinate_start[0][0]-coordinate_end[1][0])**2 + (coordinate_start[0][1] - coordinate_end[1][1])**2 + (coordinate_start[0][2]-coordinate_end[1][2])**2)
        total_length += length
    return total_length


def segmentwise_area(neuron_df):
    total_area = 0
    for tube_key in neuron_df[0]:
        coordinate_start =neuron_df[0][tube_key]
        coordinate_end = neuron_df[0][tube_key]
        radius_start = neuron_df[1][tube_key][0]
        radius_end = neuron_df[1][tube_key][1]
        area = np.pi*(radius_start + radius_end) * math.sqrt((radius_start-radius_end)**2 + abs(coordinate_start[0][2]-coordinate_end[1][2])**2)
        total_area += area
    for branch_key in neuron_df[2]:
        coordinate_start =neuron_df[2][branch_key]
        coordinate_end = neuron_df[2][branch_key]
        radius_start = neuron_df[3][branch_key][0]
        radius_end = neuron_df[3][branch_key][1]
        area = np.pi*(radius_start + radius_end) * math.sqrt((radius_start-radius_end)**2 + abs(coordinate_start[0][2]-coordinate_end[1][2])**2)
        total_area+= area
    return total_area

def segmentwise_volume(neuron_df):
    total_volume = 0
    for tube_key in neuron_df[0]:
        coordinate_start =neuron_df[0][tube_key]
        coordinate_end = neuron_df[0][tube_key]
        radius_start = neuron_df[1][tube_key][0]
        radius_end = neuron_df[1][tube_key][1]
        #Height = z, but this can and should be corrected further with rotated circles
        volume =(1/3)* np.pi*abs(coordinate_start[0][2]-coordinate_end[1][2]) * (radius_start**2 + radius_end**2 + radius_start*radius_end)
        total_volume += volume
    for branch_key in neuron_df[2]:
        coordinate_start =neuron_df[2][branch_key]
        coordinate_end = neuron_df[2][branch_key]
        radius_start = neuron_df[3][branch_key][0]
        radius_end = neuron_df[3][branch_key][1]
        volume =(1/3)* np.pi*abs(coordinate_start[0][2]-coordinate_end[1][2]) * (radius_start**2 + radius_end**2 + radius_start*radius_end)
        total_volume += volume
    return total_volume

    
def longest_length(neuron_df):
    starting_index = 0
    ending_indices = neuron_df[-2] + list(neuron_df[-1].keys())[-5:]
    list_of_lengths = []
    coordinate_start = neuron_df[0][0]
    for ending_index in ending_indices:
        if ending_index in neuron_df[0].keys():
            coordinate_end = neuron_df[0][ending_index]
        else:
            coordinate_end = neuron_df[2][ending_index]
        length = math.sqrt((coordinate_start[0][0]-coordinate_end[1][0])**2 + (coordinate_start[0][1] - coordinate_end[1][1])**2 + (coordinate_start[0][2]-coordinate_end[1][2])**2)
        list_of_lengths.append(length)
    return np.max(list_of_lengths)


def euclidean_distance(a_1,a_2):
    dst = math.sqrt((a_1[0]- a_2[0])**2 + (a_1[1]- a_2[1])**2 + (a_1[2] - a_2[2])**2)
    return dst
    
def bounding_box_overlap_check(bounds_1,bounds_2):
    x_1 = bounds_1[0]
    x_2 = bounds_1[1]
    y_1 = bounds_1[2]
    y_2 = bounds_1[3]
    z_1 = bounds_1[4]
    z_2 = bounds_1[5]
    
    a_1 = bounds_2[0]
    a_2 = bounds_2[1]
    b_1 = bounds_2[2]
    b_2 = bounds_2[3]
    c_1 = bounds_2[4]
    c_2 = bounds_2[5]
    
    x_bool = False
    y_bool = False
    z_bool = False
    if (x_1<= a_1 and a_1<x_2) or (x_1<= a_2 and a_2<x_2) or (a_1<= x_1 and x_1<a_2) or (a_1<= x_1 and x_1<a_2):
        x_bool = True
    if (y_1<= b_1 and b_1<y_2) or (y_1<= b_2 and b_2<y_2) or (b_1<= y_1 and y_1<b_2) or (b_1<=y_1 and y_1<b_2):
        y_bool = True
    if (z_1<= c_1 and c_1<z_2) or (z_1<= c_2 and c_2<z_2) or (c_1<= z_1 and z_1<c_2) or (c_1<= z_1 and z_1<c_2):
        z_bool = True
    if (x_bool == True and y_bool == True and z_bool == True):
        return True
    else:
        return False


def distance_between_point_groups(a_list_1,a_list_2):
    distance_list = []
    for a_1 in a_list_1:
        for a_2 in a_list_2:
            distance_list.append(euclidean_distance(a_1,a_2))
    return distance_list

    
    
def create_octree(pts,bounds):
    x_list = []
    for point in pts:
        x_list.append(point[0])
    #x_half = np.median(x_list)
    x_half = (bounds[1] - bounds[0])/2
    first_half = [i for i in pts if i[0] <=x_half]
    second_half = [i for i in pts if i[0] > x_half]

    first_half_y_list = []
    second_half_y_list = []
    for points in first_half:
        first_half_y_list.append(points[1])
    for points in second_half:
        second_half_y_list.append(points[1])
        
    #y_half_1 = np.median(first_half_y_list)
    #y_half_2 = np.median(second_half_y_list)
    y_half_1 = (bounds[3] - bounds[2])/2
    y_half_2 = (bounds[3] - bounds[2])/2  
    first_half_1 = [i for i in first_half  if i[1] <= y_half_1]
    first_half_2 = [i for i in first_half if i[1] > y_half_1]
    second_half_1 = [i for i in second_half if i[1] <= y_half_2 ]
    second_half_2 = [i for i in second_half  if i[1] > y_half_2 ]

    first_half_1_z_list = []
    first_half_2_z_list = []
    second_half_1_z_list = []
    second_half_2_z_list = []
    for points in first_half_1:
        first_half_1_z_list.append(points[2])
    for points in first_half_2:
        first_half_2_z_list.append(points[2])
    for points in second_half_1 :
        second_half_1_z_list.append(points[2])
    for points in second_half_2:
        second_half_2_z_list.append(points[2])
        
    #z_half_1 = np.median(first_half_1_z_list)
    #z_half_2 = np.median(first_half_2_z_list)
    #z_half_3 = np.median(second_half_1_z_list)
    #z_half_4 = np.median(second_half_2_z_list)
    
    z_half_1 = (bounds[5] - bounds[4])/2  
    z_half_2 = (bounds[5] - bounds[4])/2  
    z_half_3 = (bounds[5] - bounds[4])/2  
    z_half_4 = (bounds[5] - bounds[4])/2   
    point_list_color = []
    octree = {}
    octree['Levels'] = 1
    
    for i in range(1,9):
        octree[i] = {}
        octree[i]['Point_list'] = []
    octree[1]['Bounds'] = (bounds[0],x_half ,bounds[2] , y_half_1,bounds[4] , z_half_1)
    octree[2]['Bounds'] = (bounds[0],x_half ,bounds[2] , y_half_1,z_half_1 , bounds[5])
    octree[3]['Bounds'] = (bounds[0],x_half ,y_half_1 , bounds[3],bounds[4] , z_half_2)
    octree[4]['Bounds'] = (bounds[0],x_half ,y_half_1 , bounds[3],z_half_2 , bounds[5])
    octree[5]['Bounds'] = (x_half, bounds[1],bounds[2] , y_half_2,bounds[4] , z_half_3)
    octree[6]['Bounds'] = (x_half, bounds[1] ,bounds[2] , y_half_2,z_half_3 , bounds[5])
    octree[7]['Bounds'] = (x_half, bounds[1],y_half_2 , bounds[3],bounds[4] , z_half_4)
    octree[8]['Bounds'] = (x_half, bounds[1] ,y_half_2 , bounds[3],z_half_4 , bounds[5])
    for point in pts:
        if point[0] <= x_half:
            if point[1] <= y_half_1:
                if point[2] <= z_half_1:
                    octree[1]['Point_list'].append(point)
                else: #if > z_half_1
                    octree[2]['Point_list'].append(point)
            else: # if > y_half_1
                if point[2] <= z_half_2:
                    octree[3]['Point_list'].append(point)
                else: #if > z_half_2
                    octree[4]['Point_list'].append(point)
        else: # if > x_half
            if point[1] <= y_half_2:
                if point[2] <= z_half_3:
                    octree[5]['Point_list'].append(point)
                else: #if > z_half
                    octree[6]['Point_list'].append(point)
            else: # if > y_half
                if point[2] <= z_half_4:
                    octree[7]['Point_list'].append(point)
                else: #if > z_half                
                    octree[8]['Point_list'].append(point)
    return octree

def create_second_level_subtrees(octree):
    subtree_1  = create_octree(octree[1]['Point_list'],octree[1]['Bounds'])
    subtree_2  = create_octree(octree[2]['Point_list'],octree[2]['Bounds'])
    subtree_3 = create_octree(octree[3]['Point_list'],octree[3]['Bounds'])
    subtree_4 = create_octree(octree[4]['Point_list'],octree[4]['Bounds'])
    subtree_5 = create_octree(octree[5]['Point_list'],octree[5]['Bounds'])
    subtree_6  = create_octree(octree[6]['Point_list'],octree[6]['Bounds'])
    subtree_7  = create_octree(octree[7]['Point_list'],octree[7]['Bounds'])
    subtree_8 = create_octree(octree[8]['Point_list'],octree[8]['Bounds'])
    subtrees = [subtree_1,subtree_2,subtree_3,subtree_4,subtree_5,subtree_6,subtree_7,subtree_8]
    return subtrees

def distances_between_two_boxes(bounds_1,bounds_2):
    x_side_1 = abs(bounds_1[1] - bounds_1[0])
    x_side_2 = abs(bounds_2[1] - bounds_2[0])
    y_side_1 = abs(bounds_1[3] - bounds_1[2])
    y_side_2 = abs(bounds_2[3] - bounds_2[2])
    z_side_1 = abs(bounds_1[5] - bounds_1[4])
    z_side_2 = abs(bounds_2[5] - bounds_2[4])
    half_diag_1 = math.sqrt((x_side_1/2)**2 + (y_side_1/2)**2 + (z_side_1/2)**2) 
    half_diag_2 = math.sqrt((x_side_2/2)**2 + (y_side_2/2)**2 + (z_side_2/2)**2)
    center_1  = ((bounds_1[0] + bounds_1[1])/2,(bounds_1[2] + bounds_1[3])/2,(bounds_1[4] + bounds_1[5])/2)
    center_2  = ((bounds_2[0] + bounds_2[1])/2,(bounds_2[2] + bounds_2[3])/2,(bounds_2[4] + bounds_2[5])/2)
    center_to_center_distance = math.sqrt((center_1[0]-center_2[0])**2 + (center_1[1]-center_2[1])**2 + (center_1[2]-center_2[2])**2)
    theoretical_closest_distance = center_to_center_distance - half_diag_1 - half_diag_2
    theoretical_farthest_distance = center_to_center_distance + half_diag_1 + half_diag_2
    return {'center_to_center':center_to_center_distance,'theoretical_close_distance':theoretical_closest_distance,'theoretical_far_distance':theoretical_farthest_distance}


def euclidean_distance(point_1,point_2):
    
    return math.sqrt((point_1[0]-point_2[0])**2 + (point_1[1]-point_2[1])**2 + (point_1[2]-point_2[2])**2)

def relevant_close_boxes(starting_bounds, second_set_of_bounds):
    box_dictionary_list = []
    for bounds in second_set_of_bounds:
        box_dictionary_list.append(distances_between_two_boxes(starting_bounds,bounds))
    
    far_distances_dict = {i:box_dictionary_list[i]['theoretical_far_distance'] for i in range(0,len(box_dictionary_list))}
    close_distances_dict = {i:box_dictionary_list[i]['theoretical_close_distance'] for i in range(0,len(box_dictionary_list))}
    #sort dictionary by far, to start with best candidates
    #Closests distances from box 1')

    #Boxes which farthests distances are closer than closest distances of other boxes')
    
    lengths_far_better_than_close = {}
    indices_far_better_than_close = {}
    for i in range(0,len(far_distances_dict)):
        #start = far_distances_dict[i]
        indices_far_better_than_close[i] = []
    
        for j in range(0,len(close_distances_dict)):
            if far_distances_dict[i] < close_distances_dict[j]:
                indices_far_better_than_close[i].append(j)
        lengths_far_better_than_close[i] = len(indices_far_better_than_close[i])
        if lengths_far_better_than_close[i] == len(close_distances_dict) - 1:
            return [second_set_of_bounds[i]]
    maximum_length = max(lengths_far_better_than_close.values())
    selected_boxes = []
    for key in lengths_far_better_than_close.keys():
        if lengths_far_better_than_close[key] == maximum_length:
            selected_boxes.append(second_set_of_bounds[i])
    return selected_boxes

def relevant_boxes_two_sets(first_set_of_bounds,second_set_of_bounds):
    relevant_boxes_second_set = []

    for i in range(0,len(first_set_of_bounds)):
        relevant_boxes_second_set.append(relevant_close_boxes(first_set_of_bounds[i], second_set_of_bounds))

    relevant_boxes_first_set = []
    for i in range(0,len(second_set_of_bounds)):
        relevant_boxes_first_set.append(relevant_close_boxes(second_set_of_bounds[i], first_set_of_bounds))
    return relevant_boxes_first_set, relevant_boxes_second_set

        
def center_distances_between_two_boxes(bounds_1,bounds_2):
    x_side_1 = abs(bounds_1[1] - bounds_1[0])
    x_side_2 = abs(bounds_2[1] - bounds_2[0])
    y_side_1 = abs(bounds_1[3] - bounds_1[2])
    y_side_2 = abs(bounds_2[3] - bounds_2[2])
    z_side_1 = abs(bounds_1[5] - bounds_1[4])
    z_side_2 = abs(bounds_2[5] - bounds_2[4])
    center_1  = ((bounds_1[0] + bounds_1[1])/2,(bounds_1[2] + bounds_1[3])/2,(bounds_1[4] + bounds_1[5])/2)
    center_2  = ((bounds_2[0] + bounds_2[1])/2,(bounds_2[2] + bounds_2[3])/2,(bounds_2[4] + bounds_2[5])/2)
    center_to_center_distance = math.sqrt((center_1[0]-center_2[0])**2 + (center_1[1]-center_2[1])**2 + (center_1[2]-center_2[2])**2)
    return center_to_center_distance
def center_distances_between_two_boxes(bounds_1,bounds_2):
    x_side_1 = abs(bounds_1[1] - bounds_1[0])
    x_side_2 = abs(bounds_2[1] - bounds_2[0])
    y_side_1 = abs(bounds_1[3] - bounds_1[2])
    y_side_2 = abs(bounds_2[3] - bounds_2[2])
    z_side_1 = abs(bounds_1[5] - bounds_1[4])
    z_side_2 = abs(bounds_2[5] - bounds_2[4])
    center_1  = ((bounds_1[0] + bounds_1[1])/2,(bounds_1[2] + bounds_1[3])/2,(bounds_1[4] + bounds_1[5])/2)
    center_2  = ((bounds_2[0] + bounds_2[1])/2,(bounds_2[2] + bounds_2[3])/2,(bounds_2[4] + bounds_2[5])/2)
    center_to_center_distance = math.sqrt((center_1[0]-center_2[0])**2 + (center_1[1]-center_2[1])**2 + (center_1[2]-center_2[2])**2)
    return center_to_center_distance


def create_octree(neuron_bounds):
    if type(neuron_bounds) == tuple:
        bounds = neuron_bounds
        x_half = (bounds[1] + bounds[0])/2
        y_half = (bounds[3] + bounds[2])/2
        z_half = (bounds[5] + bounds[4])/2
        box_1 = (bounds[0],x_half ,bounds[2] , y_half,bounds[4] , z_half)
        box_2 = (bounds[0],x_half ,bounds[2] , y_half,z_half , bounds[5])
        box_3  = (bounds[0],x_half ,y_half , bounds[3],bounds[4] , z_half)
        box_4  = (bounds[0],x_half ,y_half , bounds[3],z_half ,bounds[5])
        box_5 = (x_half, bounds[1], bounds[2], y_half, bounds[4], z_half)
        box_6  = (x_half, bounds[1], bounds[2], y_half, z_half, bounds[5])
        box_7 = (x_half, bounds[1],y_half , bounds[3],bounds[4] , z_half)
        box_8  = (x_half, bounds[1] ,y_half , bounds[3],z_half , bounds[5])
        bounds_list = [box_1,box_2,box_3,box_4,box_5,box_6,box_7,box_8]
    else:
        bounds_list = []
        for bounds in neuron_bounds:
            x_half = (bounds[1] + bounds[0])/2
            y_half = (bounds[3] + bounds[2])/2
            z_half = (bounds[5] + bounds[4])/2
            box_1 = (bounds[0],x_half ,bounds[2] , y_half,bounds[4] , z_half)
            box_2 = (bounds[0],x_half ,bounds[2] , y_half,z_half , bounds[5])
            box_3  = (bounds[0],x_half ,y_half , bounds[3],bounds[4] , z_half)
            box_4  = (bounds[0],x_half ,y_half , bounds[3],z_half ,bounds[5])
            box_5 = (x_half, bounds[1], bounds[2], y_half, bounds[4], z_half)
            box_6  = (x_half, bounds[1], bounds[2], y_half, z_half, bounds[5])
            box_7 = (x_half, bounds[1],y_half , bounds[3],bounds[4] , z_half)
            box_8  = (x_half, bounds[1] ,y_half , bounds[3],z_half , bounds[5])
            box_list = [box_1,box_2,box_3,box_4,box_5,box_6,box_7,box_8]
            for box in box_list:
                bounds_list.append(box)
    return bounds_list

def total_bounds(bounding_box_list):
    x_list = []
    y_list = []
    z_list = []
    for bounds in bounding_box_list:
        x_list.append(bounds[0])
        x_list.append(bounds[1])
        y_list.append(bounds[2])
        y_list.append(bounds[3])
        z_list.append(bounds[4])
        z_list.append(bounds[5])
    return (min(x_list),max(x_list),min(y_list),max(y_list),min(z_list),max(z_list))