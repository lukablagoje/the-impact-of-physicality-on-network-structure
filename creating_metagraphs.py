import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import pickle
import timeit
import itertools

number_neurons = {'ME(R)': 3721}
roi_list = ['ME(R)']


df_dict = {}
for neuron_df_key in roi_list:
    skeleton_df = pd.read_csv(neuron_df_key+'_'+str(number_neurons[neuron_df_key]) + '.csv')
    df_dict[neuron_df_key] = skeleton_df
    
    
    
points_df_dict = {}
for region_key in df_dict.keys():
    points_df = df_dict[region_key]
    points_df_dict[region_key] = points_df[['x','y','z','radius','bodyId']].astype('float')
    
    
# Takes a single row with 3D coordinates local radius/diameter and bodyid. Returns 14 points around it, without radius.
def spherical_sample_dataframe_row(array):
        results_array = np.empty((0,4),np.float64)
        point = np.array([array[0],array[1],array[2],array[4]])
        #if diameter_or_radius == "diameter":
        # If input is diameter, divide by 2 to get radius
        r = array[3]
        #elif diameter_or_radius == "radius":
        #    r = array[3]
        points_array = [point + np.array([0,0,r,0]),point + np.array([0,0,-r,0]),point + np.array([0,r,0,0]),point + np.array([0,-r,0,0]),
                        point + np.array([r,0,0,0]),point + np.array([-r,0,0,0]),point + np.array([r * 1/2,r * 1/2,r* np.sqrt(2)/2,0]),point + np.array([r * 1/2,-r* 1/2,r* np.sqrt(2)/2,0]),
                        point + np.array([-r* 1/2,r* 1/2,r* np.sqrt(2)/2,0]),point + np.array([-r* 1/2,-r* 1/2,r* np.sqrt(2)/2,0]),point + np.array([r* 1/2,r *1/2,-r* np.sqrt(2)/2,0]),
                        point + np.array([r* 1/2,-r* 1/2,-r* np.sqrt(2)/2,0]), point + np.array([-r* 1/2,r* 1/2,-r* np.sqrt(2)/2,0]),point + np.array([-r* 1/2,-r* 1/2,-r* np.sqrt(2)/2,0])
                                    ]
        results_array =   np.append(results_array,np.array(points_array),axis=0)
        return np.array(results_array).astype(np.float64)

    
def dataframe_to_bodyid_kd_tree(points_df,spherical):
    start_total = timeit.default_timer()
    if spherical == True:
        points_numpy = points_df.to_numpy()
        points_numpy = np.apply_along_axis(spherical_sample_dataframe_row,1,points_numpy )
        points_numpy = np.concatenate(points_numpy)
    else:
        points_df_non_spherical = points_df[['x','y','z','bodyId']]
        points_numpy = points_df_non_spherical.to_numpy()
    end_total = timeit.default_timer()
    print('Loading time:',end_total-start_total)
    body_id_array, occurence_start_array = np.unique(points_numpy[:,3],return_index = True)
    occurence_start_array = np.sort(occurence_start_array)
    big_tree = KDTree(points_numpy[:,:3],compact_nodes=True,balanced_tree=True,leafsize = 400)
    end_total = timeit.default_timer()
    print('KD Tree construction time:',end_total-start_total)
    return big_tree,body_id_array,occurence_start_array,points_numpy


def query_tree_with_neuron_bodyid(big_tree,occurence_start_array,body_id_array,points_numpy,radius):
    metagraph_edgelist_dictionary = {}
    for i in range(0,occurence_start_array.shape[0]):
            body_id = int(points_numpy[:,3][occurence_start_array[i]])
            if i <occurence_start_array.shape[0]-1:
                interval_start_index = occurence_start_array[i]
                interval_end_index = occurence_start_array[i+1]
                results = big_tree.query_ball_point(points_numpy[interval_start_index :interval_end_index, :3], r = radius,return_sorted=True)
                results_array = np.fromiter(itertools.chain.from_iterable(results), np.int64)
            else:
                interval_start_index = occurence_start_array[i]
                results = big_tree.query_ball_point(points_numpy[interval_start_index :, :3], r = radius,return_sorted=True)
                results_array = np.fromiter(itertools.chain.from_iterable(results), np.int64)
            unique_results = np.unique(results_array)
            metagraph_edgelist_dictionary[body_id] = []
            # Check for each neuron if it's in the neighborhood
            for j in range(0,occurence_start_array.shape[0]):
                    # Remove the neuron we are looking at 
                    body_id_neighbor =  int(points_numpy[:,3][occurence_start_array[j]])
                    if body_id != body_id_neighbor:
                        if j <occurence_start_array.shape[0]-1:
                            interval_start_index = occurence_start_array[j]
                            interval_end_index = occurence_start_array[j+1]
                            #print('BodyID of the testing neighbor neuron:',body_id_neighbor)
                            #print('Starting index of the neuron in the large array:',interval_start_index)
                            #print('Ending index of the neuron in the large array:',interval_end_index)
                            mask_inverted = np.logical_and((interval_start_index  <= unique_results),(unique_results < interval_end_index))
                        else:
                            interval_start_index = occurence_start_array[-1]
                            #print('BodyID of the testing neighbor neuron:',body_id_neighbor)
                            #print('Starting index of the neuron in the large array:',interval_start_index)
                            #print('Ending index of the neuron in the large array:',)
                            mask_inverted = np.logical_and((interval_start_index  <= unique_results),True)
                        statement = mask_inverted.any()
                        #If the neighbor with this bodyid is present in the neighborhood
                        if mask_inverted.any():
                            metagraph_edgelist_dictionary[body_id].append(body_id_neighbor)
                        mask = np.logical_not(mask_inverted) 
                        unique_results = unique_results[mask]
                        if unique_results.size == 0: #If array is emptied, stop the process
                            break
    return metagraph_edgelist_dictionary


for region_key in list(points_df_dict.keys()):
    print(region_key)
    spherical = True
    points_df = points_df_dict[region_key]
    points_df = points_df[pd.notnull(points_df.bodyId)]
    big_tree,body_id_array,occurence_start_array,points_numpy = dataframe_to_bodyid_kd_tree(points_df,spherical)
    total_number_of_edges_spherical = []
    #step  = 0.01
    radius_list_spherical = np.arange(0,1 + step ,0.01)
    radius_list_spherical = [1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 45, 50, 60, 70, 80, 90, 100, 150,200]
    for radius in radius_list_spherical:
        start_total = timeit.default_timer()
        metagraph_edgelist_dictionary = {}
        print('Radius',radius)
        metagraph_edgelist_dictionary = query_tree_with_neuron_bodyid(big_tree,occurence_start_array,body_id_array,points_numpy,radius)
        count_edges = 0
        # Go through each neuron and count the number of it's neighbors
        for key in metagraph_edgelist_dictionary.keys():
            count_edges += len(metagraph_edgelist_dictionary[key])
        # Edges are always double counted, since i - j will always have j-i as a neighbor
        count_edges = count_edges/2
        total_number_of_edges_spherical.append(count_edges)
        end_total = timeit.default_timer()
        print('Total time taken for radius:',end_total-start_total)
        print('Count',count_edges)
        with open(region_key+"_edge_dict_spherical_" + str(radius)+".pkl", "wb") as h:
                    pickle.dump(metagraph_edgelist_dictionary, h)
