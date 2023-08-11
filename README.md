# Research overview
This project is a part of the published research work "Understanding the Impact of physicality on network structure", done in the emerging field of Physical Networks, which aims to understand the properties of three-dimensional networked systems (such as biological neural networks):
![image](https://github.com/lukablagoje/physical-networks-spatially-embedded-networks/assets/52599010/f9c2db46-4a1d-4cff-bfaf-4c9076193b91)

This research is focused on specific models and meta-graph, which encode if physical objects are in collision with each other (if they collide, there is an edge). I technically implemented a generalized meta-graph (using point clouds and k-d trees), by encoding information on how far away are the neighboring links to each other, in terms of Euclidean distance in 3D space (so it's not only a binary collision). I applied this new representation to efficiently solve a computationally challenging task: finding out how many link-to-link physical collisions are obtained if the thickness of the neurons is increased for 20 different additive factors. The results of this analysis have shown that biological neural networks are composed of highly confined objects, which have many neighbors in their local physical neighborhood, which cannot be said for mitochondrial, vascular, and plant root networks:

![image](https://github.com/lukablagoje/physical-networks-spatially-embedded-networks/assets/52599010/cb6c7b62-82d1-43e0-b486-7908f07bafc6)

I further developed the representation to allow the analysis of how each individual neuron is constrained:

![image](https://github.com/lukablagoje/physical-networks-spatially-embedded-networks/assets/52599010/0ee93cf4-eaab-42f5-a69f-9887f246c85b)

In this research, we found that indeed, the physical and network aspects are intertwined.

If you want to read more about this research, please check out the links below:
Publication link: (Coming soon)
arXiv link: https://arxiv.org/abs/2211.13265
# Technical project overview
First, I  access and download the neuron skeleton data, which was further processed to obtain a dataset composed only of points (point clouds) -  **1. obtain_neuron_skeletons.ipynb**.
# Data
To access the original data, you will need to use the Python library provided by the Janelia project, which you can install from https://pypi.org/project/neuprint-python/
If you want to understand this library more, you can read through their documentation https://connectome-neuprint.github.io/neuprint-python/docs/index.html

More specifically, I have saved data for the intermediate steps of my empirical analysis:
**neuron_regions_information** - Used for storing information about the neurons in specific regions.
**neuron_regions_points** - Used for storing point clouds of neuron data.


# Folders
**neuron_regions_information** - Used for storing information about the neurons in specific regions.
**neuron_regions_points** - Used for storing point clouds of neuron data.
