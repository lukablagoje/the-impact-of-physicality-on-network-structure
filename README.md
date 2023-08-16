# Research overview
This project is a part of the published research work "Understanding the Impact of physicality on network structure", done in the emerging field of Physical Networks, which aims to understand the properties of three-dimensional networked systems (such as biological neural networks).

A concept of "generalized meta-graph", which I implemented in Python (using point clouds and kd-trees), by encoding information on Euclidean distances between neighboring edges.

I applied this new representation to efficiently solve a computationally challenging task: finding out how many unique neuron-to-neuron physical collisions are obtained if the thickness of the neurons is increased for 20 different additive factors. The results of this analysis have shown that biological neural networks are composed of highly confined objects, which have many neighbors in their local physical neighborhood, which cannot be said for mitochondrial, vascular, and plant root networks.

Finally, I further developed the representation to analyze not the entire dataset, but each individual neuron and its physical confinement.

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
