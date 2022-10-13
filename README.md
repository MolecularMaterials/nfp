[![Build Status](https://travis-ci.org/NREL/nfp.svg?branch=master)](https://travis-ci.org/NREL/nfp)
[![PyPI version](https://badge.fury.io/py/nfp.svg)](https://badge.fury.io/py/nfp)

# Neural fingerprint (nfp)

Keras layers for end-to-end learning on molecular structure. Based on Keras, Tensorflow, and RDKit. Source code used in the study [Message-passing neural networks for high-throughput polymer screening](https://arxiv.org/abs/1807.10363)

## Related Work

1. [Convolutional Networks on Graphs for Learning Molecular Fingerprints](https://arxiv.org/abs/1509.09292)
2. [Neural Message Passing for Quantum Chemistry](https://arxiv.org/pdf/1704.01212.pdf)
3. [Relational inductive biases, deep learning, and graph networks](https://arxiv.org/abs/1806.01261)
4. [Neural Message Passing with Edge Updates for Predicting Properties of Molecules and Materials](https://arxiv.org/abs/1806.03146)

## (Main) Requirements

- [rdkit](http://www.rdkit.org/docs/Install.html)
- keras (github master, until [#11548](https://github.com/keras-team/keras/pull/11548) is included in a release)
- tensorflow

## Getting started

This library extends Keras with additional layers for handling molecular structures (i.e., graph-based inputs). There a strong familiarity with Keras is recommended.

An overview of how to build a model is shown in `examples/solubility_test_graph_output.ipynb`. Models can optionally include 3D molecular geometry; a simple example of a network using 3D geometry is found in `examples/model_3d_coordinates.ipynb`.

The current state-of-the-art architecture on QM9 (published in [4]) is included in `examples/schnet_edgeupdate.py`. This script requires qm9 preprocessing to be run before the model is evaluated with `examples/preprocess_qm9.py`.

## Specific instructions for Mo2C-GNN branch
![alt text](https://github.com/MolecularMaterials/nfp/tree/Mo2C-GNN/blob/main/Mo2C-GNN.png?raw=true)
Including code modifications and graph data for:\
***"Accelerating Catalyst Screening via Machine-learned Local Coordination Graph Representations"***\
Hieu A. Doan, Chenyang Li, Logan Ward, Mingxia Zhou, Larry A. Curtiss, and Rajeev S. Assary.

### Preprocess .graphml files and create an input dataframe
- Unzip **notebooks/Oads_Mo2C_graphml.tar.gz** that contains all .graphml files of adsorption geometries:\
`tar -xvf Oads_Mo2C_graphml.tar.gz`
- Run **notebooks/Oads_Mo2C_catalysts_PreprocessGraphStructure.ipynb** to convert graph data into a dataframe input for MPNN model

### Train the model and make predictions on the test set
- Run **notebooks/Oads_Mo2C_catalysts.ipynb**

### Additional instructions for analyzing Atomic Simulation Environment (ASE) database and generating graphs in .graphml format
If you want to experiment with building your own adsorption graphs from the VASP outputs used in the paper, please follow the steps below:
1. Download the data in ASE database format (.db) from the Materials Data Facility [here](https://petreldata.net/mdf/detail/doan_datasets_accelerating_representations_v1.1/)
2. Run **notebooks/DescriptorGen-networkxGraph-pristine-Mo2C.ipynb** and **notebooks/DescriptorGen-networkxGraph-doped-Mo2C.ipynb** with the downloaded "Oads_Mo2C_pristine_MDF.db" and "Oads_Mo2C_doped_MDF.db", respectively
