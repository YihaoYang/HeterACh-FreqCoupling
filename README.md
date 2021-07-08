# HeterACh-FreqCoupling
Code for publication: "Theta-gamma coupling emerges from spatially heterogeneous cholinergic neuromodulation" in PLoS Computational Biology.

You need Python 3.7 (or above) and MATLAB 2019a (or above) to run the scripts.

The finished data files can be found in the following google-drive link:
1. https://drive.google.com/file/d/17kjkEAFkVWqeOBKRT8flumjyNK6FhMR5/view?usp=sharing
2. https://drive.google.com/file/d/17nH25SdP55Y_lG_dJkbdHfrdXYjZw22W/view?usp=sharing

This directory consists of:

  DisplayResults.m: a MATLAB code with multiple sections corresponding to the generation of figures in the manuscript. Above data files will be needed for plotting.
  
  EI2DNet.m: a MATLAB code that containing the class for neuronal network introduced in the manuscript. It consists of multiple methods including simulation method and network dynamic detection algorithm.
  
  NeuralNet.m: a MATLAB code containing the parent class for EI2DNet.
  
  NeuralNet.py: a Python code containing the method to generate the adjacency matrix needed for simulations.
  
  GenerateAdjMat.py: a Python code to generate the adjacency matrix needed for simulations.
  
  Simulations.m: a MATLAB code with multiple sections corresponding to simulations described in the manuscript.
  
  build2DAdjMatrix.m: a MATLAB code containing methods to generate matrix needed for gks map generation.
  
  showColorMap.m: a MATLAB code containing methods used in DisplayResults.m
  
  adjMat.mat, adjMatForNearest12.mat: datafiles containing the adjacency matrix needed for simulations, which can be generated using GenerateAdjMat.py

