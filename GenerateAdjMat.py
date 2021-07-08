# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 16:18:22 2021

@author: yihao yang
"""

# scripts to generate adjacency matrix

from NeuroNet import NeuroNet
import numpy as np
import scipy.io

net = NeuroNet()
net.mexicanHat()
scipy.io.savemat("adjMat.mat", mdict={"adjMat":net.adjMat})

net.mexicanHat(degreeEE=12,degreeEI=12,weightEE=1,weightEI=1)
scipy.io.savemat("adjMatForNearest12.mat", mdict={"adjMatForNearest12":net.adjMat + np.eye(500)})
