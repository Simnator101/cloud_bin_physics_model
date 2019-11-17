#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 00:06:46 2019

@author: simon
"""

import numpy as np
import matplotlib.pyplot as  plt

u,v = np.genfromtxt("../data/test_uv.txt")
x = np.arange(u.size)
y = np.arange(v.size)

XV, YV = np.meshgrid(x,y)
u = u.reshape((101, 101))
v = v.reshape((101, 101))

f, ax = plt.subplots()
ax.quiver(XV, YV, u, v)
plt.show()