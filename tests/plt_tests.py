#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:45:36 2019

@author: simon
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


# Droplet Advection
qd = np.genfromtxt("./test_dadvec.txt")
f, ax = plt.subplots()
ax.plot(qd * 1e3)
ax.grid()
ax.set_xlabel("Time (index)")
ax.set_ylabel("$q$ (g/kg)")
ax.set_title("Droplet Advection with MPDATA")
plt.show()

r = [1.000000E-06, 1.201465E-06, 1.408777E-06, 1.622384E-06, 1.842765E-06, 2.070440E-06,
     2.305966E-06, 2.549943E-06, 2.803018E-06, 3.065886E-06, 3.339296E-06,
     3.624055E-06, 3.921029E-06, 4.231154E-06, 4.555434E-06, 4.894952E-06,
     5.250873E-06, 5.624452E-06, 6.017038E-06, 6.430085E-06, 6.865158E-06,
     7.323941E-06, 7.808247E-06, 8.320027E-06, 8.861382E-06, 9.434573E-06,
     1.004204E-05, 1.068639E-05, 1.137046E-05, 1.209727E-05, 1.287011E-05,
     1.369248E-05, 1.456818E-05, 1.550127E-05, 1.649616E-05, 1.755757E-05,
     1.869058E-05, 1.990066E-05, 2.119372E-05, 2.257609E-05, 2.405461E-05,
     2.563662E-05, 2.733005E-05, 2.914340E-05, 3.108586E-05, 3.316729E-05,
     3.539831E-05, 3.779038E-05, 4.035579E-05, 4.310782E-05, 4.606072E-05,
     4.922985E-05, 5.263176E-05, 5.628423E-05, 6.020644E-05, 6.441899E-05,
     6.894411E-05, 7.380568E-05, 7.902943E-05, 8.464306E-05, 9.067638E-05,
     9.716148E-05, 0.000104, 0.000112, 0.000120, 0.000128,
     0.000138, 0.000148, 0.000158, 0.000170, 0.000183,
     0.000196, 0.000210, 0.000226, 0.000243, 0.000261,
     0.000280, 0.000301, 0.000323, 0.000347, 0.000373,
     0.000401, 0.000431, 0.000463, 0.000498, 0.000535,
     0.000576, 0.000619, 0.000666, 0.000716, 0.000770,
     0.000828, 0.000891, 0.000958, 0.001030, 0.001108,
     0.001192, 0.001283, 0.001380, 0.001485, 0.001597,
     0.001719, 0.001849, 0.001990, 0.002141, 0.002304,
     0.002479, 0.002668, 0.002871, 0.003090, 0.003325,
     0.003578, 0.003851, 0.004145, 0.004461, 0.004801,
     0.005167, 0.005561, 0.005985, 0.006442]
r = np.array(r) * 1e6
K0, K30, K60 = np.genfromtxt("./test_golv.txt")
f, ax = plt.subplots()
ax.semilogx(r, K0 * 1e3)
ax.semilogx(r, K30 * 1e3)
ax.semilogx(r, K60 * 1e3)
ax.set_title("Collision Kernel Test")
ax.grid()
ax.set_xlabel("$\\mu$m (r)")
ax.set_ylabel("g(g/kg)")
ax.xaxis.set_major_formatter(ScalarFormatter())

# Collision
ima = np.genfromtxt("./test_ima.txt")
col_c = np.genfromtxt("./test_colC.txt")

ckern = np.genfromtxt("./test_ckern.txt") * 1e6
XV, YV = np.meshgrid(r, r)
f, ax = plt.subplots()
cb = ax.contourf(XV, YV, ckern)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("r ($\\mu$m)")
ax.set_ylabel("r ($\\mu$m)")
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
f.colorbar(cb)
plt.show()

# Mixing
mix = np.genfromtxt("./test_eddy_mixing.txt")
XV, YV = np.meshgrid(range(mix.shape[0]), range(mix.shape[1]))
f, ax = plt.subplots()
cb = ax.contourf(XV, YV, mix)
ax.set_title("Test Mixing Av. should be 150")
f.colorbar(cb)
plt.show()