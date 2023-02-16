#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file', nargs="*")
args = par.parse_args()

fig, ax = plt.subplots()
for i in args.file:
    fr = np.loadtxt(i)
    ax.plot(fr[:, 0], fr[:, 1]-np.min(fr[:, 1]), markersize=1,
            label=i.split(".")[0])
    plt.xlabel('CV', fontsize=14)
    plt.ylabel('kcal/mol', fontsize=14)

plt.legend()
plt.show()
