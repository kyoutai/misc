#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
mini, maxi = 0, 300
m, M = 150, 300
cal_to_kj = 4.184

par = argparse.ArgumentParser(description="hills_plot1.py fes.dat")
par.add_argument('file', help='one step FESdat file')
par.add_argument('-m', '--mode', choices=["kj", "kcal"], default="kj",
                 help="choose y axis unit")
args = par.parse_args()

fig, ax = plt.subplots()

fr = np.loadtxt(args.file)
mask1 = mini < fr[:, 0]
mask2 = fr[:, 0] < maxi
mask = mask1 * mask2
# for idx, i in enumerate(fr[:, 0]):
#     print(i, mask[idx])
mask1 = m < fr[:, 0]
mask2 = fr[:, 0] < M
minimum = np.min(fr[mask1*mask2, 1])
print(minimum)
X = fr[mask, 0]
Y = fr[mask, 1] - minimum
if args.mode == "kj":
    Y *= cal_to_kj

ax.plot(X, Y)
ax.set_xlim(mini, maxi)
# ax.set_ylim(0, np.max(Y[:-5]))
ax.set_ylim(0, 800)
plt.xlabel('CV', fontsize=14)
plt.ylabel(args.mode+'/mol', fontsize=14)
plt.show()
