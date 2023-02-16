#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file', nargs="+")
par.add_argument('-o', '--out', help='output file')
args = par.parse_args()


def theta(q):
    return np.arcsin(1.5406/4/np.pi*q)


fig, ax = plt.subplots(nrows=2, ncols=2)

for idx, i in enumerate(args.file):
    foo = np.loadtxt(i, dtype=float)
    # ax[0][0].plot(np.degrees(theta(foo[:, 0]))/2, foo[:, 1]+idx, "r-")
    # ax[0][1].plot(np.degrees(theta(foo[:, 0]))/2, foo[:, 2]+idx, "b-")
    # ax[1][0].plot(np.degrees(theta(foo[:, 0]))/2, foo[:, 3]+idx, "g-")
    # ax[1][1].plot(np.degrees(theta(foo[:, 0]))/2, foo[:, 4]+idx, "-")
    ax[0][0].plot(foo[:, 0], foo[:, 1]+idx, "r-")
    ax[0][1].plot(foo[:, 0], foo[:, 2]+idx, "b-")
    ax[1][0].plot(foo[:, 0], foo[:, 3]+idx, "g-")
    ax[1][1].plot(foo[:, 0], foo[:, 4]+idx, "-")
# ax[0][0].set_xlim(0, 15)
# ax[0][1].set_xlim(0, 15)
# ax[1][0].set_xlim(0, 15)
# ax[1][1].set_xlim(0, 15)
ax[0][0].set_xlim(0, 60)
ax[0][1].set_xlim(0, 60)
ax[1][0].set_xlim(0, 60)
ax[1][1].set_xlim(0, 60)
ax[0][0].set_ylim(-2, 2+idx)
ax[0][1].set_ylim(-2, 2+idx)
ax[1][0].set_ylim(-2, 2+idx)
ax[1][1].set_ylim(-0.6, 0.6+idx)
ax[0][0].set_title("Si-Si")
ax[0][1].set_title("Si-O")
ax[1][0].set_title("O-O")
ax[1][1].set_title("all")
# ax.set_xlim(0, 12)
# ax.set_ylabel('y_axis')
plt.tight_layout()
if args.out:
    fig.savefig(args.out)
else:
    plt.show()
