#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input SQ file', nargs="+")
par.add_argument('-o', '--out', help='output file')
# par.add_argument('-c', '--csv', help='additional csv file')
args = par.parse_args()

# plt.rcParams.update({'legend.fontsize': 12,
#                      'axis.fontsize': 12})

fig, ax = plt.subplots(nrows=2, ncols=2)
color = ["black", "blue", "red", "purple", "orange", "brown"]

for idx, fr in enumerate(args.file):
    foo = np.loadtxt(fr, dtype=float)
    label = fr.split(".")[0]
    ax[0][0].plot(foo[:, 0], foo[:, 1]+idx, "-", color=color[idx])
    ax[0][1].plot(foo[:, 0], foo[:, 2]+idx, "-", color=color[idx])
    ax[1][0].plot(foo[:, 0], foo[:, 3]+idx, "-", color=color[idx])
    ax[1][1].plot(foo[:, 0], foo[:, 4]+idx, "-", color=color[idx], label=label)
    # ax[0][0].set_ylim(-2, 2)
    # ax[0][1].set_ylim(-2, 2)
    # ax[1][0].set_ylim(-2, 2)
    # ax[1][1].set_ylim(-0.6, 0.6)
    ax[0][0].set_title("Si-Si")
    ax[0][1].set_title("Si-O")
    ax[1][0].set_title("O-O")
    ax[1][1].set_title("all")
    ax[0][0].set_xlim(0, 20)
    ax[0][1].set_xlim(0, 20)
    ax[1][0].set_xlim(0, 20)
    ax[1][1].set_xlim(0, 20)
    ax[0][0].set_ylim(-10, 10)
    ax[0][1].set_ylim(-5, 8)
    ax[1][0].set_ylim(-5, 8)
    ax[1][1].set_ylim(-5, 8)
    ax[0][0].set_xlabel("Q")
    ax[0][1].set_xlabel("Q")
    ax[1][0].set_xlabel("Q")
    ax[1][1].set_xlabel("Q")
    ax[0][0].set_xlabel("S(Q)")
    ax[0][1].set_xlabel("S(Q)")
    ax[1][0].set_xlabel("S(Q)")
    ax[1][1].set_xlabel("S(Q)")
    # ax.set_ylabel('y_axis')
    ax[1][1].legend()
plt.tight_layout()
if args.out:
    fig.savefig(args.out)
else:
    plt.show()
