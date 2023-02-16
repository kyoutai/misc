#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

par = argparse.ArgumentParser(description="bar")
par.add_argument('reffile', help='input reference SQ file')
par.add_argument('myfile', help='input my SQ file', nargs="+")
par.add_argument('-o', '--out', help='output file')
# par.add_argument('-c', '--csv', help='additional csv file')
args = par.parse_args()

# plt.rcParams.update({'legend.fontsize': 12,
#                      'axis.fontsize': 12})

fig, ax = plt.subplots(nrows=2, figsize=(8, 6))
# fig, ax = plt.subplots()
color = ["black", "blue", "purple", "orange", "brown"]

for idx, i in enumerate(args.myfile):
    f = np.loadtxt(i, dtype=float)
    label = i.split(".")[-2]
    ax[0].plot(f[:, 0], f[:, 4]+1, "-", color=color[idx], label=label)
    ax[1].plot(f[:, 0], f[:, 5]+1, "-", color=color[idx], label=label)

label = "experiment"
if args.reffile.split(".")[-1] == "csv":
    f = np.loadtxt(args.reffile, dtype=float, skiprows=2, delimiter=',')
    ax[0].plot(f[:, 0], f[:, 1], ".", color="red", label=label, markersize=1)
    ax[1].plot(f[:, 2], f[:, 3], ".", color="red", label=label, markersize=1)
else:
    f = np.loadtxt(args.reffile, dtype=float)
    ax[0].plot(f[:, 0], f[:, 4], ".", color="red", label=label, markersize=1)
    ax[1].plot(f[:, 0], f[:, 5], ".", color="red", label=label, markersize=1)

ax[0].set_xlim(0, 15); ax[1].set_xlim(0, 15)
# ax[0].set_ylim(0, 3); ax[1].set_ylim(0, 3)
ax[0].set_xlabel("Q"); ax[1].set_xlabel("Q")
ax[0].set_xlabel("S(Q)"); ax[1].set_xlabel("S(Q)")
ax[0].set_title("S(Q)-XRD"); ax[1].set_title("S(Q)-ND")
ax[0].legend(loc="best"); ax[1].legend(loc="best")
plt.tight_layout()
if args.out:
    fig.savefig(args.out)
else:
    plt.show()
