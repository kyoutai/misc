#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input SQ file', nargs="+")
par.add_argument('-o', '--out', help='output file')
par.add_argument('-c', '--csv', help='additional csv file', nargs="*")
args = par.parse_args()

# plt.rcParams.update({'legend.fontsize': 12,
#                      'axis.fontsize': 12})

fig, ax = plt.subplots(nrows=2, figsize=(6,6))
cnt = 0
co = ["black", "blue", "red", "purple", "orange", "brown"]

for idx, fr in enumerate(args.file):
    foo = np.loadtxt(fr, dtype=float)
    label = fr.split(".")[0]
    ax[0].plot(foo[:, 0], foo[:, 4]+cnt, "-", color=co[cnt], label=label)
    ax[1].plot(foo[:, 0], foo[:, 5]+cnt, "-", color=co[cnt], label=label)
    cnt += 1
if args.csv:
    for fr in args.csv:
        lab = "experiment"
        foo = np.loadtxt(fr, dtype=float, skiprows=2, delimiter=',')
        ax[0].plot(foo[:, 2], foo[:, 3]+cnt, "-", color=co[cnt], label=lab)
        ax[1].plot(foo[:, 0], foo[:, 1]+cnt, "-", color=co[cnt], label=lab)
        cnt += 1
ax[0].set_xlim(0, 15); ax[1].set_xlim(0, 15)
ax[0].set_ylim(-1, 6); ax[1].set_ylim(-1, 6)
ax[0].set_xlabel("Q"); ax[1].set_xlabel("Q")
ax[0].set_xlabel("S(Q)"); ax[1].set_xlabel("S(Q)")
ax[0].set_title("S(Q)-XRD"); ax[1].set_title("S(Q)-ND")
ax[0].legend(loc="best"); ax[1].legend(loc="best")
plt.tight_layout()
if args.out:
    fig.savefig(args.out)
else:
    plt.show()
