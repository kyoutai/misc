#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('files', help='input file', nargs="+")
par.add_argument('-o', "--outfile", help='output file')
args = par.parse_args()


fig, ax = plt.subplots()
mark = ["r-", "b-", "g-"]

for idx, name in enumerate(args.files):
    data = np.loadtxt(name)
    label = name.split("SOS")[0]
    ax.plot(data[:, 0], data[:, 1]/np.sum(data[:, 1]), mark[idx], label=label)
ax.set_xlabel('Si-O-Si bond angle (theta)')
ax.set_ylabel("Probability (arb.units)")
ax.set_xlim(75, 180)
plt.legend()
plt.tight_layout()
if args.outfile:
    plt.savefig(args.outfile)
else:
    plt.show()
