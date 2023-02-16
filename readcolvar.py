#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

#print(plt.rcParams.keys())

plt.rcParams.update({'legend.fontsize': 12,
                     'axes.labelsize': 20,
                     "figure.figsize": [7, 3],
                     'xtick.labelsize': 16,
                     'ytick.labelsize': 16})


par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-o', '--outfile', help="outfile name")
par.add_argument('--switch', help="ステップ数削減の無効化", action="store_true")
args = par.parse_args()

# with open(args.file) as f:
#     line = f.readline().strip().split()[2:]
#     print(line)


fr = np.loadtxt(args.file)
if args.switch:
    fr = fr[(fr[:, 0] % 50000 == 0)]  # 要素数を削減して軽量化
if fr.shape[1] == 2:
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(fr[:, 0]/10e8, fr[:, 1], ".", markersize=1)
    ax.set_xlabel('steps (micro sec.)')
    ax.set_ylabel('CV (XRD Si-Si{111})')
    ax.plot([np.min(fr[:, 0])/10e8, np.max(fr[:, 0])/10e8], [150, 150], "--r")
    # ax.set_xlim(-1, 2)
    name = args.file.split("/")
    if len(name) == 2:
        ax.set_title(name[0])
elif fr.shape[1] == 3:
    fig = plt.figure(figsize=(7, 3))
    ax1 = plt.subplot2grid((2, 2), (0, 0))
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    ax3 = plt.subplot2grid((2, 2), (0, 1), rowspan=2)
    cm = plt.cm.get_cmap("gnuplot")
    ax1.plot(fr[:, 0]/1e6, fr[:, 1])
    ax2.plot(fr[:, 0]/1e6, fr[:, 2])
    mappable = ax3.scatter(fr[:, 1], fr[:, 2], c=fr[:, 0]/1e6)
    ax3.set_xlabel('CV 1')
    ax3.set_ylabel('CV 2')
    ax1.set_ylabel('CV 1')
    ax2.set_ylabel('CV 2')
    ax1.set_xlabel('time (ps)')
    ax2.set_xlabel('time (ps)')
    # ax.set_xlim(-1, 2)
    fig.colorbar(mappable, ax=ax3, label="time (ps)")
    plt.tight_layout()
    print(np.min(fr[:, 1]), np.max(fr[:, 1]))
    print(np.min(fr[:, 2]), np.max(fr[:, 2]))
else:
    fig, ax = plt.subplots()
    ax.plot(fr[:, 0]/10e8, fr[:, -1], ".", markersize=1)
    ax.set_xlabel('steps (micro sec.)')
    ax.set_ylabel('CV (XRD Si-Si{111})')
    ax.plot([np.min(fr[:, 0])/10e8, np.max(fr[:, 0])/10e8], [150, 150], "--r")

plt.tight_layout()
if args.outfile:
    fig.savefig(args.outfile, dpi=300)
else:
    plt.show()
