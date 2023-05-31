#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

plt.rcParams.update({'legend.fontsize': 12,
                     'axes.labelsize': 20,
                     'axes.titlesize': 20,
                     "figure.figsize": [7, 3],
                     'xtick.labelsize': 16,
                     'ytick.labelsize': 16})


par = argparse.ArgumentParser(description="bar")
par.add_argument('files', help='input file', nargs="+")
par.add_argument('-o', '--outfile', help="outfile name")
par.add_argument('--switch', help="ステップ数削減の無効化", type=int)
args = par.parse_args()

fig, ax = plt.subplots(figsize=(14, 6*len(args.files)), nrows=len(args.files))

for idx, name in enumerate(args.files):
    fr = np.loadtxt(name)
    ffname = name.split("/")
    if args.switch:
        fr = fr[(fr[:, 0] % args.switch == 0)]  # 要素数を削減して軽量化
    if fr.shape[1] == 2:  # 1次元集団変数プロット
        ax[idx].plot(fr[:, 0]/10e5, fr[:, 1], ".", markersize=1)
        ax[idx].set_xlabel('steps (nano sec.)')
        ax[idx].set_ylabel('CV (XRD peak)')
        ax[idx].plot([np.min(fr[:, 0])/10e5, np.max(fr[:, 0])/10e5], [130, 130], "--r")
        # ax[idx].plot([np.min(fr[:, 0])/10e5, np.max(fr[:, 0])/10e5], [25, 25], "--r")
    elif fr.shape[1] == 3:  # 2次元集団変数プロット
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
        # ax[idx].set_xlim(-1, 2)
        fig.colorbar(mappable, ax=ax3, label="time (ps)")
        plt.tight_layout()
        print(np.min(fr[:, 1]), np.max(fr[:, 1]))
        print(np.min(fr[:, 2]), np.max(fr[:, 2]))
    else:  # 多次元集団変数の場合。
        # 例えば、Well-Tempered で再重み付けに必要なCVを含むとき。
        # ex) XRD, VOLUME, ENERGY, etc...
        fig, ax = plt.subplots()
        ax[idx].plot(fr[:, 0]/10e5, fr[:, -1], ".", markersize=1)
        ax[idx].set_xlabel('steps (nano sec.)')
        ax[idx].set_ylabel('CV (XRD peak)')
        ax[idx].plot([np.min(fr[:, 0])/10e5, np.max(fr[:, 0])/10e5], [130, 130], "--r")
        ax[idx].plot([np.min(fr[:, 0])/10e5, np.max(fr[:, 0])/10e5], [25, 25], "--r")
    if len(ffname) == 2:
        ax[idx].set_title(ffname[0])

plt.tight_layout()
if args.outfile:
    fig.savefig(args.outfile, dpi=300)
else:
    plt.show()
