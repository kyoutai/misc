#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file', nargs="+")
par.add_argument('-o', '--out', help='output file')
args = par.parse_args()
# 図の体裁
# plt.style.use('classic')
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
plt_dic['axes.grid'] = True
plt_dic['font.size'] = 18
plt_dic['legend.fontsize'] = 18
plt_dic['axes.labelsize'] = 18
plt_dic['xtick.major.size'] = 5
plt_dic['xtick.minor.size'] = 3
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt_dic['savefig.transparent'] = True
# plt_dic['font.sans-serif'] = ['Hiragino Maru Gothic Pro']
plt.rcParams.update(plt_dic)

fig, ax = plt.subplots()
for idx, i in enumerate(args.file):
    data = np.loadtxt(i)
    ax.plot(data[:, 0], data[:, 1], label=i)

    

# fig, ax = plt.subplots(nrows=2, ncols=2)
# ax[0][0].plot(foo[:, 0], foo[:, 1], "r-")
# ax[0][1].plot(foo[:, 0], foo[:, 2], "b-")
# ax[1][0].plot(foo[:, 0], foo[:, 3], "g-")
# ax[1][1].plot(foo[:, 0], foo[:, 4], "-")
# ax[0][0].set_ylim(-2, 2)
# ax[0][1].set_ylim(-2, 2)
# ax[1][0].set_ylim(-2, 2)
# ax[1][1].set_ylim(-0.6, 0.6)
# ax[0][0].set_title("Si-Si")
# ax[0][1].set_title("Si-O")
# ax[1][0].set_title("O-O")
# ax[1][1].set_title("all")
ax.set_xlim(1, 5)
ax.set_ylabel("g(r)")
ax.set_xlabel("Å")
plt.legend()
plt.tight_layout()
if args.out:
    fig.savefig(args.out)
else:
    plt.show()
