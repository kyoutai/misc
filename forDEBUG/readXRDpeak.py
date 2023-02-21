#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file', nargs="+")
par.add_argument('-o', '--outfile', help='output file')
args = par.parse_args()

# 図の体裁
# plt.style.use('classic')
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
plt_dic['axes.grid'] = True
plt_dic['font.size'] = 12
plt_dic['legend.fontsize'] = 12
plt_dic['axes.labelsize'] = 14
plt_dic['xtick.major.size'] = 5
plt_dic['xtick.minor.size'] = 3
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt_dic['savefig.transparent'] = True
# plt_dic['font.sans-serif'] = ['Hiragino Maru Gothic Pro']
plt.rcParams.update(plt_dic)

fig, ax = plt.subplots()


x = np.linspace(15, 50, 36)
for idx, i in enumerate(args.file):
    data = np.loadtxt(i)
    data = data[-5:, :]
    y = np.average(data[:, 1:], axis=0)
    ax.plot(x, y, ".")
# ax.set_xlabel('x_axis')
# ax.set_ylabel('y_axis')
plt.tight_layout()
if args.outfile:
    fig.savefig(args.outfile)
else:
    plt.show()
