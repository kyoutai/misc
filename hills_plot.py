#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# plt.rcParams.update({'legend.fontsize': 12,
#                      'axis.labelsize': 12})
# plt.rcParams.update({'legend.fontsize': 12})
fig, ax = plt.subplots()
par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-o', '--outfile', help="outfile name")
par.add_argument('temp', help='system temperature', type=int)
args = par.parse_args()
# 図の体裁
# plt.style.use('classic')
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
# plt_dic['axes.grid'] = True
plt_dic['font.size'] = 10
# plt_dic['legend.fontsize'] = 12
plt_dic['axes.labelsize'] = 14
plt_dic['xtick.major.size'] = 5
plt_dic['xtick.minor.size'] = 3
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt_dic['savefig.transparent'] = True
# plt_dic['font.sans-serif'] = ['Hiragino Maru Gothic Pro']
plt.rcParams.update(plt_dic)

# # s_size = int(tupp[0])
# array = np.array([23, 45, 67, 78])
# tupp = np.unique(array)
# print(tupp)
# print(tupp.shape)
# print(type(tupp.shape[0]))
# # print(s_size)
# exit()
cal_to_kj = 4.184
kB = 1.380649e-23               # J/K
mol = 6.0221407e23
kBT = kB * mol * args.temp * 1e-3  # kj/mol

fr = np.loadtxt(args.file)
s_unique, t_unique = np.unique(fr[:, 0]), np.unique(fr[:, 1])
s_size, t_size = s_unique.shape[0], t_unique.shape[0]
XX, YY = np.meshgrid(s_unique, t_unique, indexing='xy')
fr = fr.reshape(t_size, s_size, 5)
print("resize array into {} x {} x 5".format(t_size, s_size))
print(fr.shape)
fr[:, :, 2] = fr[:, :, 2] / kBT
# print(kBT)
# exit()

mapp = ax.pcolormesh(XX, YY, fr[:, :, 2], cmap='gnuplot2',
                     norm=Normalize(vmin=0, vmax=15))
pp = fig.colorbar(mapp, ax=ax, orientation="vertical")
ax.set_xlim(-10, 230)
ax.set_ylim(20, 90)
plt.text(0, 85, str(args.temp)+" K", fontsize=14)
plt.xlabel('{111}', fontsize=14)
plt.ylabel('{022}', fontsize=14)
if args.outfile:
    fig.savefig(args.outfile, dpi=300)
else:
    plt.show()
