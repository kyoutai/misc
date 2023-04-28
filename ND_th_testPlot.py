#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import re
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
plt_dic['savefig.dpi'] = 300
plt_dic['savefig.transparent'] = True
# plt_dic['font.sans-serif'] = ['Hiragino Maru Gothic Pro']
plt.rcParams.update(plt_dic)

now = os.getcwd() + "/"
dir_list = os.listdir(now)
frdic, valdic = {}, {}
for name in dir_list:
    nlist = name.split(".")
    if nlist[-1] == "out":
        key = nlist[0][2]
        if key not in frdic.keys():
            frdic[key] = [int(nlist[0][5:7])]
            valdic[key] = []
        else:
            frdic[key].append(int(nlist[0][5:7]))
        fr = open("ND{}th{:02d}test.out".format(key, frdic[key][-1])).read()
        tmp = 0
        for v in re.findall("Perform.+?rs/ns, (\d+\.\d+) time", fr, re.DOTALL):
            tmp += float(v)
        valdic[key].append(tmp/3)

fig, ax = plt.subplots()
for key in frdic.keys():
    label = str(key) + "_node"
    ax.plot(frdic[key], valdic[key], "o-", markersize=5, label=label)
# print(array)
ax.set_xlabel('number of threads')
ax.set_ylabel('timesteps/s')
plt.legend()
# plt.tight_layout()
plt.show()
