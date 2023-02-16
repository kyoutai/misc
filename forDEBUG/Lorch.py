#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

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
space = np.linspace(0, 100, 101)/100
hoge = np.sin(space)/space
ax.plot(space, hoge, "o-")
ax.set_xlabel('x')
ax.set_ylabel('sin(x)')
fig.savefig("lorch.png")
# plt.tight_layout()
# plt.show()
