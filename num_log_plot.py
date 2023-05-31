#!/usr/bin/env python
import argparse
import numpy as np
import re
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
plt_dic['savefig.dpi'] = 300
plt_dic['savefig.transparent'] = True
# plt_dic['font.sans-serif'] = ['Hiragino Maru Gothic Pro']
plt.rcParams.update(plt_dic)


description = """何もプロットされないならば、正規表現に間違いあるのでプログラムを書き換えましょう。"""
par = argparse.ArgumentParser(description=description)
par.add_argument("ev", help="irradiation energy", type=int)
par.add_argument('files', nargs='+')
par.add_argument('-n', "--natoms", help="number of atoms", type=int)
par.add_argument('-o', "--out", help="outfile name (hoge.png)")
args = par.parse_args()


def Density(files, ev1st):
    for idx, i in enumerate(files):
        body = open(i).read()
        data = re.findall(
           "Step           Dt.+?c_2      \n(.+?)\nLoop time", body, re.DOTALL)
        # data = re.findall(
        #     "Step Dt.+?c_2 \n(.+?)\nLoop time", body, re.DOTALL)
        data = np.array(" ".join(data).split()).reshape((-1, 11)
                                                        ).astype(float)
        enthalpy = np.mean(data[-50:, 8])
        density = np.mean(data[-50:, 9])
        x = idx * ev1st
        bx.plot(x, enthalpy*conv/args.natoms, "r.")
        ax.plot(x, density, "b.")


conv = 1/23.060  # 1 eV = 23.060 kcal/mol
fig, ax = plt.subplots()
bx = ax.twinx()
Density(args.files, args.ev / args.natoms)
bx.set_ylabel('${\\rm Enthalpy\ eV/atom}$')
ax.set_xlabel('${\\rm Deposit\ energy\ eV/atom}$')
ax.set_ylabel('${\\rm Density\ g/cm^3}$')
ax.set_ylim(2.0, 2.7)
ax.set_title(args.files[0].split("_")[0]+"-"+str(args.ev)+" ev/time")
plt.tight_layout()
if args.out:
    plt.savefig(args.out)
    print(args.out, "was created.")
else:
    plt.show()
