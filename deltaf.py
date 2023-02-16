#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt
# 図の体裁
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
plt_dic['font.size'] = 12
plt_dic['legend.fontsize'] = 12
plt_dic['axes.labelsize'] = 14
plt_dic['xtick.major.size'] = 5
plt_dic['xtick.minor.size'] = 3
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt_dic['savefig.transparent'] = True
plt.rcParams.update(plt_dic)

text = """
"""


class HelpFormat(argparse.HelpFormatter):
    # https://qiita.com/moshi/items/f354a2e24424244c0451
    # 'max_help_position'のデフォルト値を「24」から「36」へ変更
    def __init__(self, prog,
                 indent_increment=2, max_help_position=36, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)


par = argparse.ArgumentParser(description=text, formatter_class=HelpFormat)
par.add_argument('end', help='読み込む fes ファイルの上限番号', type=int)
par.add_argument('time', help='fes_.dat ファイル間の時間間隔', type=int)
par.add_argument('-t', '--temp', type=int, default=2700,
                 help="系の温度")
par.add_argument('--threshold', type=int, default=150,
                 help="2相の区別境界")
par.add_argument('--unit', choices=["kj", "kcal"], default="kj",
                 help="FES エネルギー単位を指定")
par.add_argument('-o', '--outfile', help="outfile name")
args = par.parse_args()
fig, ax = plt.subplots()
cm = plt.get_cmap("gnuplot")

# set constant
cal_to_kj = 4.184
kB = 1.380649e-23               # J/K
mol = 6.0221407e23
kBT = kB * mol * args.temp * 1e-3  # kj/mol
cnt, cum, mini, maxi = 0, 0, 0, 300
if args.unit == "kj":
    coeff = 4.184
else:
    kBT /= cal_to_kj            # kcal/mol
    coeff = 1.0

time = np.linspace(0, args.end, args.end + 1) * args.time
deltaF = np.empty((args.end + 1))

for i in range(0, args.end + 1):
    fglass, fcryst = 0.0, 0.0
    buff = np.loadtxt("fes_{}.dat".format(i))
    cv, fe = buff[:, 0], buff[:, 1]
    minf = np.min(fe)
    fe = np.exp(-(fe+minf) / kBT)
    mask_g = np.where((0 < cv) & (cv < args.threshold), True, False)
    mask_c = np.where((args.threshold < cv) & (cv < 300), True, False)
    fglass = np.sum(fe * mask_g)
    fcryst = np.sum(fe * mask_c)
    dF = -kBT * np.log(fcryst/fglass)
    deltaF[i] = dF
ax.plot(time, deltaF)
ax.set_xlabel("time (ps)")
ax.set_ylabel("deltaF ("+args.unit+"/mol)", fontsize=14)
# ax.set_xlim(0, 300)
# ax.set_ylim(0, 1000)

if args.outfile:
    fig.savefig(args.outfile, dpi=400)
else:
    plt.show()
