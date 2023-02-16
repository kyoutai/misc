#!/usr/bin/env python
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib.colors as colors
# 図の体裁
# plt.style.use('classic')
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
# plt_dic['axes.grid'] = True
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

# set constant
mini, maxi = 0, 300
cal_to_kj = 4.184
text = """
メタダイナミクス計算の出力ファイル HILLS を後処理した fes_*.dat ファイルを読み、自由エネルギー表面 FES を描画するプログラム。後処理コードは、 plumed sum_hills --minto_zero --stride 31250 を想定し、1ファイル毎に10 ns ずつ進むことを考慮してコーディングした。args.stride*10 ns 毎の FES と1 μs以後の平均 FES をそれぞれプロットした図を出力
"""


class HelpFormat(argparse.HelpFormatter):
    # https://qiita.com/moshi/items/f354a2e24424244c0451
    # 'max_help_position'のデフォルト値を「24」から「36」へ変更
    def __init__(self, prog,
                 indent_increment=2, max_help_position=36, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)


par = argparse.ArgumentParser(description=text, formatter_class=HelpFormat)
par.add_argument('end', help='読み込む fes ファイルの上限番号', type=int)
par.add_argument('-s', '--stride', type=int, default=10,
                 help="FES のプロット頻度")
par.add_argument('-a', '--average', type=int, default=100,
                 help="FES のプロット頻度")
par.add_argument('-m', '--maximum', type=int, default=1000,
                 help="FES の縦軸(エネルギー)上限")
par.add_argument('--unit', choices=["kj", "kcal"], default="kj",
                 help="FES エネルギー単位を指定")
par.add_argument('-o', '--outfile', help="outfile name")
par.add_argument("--mintozero", action="store_true")
args = par.parse_args()
fig, ax = plt.subplots(ncols=2, figsize=(12, 4))
cm = plt.get_cmap("gnuplot")
num_data = round((args.end) / args.stride)


cnt, cum, mini, maxi = 0, 0, 0, 300
hilmax, hilmin = 1000, 0
if args.unit == "kj":
    coeff = 4.184
else:
    coeff = 1.0

buff = np.loadtxt("fes_{}.dat".format(0))
mask1 = mini < buff[:, 0]
mask2 = buff[:, 0] < maxi
mask = mask1 * mask2
X = buff[mask, 0]
Y = np.empty_like(X)

for i in range(0, args.end, args.stride):
    # print(i, end=" ")
    cnt += 1
    ax[0].plot(X, np.loadtxt("fes_{}.dat".format(i))[mask, 1]*coeff,
               color=cm(cnt/num_data))
    if i >= args.average:
        cum += 1
        Y += np.loadtxt("fes_{}.dat".format(i))[mask, 1]
        ax[1].plot(X, Y*coeff/cum, color=cm(cnt/num_data))
ax[0].set_ylabel(args.unit+"/mol", fontsize=14)
ax[1].set_ylabel(args.unit+"/mol", fontsize=14)
ax[0].set_title("Cumulative FES", fontsize=16)
ax[1].set_title("Cumulative average FES", fontsize=16)
hilmin = 0
if not args.mintozero:
    hilmin = np.min(np.loadtxt("fes_{}.dat".format(args.end))[mask, 1]*coeff)
    print(hilmin)
    hilmax = 0
ax[0].set_xlim(0, 300)
ax[1].set_xlim(0, 300)
ax[0].set_ylim(0, args.maximum)
ax[1].set_ylim(-30, args.maximum)
# カラーバー凡例をプロット
# https://qiita.com/tmoooomooooo/items/fa60b4c726d7b1feffb5
axpos = ax[1].get_position()
cbar_ax = fig.add_axes([0.87, axpos.y0, 0.02, axpos.height])
norm = colors.Normalize(vmin=0, vmax=args.end)
mappable = ScalarMappable(cmap='gnuplot', norm=norm)
mappable._A = []
fig.colorbar(mappable, cax=cbar_ax, label="time(10 ps)")
# カラーバー位置を調整
plt.subplots_adjust(right=0.85)
plt.subplots_adjust(wspace=0.1)
# 軸ラベルを消す
# https://www.delftstack.com/ja/howto/matplotlib/how-to-hide-axis-text-ticks-and-or-tick-labels-in-matplotlib/
ax[1].axes.yaxis.set_visible(False)

if args.outfile:
    fig.savefig(args.outfile, dpi=400)
else:
    plt.show()
