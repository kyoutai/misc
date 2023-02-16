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
メタダイナミクス計算の出力ファイル HILLS を後処理した fes_*.dat ファイルを読み、自由エネルギー表面 FES を描画するプログラム。後処理コードは、 plumed sum_hills --minto_zero --stride 10000 を想定し、1ファイル毎に10 ps ずつ進むことを考慮してコーディングした。args.stride*10 ps 毎の FES と1 μs以後の平均 FES をそれぞれプロットした図を出力
"""


class HelpFormat(argparse.HelpFormatter):
    # https://qiita.com/moshi/items/f354a2e24424244c0451
    # 'max_help_position'のデフォルト値を「24」から「30」へ変更
    def __init__(self, prog,
                 indent_increment=2, max_help_position=30, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)


par = argparse.ArgumentParser(description=text, formatter_class=HelpFormat)
par.add_argument('file', help='input file name')
args = par.parse_args()
fig, ax = plt.subplots()
cm = plt.get_cmap("gnuplot")

# ------------------------------
# ここにプログラムを作成       |
# ------------------------------

# カラーバー凡例をプロット
# https://qiita.com/tmoooomooooo/items/fa60b4c726d7b1feffb5
axpos = ax.get_position()
cbar_ax = fig.add_axes([0.87, axpos.y0, 0.02, axpos.height])
norm = colors.Normalize(vmin=0, vmax=100)
mappable = ScalarMappable(cmap='gnuplot', norm=norm)
mappable._A = []
fig.colorbar(mappable, cax=cbar_ax, label="time")
# カラーバー位置を調整
plt.subplots_adjust(right=0.85)
plt.subplots_adjust(wspace=0.1)
# 軸ラベルを消す
# https://www.delftstack.com/ja/howto/matplotlib/how-to-hide-axis-text-ticks-and-or-tick-labels-in-matplotlib/
ax.axes.yaxis.set_visible(False)

# output figurep
if args.outfile:
    fig.savefig(args.outfile, dpi=400)
else:
    plt.show()
