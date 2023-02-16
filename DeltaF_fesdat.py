#!/usr/bin/env python
import numpy as np
import argparse
text = """
"""


class HelpFormat(argparse.HelpFormatter):
    # https://qiita.com/moshi/items/f354a2e24424244c0451
    # 'max_help_position'のデフォルト値を「24」から「36」へ変更
    def __init__(self, prog,
                 indent_increment=2, max_help_position=36, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)


par = argparse.ArgumentParser(description=text, formatter_class=HelpFormat)
par.add_argument("file", help="input Free Energy Surface data")
par.add_argument('-t', '--temp', type=int, default=2700,
                 help="系の温度")
par.add_argument('--threshold', type=int, default=150,
                 help="2相の区別境界")
par.add_argument('--unit', choices=["kj", "kcal"], default="kj",
                 help="FES エネルギー単位を指定")
par.add_argument('-o', '--outfile', help="outfile name")
args = par.parse_args()

# set constant
cal_to_kj = 4.184
kB = 1.380649e-23               # J/K
mol = 6.0221407e23
kBT = kB * mol * args.temp * 1e-3  # kj/mol

buff = np.loadtxt(args.file)
cv, fe = buff[:, 0], buff[:, 1]
fe = np.exp((-fe+np.min(fe)) / kBT)
mask_g = np.where((0 < cv) & (cv < args.threshold), True, False)
mask_c = np.where((args.threshold < cv) & (cv < 300), True, False)
fglass = np.sum(fe * mask_g)
fcryst = np.sum(fe * mask_c)
deltaF = -kBT * np.log(fglass/fcryst)
print(fglass/fcryst)
print(kBT)
exit(deltaF)
