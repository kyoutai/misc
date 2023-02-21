#!/usr/bin/env python
import numpy as np
import argparse
# www.pnas.org/cgi/doi/10.1073/pnas.1803919115
# cv2 metadynamics inputfile maker

# description = "cv値を2つのXRDピークとし、Si,O両方を考慮するPLUMED実行ファイルを生成するプログラム。本プログラムは、18原子のlammpsdataファイルを使用することを想定して作成した。"
description = "make2dimDat.py 324"
par = argparse.ArgumentParser(description=description)
par.add_argument('atoms', help='number of atoms', default="192", type=int)
par.add_argument("--alpha", action='store_true',
                 help="alpha quartz 専用オプション。group を変えます")
par.add_argument("--lorch", action='store_true', help="use lorch function")

args = par.parse_args()
theta = np.linspace(15, 50, 36)


def AtomFactorS(ang):            # ang: 2theta
    lbd = 1.5406
    q1 = (np.sin(ang/360*np.pi)/lbd)**2
    fs1 = (5.275329*np.exp(-2.631338*q1) + 3.191038*np.exp(-33.73073*q1) +
           1.511514*np.exp(-0.081119*q1) + 1.356849*np.exp(-86.28864*q1) +
           2.519114*np.exp(-1.170087*q1) + 0.145073)
    return fs1


def AtomFactorO(ang):            # ang: 2theta
    lbd = 1.5406
    q = (np.sin(ang/360*np.pi)/lbd)**2
    fs = (2.960427*np.exp(-14.18226*q) + 2.508818*np.exp(-5.936858*q) +
          0.637853*np.exp(-0.112726*q) + 0.722838*np.exp(-34.95848*q) +
          1.142756*np.exp(-0.390240*q) + 0.027014)
    return fs


Q, fss = np.empty_like(theta), np.empty_like(theta)
fso, foo = np.empty_like(theta), np.empty_like(theta)
for idx, i in enumerate(theta):
    Q[idx] = 4*np.pi/1.5406 * np.sin((i/2)/180*np.pi)
    otmp = AtomFactorO(i)
    stmp = AtomFactorS(i)
    fss[idx], fso[idx], foo[idx] = stmp**2, stmp*otmp, otmp**2

# 散乱ベクトルQ導出。Q:{111}, P{022} それぞれのピーク
x = open("xrd.dat", "w")
x.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\n")

# COORDINATION を使用して、従来のDISTANCE, CUSTOM, COMBINE, CUSTOM を1行で実現
# 1. groupの定義
x.write("SI: GROUP ATOMS=1-{:d}".format(int(args.atoms/3)))

# 3. COORDINATION 本体の定義
# ここでは、sin(x)/x*Lorch を定義する
coo = "coo{}: COORDINATION GROUPA=SI SWITCH={}"
if args.lorch:
    x.write("#MUST USE lmp_shik_nolorch\n")
    FUNC = "sin({:.4f}*x)*sin({:.4f}*x)/x/x*{:.4f} R_0=1"
    Rmax, PI = 25, np.pi
    for idx, i in enumerate(Q):
        # switch = "{CUSTOM FUNC=" + FUNC.format(
        #     i, PI/Rmax, 2*fss[idx]*Rmax/PI/args.atoms/i) + "}\n" # old ver.
        switch = "{CUSTOM FUNC=" + FUNC.format(
            i, PI/Rmax, Rmax/PI/i) + "}\n"  # Q, PI/Rmax, Rmax/PI/Q
        x.write(coo.format(idx, switch))
else:
    exit("under construction...")
    # FUNC = "sin({:.3f}*x)/x*{:.4f} R_0=1"  # sin(Qx)/Qx * f^2/N
    # hoge = [fss0/Q0*2, fso0/Q0*2, foo0/Q0*2, fss1/Q1*2, fso1/Q1*2, foo1/Q1*2]
    # # lstAtom = [atoms/3, atoms, atoms/3*2, atoms/3, atoms, atoms/3*2]
    # # SWITCH = "{CUSTOM FUNC=" + FUNC.format(Q0, hoge[0]) + "}\n"
    # for idx, i in enumerate(hoge):
    #     # tmp = FUNC.format(Q[idx//2], hoge[idx])
    #     switch = "{CUSTOM FUNC=" + FUNC.format(
    #         Q[idx//3], hoge[idx]/atoms) + "}\n"
    #     w.write(coord[idx % 3].format(idx//3, switch))
    #     x.write(coord[idx % 3].format(idx//3, switch))

# COORD = fsi(Q)**2 + coord * 2*fsi(Q)**2/N の定義
COO = "COO{a}: CUSTOM ARG=coo{a} FUNC=x*{b}+{c} PERIODIC=NO\n"
for idx, i in enumerate(Q):
    coeff = fss[idx]
    x.write(COO.format(a=idx, b=2*coeff/args.atoms*3, c=coeff))

# output colvar
PRINT = "PRINT ARG=COO0"
for idx, i in enumerate(Q[1:]):
    PRINT += ",COO{}".format(idx+1)
PRINT += " FILE=COLVAR STRIDE=500"
x.write(PRINT)
x.close()
print("init.dat created.")
