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


# 23/01/16 okagiriで計算していたときの、1CVパラメータ導出方法
# print(AtomFactorS(22))
# fss = AtomFactorS(22)
# Q = 4*np.pi/1.5406 * np.sin((22/2)/180*np.pi)
# print(Q)
# print(fss**2/Q/192*3*2)
# exit()

Q, fss = np.empty_like(theta), np.empty_like(theta)
fso, foo = np.empty_like(theta), np.empty_like(theta)
for idx, i in enumerate(theta):
    Q[idx] = 4*np.pi/1.5406 * np.sin((i/2)/180*np.pi)
    otmp = AtomFactorO(i)
    stmp = AtomFactorS(i)
    fss[idx], fso[idx], foo[idx] = stmp**2, stmp*otmp, otmp**2

# 散乱ベクトルQ導出。Q:{111}, P{022} それぞれのピーク
x = open("xrd.dat", "w")
# x.write("#Natoms={}, alpha_quartz={}, include_Lorch={}, theta={}, {}\n".format(
#     args.atoms, str(args.alpha), str(args.lorch), args.theta[0], args.theta[1]))
x.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\n")

# COORDINATION を使用して、従来のDISTANCE, CUSTOM, COMBINE, CUSTOM を1行で実現
# 1. groupの定義
if args.alpha:
    sihead = "SI: GROUP ATOMS=1-8"
    ohead = "O: GROUP ATOMS=9-24"
    loops = int(args.atoms/18)
    for i in range(loops):
        if i:
            sihead += ",{}-{}".format(i*18+1, i*18+6)
            ohead += ",{}-{}".format(i*18+7, i*18+18)
    x.write(sihead + "\n" + ohead + "\n")
else:
    sihead = "SI: GROUP ATOMS=1-8"
    ohead = "O: GROUP ATOMS=9-24"
    loops = args.atoms/24
    for i in range(int(loops)):
        if i:
            sihead += ",{}-{}".format(i*24+1, i*24+8)
            ohead += ",{}-{}".format(i*24+9, i*24+24)
    x.write(sihead + "\n" + ohead + "\n")

# 3. COORDINATION 本体の定義
cooSS = "ss{}: COORDINATION GROUPA=SI SWITCH={}"
cooSO = "so{}: COORDINATION GROUPA=SI GROUPB=O SWITCH={}"
cooOO = "oo{}: COORDINATION GROUPA=O SWITCH={}"
# Q, coord = [Q0, Q1], [cooSS, cooSO, cooOO]
if args.lorch:
    x.write("#MUST USE lmp_shik_nolorch\n")
    FUNC = "sin({:.4f}*x)*sin({:.4f}*x)/x/x*{:.4f} R_0=1"
    Rmax, PI = np.cbrt(args.atoms/24)*7.5*np.sqrt(3), np.pi
    for idx, i in enumerate(Q):
        switch = "{CUSTOM FUNC=" + FUNC.format(
            i, PI/Rmax, 2*fss[idx]*Rmax/PI/args.atoms/i) + "}\n"
        x.write(cooSS.format(idx, switch))
        switch = "{CUSTOM FUNC=" + FUNC.format(
            i, PI/Rmax, 2*fso[idx]*Rmax/PI/args.atoms/i) + "}\n"
        x.write(cooSO.format(idx, switch))
        switch = "{CUSTOM FUNC=" + FUNC.format(
            i, PI/Rmax, 2*foo[idx]*Rmax/PI/args.atoms/i) + "}\n"
        x.write(cooOO.format(idx, switch))
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

cd = "cd{a}: CUSTOM ARG=ss{a},so{a},oo{a} FUNC=x+y+z+{b} PERIODIC=NO\n"
for idx, i in enumerate(Q):
    x.write(cd.format(a=idx, b=fss[idx]/3+foo[idx]*2/3))

# output colvar
PRINT = "PRINT ARG=cd0"
for idx, i in enumerate(Q[1:]):
    PRINT += ",cd{}".format(idx+1)
PRINT += " FILE=COLVAR STRIDE=500"
x.write(PRINT)
x.close()
print("init.dat created.")
