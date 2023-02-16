#!/usr/bin/env python
import numpy as np
import argparse
# www.pnas.org/cgi/doi/10.1073/pnas.1803919115
# cv2 metadynamics inputfile maker

# description = "cv値を2つのXRDピークとし、Si,O両方を考慮するPLUMED実行ファイルを生成するプログラム。本プログラムは、18原子のlammpsdataファイルを使用することを想定して作成した。"
description = "make2dimDat.py 324 -t 22 34"
par = argparse.ArgumentParser(description=description)
par.add_argument('atoms', help='number of atoms', default="192")
par.add_argument("-P", "--PACE", default=1000)
par.add_argument("-H", "--HEIGHT", default=14.4)
par.add_argument("-S", "--SIGMA", default=2)
par.add_argument("-B", "--BIAS", default=50)
par.add_argument("-T", "--TEMP", default=2700)
par.add_argument("-t", "--theta",
                 nargs=2, default=[22, 34], help="2thata", type=float)
par.add_argument("--alpha", action='store_true',
                 help="alpha quartz 専用オプション。group を変えます")
par.add_argument("--lorch", action='store_true', help="use lorch function")

args = par.parse_args()
atoms = int(args.atoms)


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

# 散乱ベクトルQ導出。Q:{111}, P{022} それぞれのピーク
Q0 = 4*np.pi/1.5406 * np.sin((args.theta[0]/2)/180*np.pi)
Q1 = 4*np.pi/1.5406 * np.sin((args.theta[1]/2)/180*np.pi)
# w = open(args.dir+"/plumed.dat", "w")
# x = open(args.dir+"/init.dat", "w")
w = open("plumed2dim.dat", "w")
x = open("init2dim.dat", "w")
w.write("#Natoms={}, alpha_quartz={}, include_Lorch={}, theta={}, {}\n".format(
    args.atoms, str(args.alpha), str(args.lorch), args.theta[0], args.theta[1]))
x.write("#Natoms={}, alpha_quartz={}, include_Lorch={}, theta={}, {}\n".format(
    args.atoms, str(args.alpha), str(args.lorch), args.theta[0], args.theta[1]))
w.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\nRESTART\n")
x.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\n")

# COORDINATION を使用して、従来のDISTANCE, CUSTOM, COMBINE, CUSTOM を1行で実現
# 1. groupの定義
if args.alpha:
    sihead = "SI: GROUP ATOMS=1-8"
    ohead = "O: GROUP ATOMS=9-24"
    loops = int(atoms/18)
    for i in range(loops):
        if i:
            sihead += ",{}-{}".format(i*18+1, i*18+6)
            ohead += ",{}-{}".format(i*18+7, i*18+18)
    w.write(sihead + "\n" + ohead + "\n# Q{111} = " + "{:5f}, Q".format(Q0)
            + "{022} = " + "{:5f}\n".format(Q1))
    x.write(sihead + "\n" + ohead + "\n# Q{111} = " + "{:5f}, Q".format(Q0)
            + "{022} = " + "{:5f}\n".format(Q1))
else:
    group = """SI: GROUP ATOMS=1-8,25-32,49-56,73-80,97-104,121-128,145-152,169-176
O: GROUP ATOMS=9-24,33-48,57-72,81-96,105-120,129-144,153-168,178-192
# Q{111} = """ + "{:5f}, Q".format(Q0) + "{022} = " + "{:5f}\n".format(Q1)
    w.write(group)
    x.write(group)

# 2. atom scattering factor 導出
fss0 = AtomFactorS(args.theta[0]) * AtomFactorS(args.theta[0])
fso0 = AtomFactorS(args.theta[0]) * AtomFactorO(args.theta[0])
foo0 = AtomFactorO(args.theta[0]) * AtomFactorO(args.theta[0])
fss1 = AtomFactorS(args.theta[1]) * AtomFactorS(args.theta[1])
fso1 = AtomFactorS(args.theta[1]) * AtomFactorO(args.theta[1])
foo1 = AtomFactorO(args.theta[1]) * AtomFactorO(args.theta[1])

x.write("# fsiQ{111} = " + "{:5f},".format(AtomFactorS(args.theta[0])))
x.write(" fsiQ{022} = " + "{:5f}\n".format(AtomFactorS(args.theta[1])))
x.write("# foQ{111} = " + "{:5f},".format(AtomFactorO(args.theta[0])))
x.write(" foQ{022} = " + "{:5f}\n".format(AtomFactorO(args.theta[1])))
# 3. COORDINATION 本体の定義
cooSS = "ss{}: COORDINATION GROUPA=SI SWITCH={}"
cooSO = "so{}: COORDINATION GROUPA=SI GROUPB=O SWITCH={}"
cooOO = "oo{}: COORDINATION GROUPA=O SWITCH={}"
Q, coord = [Q0, Q1], [cooSS, cooSO, cooOO]
if args.lorch:
    w.write("#MUST USE lmp_shik_nolorch\n")
    x.write("#MUST USE lmp_shik_nolorch\n")
    FUNC = "sin({:.4f}*x)*sin({:.4f}*x)/x/x*{:.4f} R_0=1"
    Rmax, PI = 16*3**0.5, np.pi
    hoge = [2*fss0*Rmax/PI/atoms/Q0, 2*fso0*Rmax/PI/atoms/Q0,
            2*foo0*Rmax/PI/atoms/Q0, 2*fss0*Rmax/PI/atoms/Q1,
            2*fso0*Rmax/PI/atoms/Q1, 2*foo0*Rmax/PI/atoms/Q1]
    for idx, i in enumerate(hoge):
        switch = "{CUSTOM FUNC=" + FUNC.format(
            Q[idx//3], PI/Rmax, hoge[idx]) + "}\n"
        w.write(coord[idx % 3].format(idx//3, switch))
        x.write(coord[idx % 3].format(idx//3, switch))
else:
    FUNC = "sin({:.3f}*x)/x*{:.4f} R_0=1"  # sin(Qx)/Qx * f^2/N
    hoge = [fss0/Q0*2, fso0/Q0*2, foo0/Q0*2, fss1/Q1*2, fso1/Q1*2, foo1/Q1*2]
    # lstAtom = [atoms/3, atoms, atoms/3*2, atoms/3, atoms, atoms/3*2]
    # SWITCH = "{CUSTOM FUNC=" + FUNC.format(Q0, hoge[0]) + "}\n"
    for idx, i in enumerate(hoge):
        # tmp = FUNC.format(Q[idx//2], hoge[idx])
        switch = "{CUSTOM FUNC=" + FUNC.format(
            Q[idx//3], hoge[idx]/atoms) + "}\n"
        w.write(coord[idx % 3].format(idx//3, switch))
        x.write(coord[idx % 3].format(idx//3, switch))

cd = "cd{a}: CUSTOM ARG=ss{a},so{a},oo{a} FUNC=x+y+z+{b} PERIODIC=NO\n"
w.write(cd.format(a=0, b=fss0/3+foo0*2/3))
x.write(cd.format(a=0, b=fss0/3+foo0*2/3))
w.write(cd.format(a=1, b=fss1/3+foo1*2/3))
x.write(cd.format(a=1, b=fss1/3+foo1*2/3))

# metadynamics
metad = """# 2cv metadynamics
metad: METAD ARG=cd0,cd1 ...
    PACE={} HEIGHT={} BIASFACTOR={} SIGMA={},{} FILE=HILLS TEMP={}
    GRID_MIN=-150,-40 GRID_MAX=+350,+160 GRID_BIN=1000,400
    # GRID_MIN=-100,-40 GRID_MAX=+300,+120 GRID_BIN=800
...
""".format(args.PACE, args.HEIGHT, args.BIAS, args.SIGMA, args.SIGMA,
           args.TEMP)
w.write(metad)
x.write(metad)
# output colvar
PRINT = "PRINT ARG=cd0,cd1 FILE=COLVAR STRIDE=500"
w.write(PRINT.format(args.PACE))
x.write(PRINT.format(args.PACE))
w.close()
x.close()
print("init.dat & plumed.dat created.")
