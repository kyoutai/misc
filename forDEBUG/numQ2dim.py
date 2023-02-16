#!/usr/bin/env python
import numpy as np
import argparse

description = """www.pnas.org/cgi/doi/10.1073/pnas.1803919115
cv2 XRDピーク計算が正しいかどうか検討するための実行ファイルを出力
COLVARcd, COLVARss, COLVARso, COLVARoo の4種類を出力できる
"""
par = argparse.ArgumentParser(description=description)
par.add_argument('atoms', help='number of atoms', default="192")
# par.add_argument("-P", "--PACE", default=1000)
# par.add_argument("-H", "--HEIGHT", default=14.4)
# par.add_argument("-S", "--SIGMA", default=2)
# par.add_argument("-B", "--BIAS", default=50)
# par.add_argument("-T", "--TEMP", default=2700)
# par.add_argument("-t", "--theta",
#                  nargs=2, default=[22, 34], help="2thata", type=float)

args = par.parse_args()


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


# 散乱ベクトルQ導出。Q:{111}, P{022} それぞれのピーク
lam = 1.5406
Ttheta = np.linspace(15, 50, 36)
Q = 4 * np.pi * np.sin(Ttheta / 360 * np.pi) / lam
# Q0 = 4*np.pi/1.5406 * np.sin((args.theta[0]/2)/180*np.pi)
# Q1 = 4*np.pi/1.5406 * np.sin((args.theta[1]/2)/180*np.pi)
# w = open(args.dir+"/plumed.dat", "w")
# x = open(args.dir+"/init.dat", "w")
w = open("plumed2dim.dat", "w")
x = open("init2dim.dat", "w")
w.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\nRESTART\n")
x.write("UNITS ENERGY=kcal/mol LENGTH=A TIME=fs\n")

# COORDINATION を使用して、従来のDISTANCE, CUSTOM, COMBINE, CUSTOM を1行で実現
# 1. groupの定義
group = """SI: GROUP ATOMS=1-8,25-32,49-56,73-80,97-104,121-128,145-152,169-176
O: GROUP ATOMS=9-24,33-48,57-72,81-96,105-120,129-144,153-168,178-192
"""
# Q{111} = 1.556, Q{022} = 2.385
w.write(group)
x.write(group)

#### start loop progress ####
for idx, q in enumerate(Q):
    w.write("# Q = {}\n".format(q))
    # 2. atom scattering factor 導出
    fss = AtomFactorS(Ttheta[idx]) * AtomFactorS(Ttheta[idx])
    fso = AtomFactorS(Ttheta[idx]) * AtomFactorO(Ttheta[idx])
    foo = AtomFactorO(Ttheta[idx]) * AtomFactorO(Ttheta[idx])
    # print(fss0, fso0, foo0, fss1, fso1, foo1)
    # exit()

    # 3. COORDINATION 本体の定義
    cooSS = "ss{}: COORDINATION GROUPA=SI SWITCH={}"
    cooSO = "so{}: COORDINATION GROUPA=SI GROUPB=O SWITCH={}"
    cooOO = "oo{}: COORDINATION GROUPA=O SWITCH={}"
    coord = [cooSS, cooSO, cooOO]
    CUSTOM0, CUSTOM1 = "{CUSTOM FUNC=", "}\n"
    FUNC = "sin({:.3f}*x)/x*{:.3f} R_0=1"
    # hoge = [fss0/Q0*2, fso0/Q0*2, foo0/Q0*2, fss1/Q1*2, fso1/Q1*2, foo1/Q1*2]
    hoge = [fss/q*2, fso/q*2, foo/q*2]
    SWITCH = CUSTOM0 + FUNC.format(q, hoge[0]) + CUSTOM1
    for IDX, i in enumerate(hoge):
        # tmp = FUNC.format(Q[idx//2], hoge[idx])
        switch = CUSTOM0 + FUNC.format(q, hoge[IDX]) + CUSTOM1
        w.write(coord[IDX % 3].format(idx, switch))
        x.write(coord[IDX % 3].format(idx, switch))

    cd = "cd{a}: CUSTOM ARG=ss{a},so{a},oo{a} FUNC=(x+y+z+{b})/192 PERIODIC=NO\n"
    w.write(cd.format(a=idx, b=int(fss*64+foo*128)))
    x.write(cd.format(a=idx, b=int(fss*64+foo*128)))

# metadynamics
# metad = """# 2cv metadynamics
# metad: METAD ARG=cd0,cd1 ...
#     PACE=500 HEIGHT=9.56 BIASFACTOR=100 SIGMA=5,5 FILE=HILLS TEMP=2800
#     GRID_MIN=-150,-40 GRID_MAX=+350,+160 GRID_BIN=1000,400
#     # GRID_MIN=-100,-40 GRID_MAX=+300,+120 GRID_BIN=800
# ...
# """
# w.write(metad)
# x.write(metad)
# output colvar
txtcd, txtss, txtso, txtoo = "", "", "", ""
for idx, i in enumerate(Q):
    txtcd += ",cd{}".format(idx)
    txtss += ",ss{}".format(idx)
    txtso += ",so{}".format(idx)
    txtoo += ",oo{}".format(idx)
PRINT = "PRINT ARG={} FILE=COLVAR{} STRIDE=50\n"
w.write(PRINT.format(txtcd[1:], "cd"))
x.write(PRINT.format(txtcd[1:], "cd"))
w.write(PRINT.format(txtss[1:], "ss"))
x.write(PRINT.format(txtss[1:], "ss"))
w.write(PRINT.format(txtso[1:], "so"))
x.write(PRINT.format(txtso[1:], "so"))
w.write(PRINT.format(txtoo[1:], "oo"))
x.write(PRINT.format(txtoo[1:], "oo"))
w.close()
x.close()
print("init.dat & plumed.dat created.")
