#!/usr/bin/env python
import numpy as np
import argparse
# https://aip.scitation.org/doi/10.1063/1.5081040

description = "Na金属のXRDピーク導出プログラム。plumedでpythonと同じ計算することができるか？デバッグ用"
description = "NaXRD.py 250"
par = argparse.ArgumentParser(description=description)
par.add_argument('atoms', help='number of atoms', default="250", type=int)
# par.add_argument("--alpha", action='store_true',
#                  help="alpha quartz 専用オプション。group を変えます")
# par.add_argument("--lorch", action='store_true', help="use lorch function")

args = par.parse_args()
theta = np.linspace(20, 80, 61)


def AtomFactorNa(ang):            # ang: 2theta
    lbd = 1.5406
    q1 = (np.sin(ang/360*np.pi)/lbd)**2
    fna = (4.910127*np.exp(-3.281434*q1)+3.081783*np.exp(-9.119178*q1) +
           1.262067*np.exp(-0.102763*q1)+1.098938*np.exp(-132.0139*q1) +
           0.560991*np.exp(-0.405878*q1)+0.0797120)
    return fna


Q, fna = np.empty_like(theta), np.empty_like(theta)
for idx, i in enumerate(theta):
    Q[idx] = 4*np.pi/1.5406 * np.sin((i/2)/180*np.pi)
    # Ntmp = AtomFactorNa(i)
    # fss[idx] = Ntmp**2
    fna[idx] = AtomFactorNa(i)**2

# 単位を units metal に揃えた。eV, A, ps
x = open("xrd.dat", "w")
x.write("UNITS ENERGY=eV LENGTH=A\n")

# COORDINATION を使用して、従来のDISTANCE, CUSTOM, COMBINE, CUSTOM を1行で実現
# 1. groupの定義
x.write("#MUST USE lmp_eam\n")
x.write("NA: GROUP ATOMS=1-{:d}\n".format(int(args.atoms)))

# 3. COORDINATION 本体の定義
# ここでは、sin(x)/x*Lorch を定義する
# cutoff = Rmax
coo = "coo{}: COORDINATION GROUPA=NA SWITCH={}"
FUNC = "sin({:.4f}*x)*sin({:.4f}*x)/x/x*{:.4f} R_0=1 D_MAX={}"
Rmax, PI = 11, np.pi
for idx, i in enumerate(Q):
    switch = "{CUSTOM FUNC=" + FUNC.format(
        i, PI/Rmax, Rmax/PI/i, Rmax) + "}\n"  # Q, PI/Rmax, Rmax/PI/Q
    x.write(coo.format(idx, switch))

# COORD = fna(Q)**2 + coord * 2*fsi(Q)**2/N の定義
COO = "COO{a}: CUSTOM ARG=coo{a} FUNC=x*{b}+{c} PERIODIC=NO\n"
for idx, i in enumerate(Q):
    coeff = fna[idx]            # coeff = fna ** 2
    x.write(COO.format(a=idx, b=2*coeff/args.atoms, c=coeff))

# output colvar
PRINT = "PRINT ARG=COO0"
for idx, i in enumerate(Q[1:]):
    PRINT += ",COO{}".format(idx+1)
PRINT += " FILE=COLVAR STRIDE=500"
x.write(PRINT)
x.close()
print("xrd.dat created.")
