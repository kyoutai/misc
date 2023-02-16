#!/usr/bin/env python
import argparse
import numpy as np
import re
import matplotlib.pyplot as plt

description = """describe"""
par = argparse.ArgumentParser(description=description)
par.add_argument('files', nargs='+')
par.add_argument('-o', "--out", help="outfile name (hoge.png)")
args = par.parse_args()
print(args.files)


def Displacement():
    SiO2 = 60.08                # g/mol
    mol = 6.022e23
    gram = SiO2 * 64 / mol
    cm3 = 1e-24                 # Å^-3 → cm^-3
    vol_ex = gram / cm3
    print(vol_ex)
    for f in args.files:
        body = open(f).read()
        # 0:Step 1:Dt 2:Time 3:Temp 4:Press 5:PotEng 6:KinEng 7:TotEng 8:Enthalpy 9:Density
        # 10:c_2(PKA以外の変位) 11:c_4(100ステップ毎のPKAの変位)
        data = re.findall("Press          Volume    \n(.+?)\nLoop",
                          body, re.DOTALL)
        data = np.array(" ".join(data).split()).reshape((-1, 7)).astype(float)
        # ax.plot(data[:, 1]/1e9, vol_ex / data[:, -1], "-")
        ax.plot(data[:, 1]/1e9, data[:, -1], "-")  # Volume


fig, ax = plt.subplots()
Displacement()
ax.set_xlabel('time microsec')
ax.set_ylabel('Density')
ax.set_ylabel('Volume')
plt.tight_layout()
if args.out:
    plt.savefig(args.out)
    print(args.out, "was created.")
else:
    plt.show()
