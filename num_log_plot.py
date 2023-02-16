#!/usr/bin/env python
import argparse
import numpy as np
import re
import matplotlib.pyplot as plt

description = """何もプロットされないならば、正規表現に間違いあるのでプログラムを書き換えましょう。"""
par = argparse.ArgumentParser(description=description)
par.add_argument("ev", help="irradiation energy", type=int)
par.add_argument('files', nargs='+')
par.add_argument('-n', "--natoms", help="number of atoms", type=int)
par.add_argument('-o', "--out", help="outfile name (hoge.png)")
args = par.parse_args()


def Density(f, ev):
    for idx, i in enumerate(f):
        body = open(i).read()
        data = re.findall(
            "Step           Dt.+?c_2      \n(.+?)\nLoop time", body, re.DOTALL)
        # "Step Dt.+?c_2 \n(.+?)\nLoop time", body, re.DOTALL)
        data = np.array(" ".join(data).split()).reshape((-1, 11)
                                                        ).astype(float)
        enthalpy = np.mean(data[-50:, 8])
        density = np.mean(data[-50:, 9])
        x = idx * ev
        # print(density)
        bx.plot(x, enthalpy*conv/args.natoms, "r.")
        ax.plot(x, density, "b.")


conv = 1/23.060  # 1 eV = 23.060 kcal/mol
fig, ax = plt.subplots()
bx = ax.twinx()
onestep_energy = args.ev / args.natoms
Density(args.files, onestep_energy)
bx.set_ylabel('${\\rm Enthalpy\ eV/atom}$')
ax.set_xlabel('${\\rm Deposit\ energy\ eV/atom}$')
ax.set_ylabel('${\\rm Density\ g/cm^3}$')
ax.set_ylim(2.0, 2.7)
ax.set_title(args.files[0].split("_")[0]+"-"+str(args.ev)+" ev/time")
plt.tight_layout()
if args.out:
    plt.savefig(args.out)
    print(args.out, "was created.")
else:
    plt.show()
