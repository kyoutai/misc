#!/usr/bin/env python
import argparse
import numpy as np
import re
import matplotlib.pyplot as plt

description = """describe"""
par = argparse.ArgumentParser(description=description)
par.add_argument('files', nargs='+')
par.add_argument('-n', "--natoms", help="number of atoms", type=int)
par.add_argument('-o', "--out", help="outfile name (hoge.png)")
par.add_argument(
    '-m', '--mode', choices=['density', 'displacement'], default='density')
args = par.parse_args()
print(args.files)


def Displacement(files):
    disp_max = []
    for f in args.files:
        body = open(f).read()
        # 0:Step 1:Dt 2:Time 3:Temp 4:Press 5:PotEng 6:KinEng 7:TotEng 8:Enthalpy 9:Density
        # 10:c_2(PKA以外の変位) 11:c_4(100ステップ毎のPKAの変位)
        data = re.findall(
            "Step Dt.+?c_2 c_4 \n(.+?)\nLoop time", body, re.DOTALL)
        n = len(data)  # 1000ステップの何回繰り返したか
        data = np.array(" ".join(data).split()).reshape((-1, 12)).astype(float)
        dmax = np.max(data[:, 10])
        disp_max.append(dmax)
        print("{}: {:5.2f} Ang.".format(f, dmax))
    bins = np.arange(5, 80, 1)
    height, bins = np.histogram(disp_max, bins=bins)
    x = bins[0:-1] + 0.5*(bins[1:] - bins[:-1])
    ax.bar(x, height, align="center")
    ax.plot(x, height, "ro-", ms=5)
    ax.set_xlabel("${\\rm Displacement/\AA}$")
    ax.set_ylabel("${\\rm Frequency}$")
    plt.show()


def Density(f):
    if len(f) == 1:
        body = open(f[0]).read()
        n = body.count("irradiation energy: 30 ev")
        print(n)
        data = re.findall(
            "Step           Dt.+?c_2      \n(.+?)\nLoop time", body, re.DOTALL)
        n = len(data)
        print("{} times irradiation.".format(n/2))
        data = np.array(" ".join(data).split()).reshape((n, -1, 11)
                                                        ).astype(float)
        enthalpy = np.mean(data[:, :, 8], axis=1)
        density = np.mean(data[:, :, 9], axis=1)
        density_ = np.std(data[:, :, 9], axis=1)
        x = np.linspace(1, n, n)/2
        bx = ax.twinx()
        ax.errorbar(x, density, yerr=density_, fmt="bo-")
        bx.plot(x, enthalpy/conv[0], "r-")
        ax.set_xlabel('${\\rm Number\ of\ irradiation}$')
        ax.set_ylabel('${\\rm Density\ g/cm^3}$')
        bx.set_ylabel('${\\rm Enthalpy\ eV/atom}$')
        ax.set_ylim(1.8, 2.8)
        ax.set_title(args.files[0].split()[0])
    else:
        for idx, i in enumerate(f):
            body = open(i).read()
            # data = re.findall(
            #         "Step Dt.+?c_2 c_4 \n(.+?)\nLoop time", body, re.DOTALL)
            data = re.findall(
                "Step           Dt.+?c_2      \n(.+?)\nLoop time", body, re.DOTALL)
            n = len(data)
            print("{} times irradiation.".format(n/2))
            data = np.array(" ".join(data).split()).reshape((n, -1, 11)
                                                            ).astype(float)
            enthalpy = np.mean(data[:, :, 8], axis=1)
            # enthalpy_ = np.std(data[:, :, 8], axis=1)
            density = np.mean(data[:, :, 9], axis=1)
            density_ = np.std(data[:, :, 9], axis=1)
            x = np.linspace(1, n, n)/2
            # fig, ax = plt.subplots()
            bx = ax[idx].twinx()
            ax[idx].errorbar(x, density, yerr=density_, fmt="bo-")
            print(conv[idx])
            # bx.errorbar(x, enthalpy, yerr=enthalpy_, fmt="rs-")
            bx.plot(x, enthalpy/conv[idx], "r-")
            # bx.grid()
            ax[idx].set_xlabel('${\\rm Number\ of\ irradiation}$')
            ax[idx].set_ylabel('${\\rm Density\ g/cm^3}$')
            bx.set_ylabel('${\\rm Enthalpy\ eV/atom}$')
            ax[idx].set_ylim(1.8, 2.8)
            ax[idx].set_title(args.files[idx].split()[0])


conv = 1.60217e-19 * 6.023e23 / 4186.05
# conv = [conv*1080, conv*192, conv*648]
conv = [conv * args.natoms]
fig, ax = plt.subplots(nrows=len(args.files), figsize=(5, 3*len(args.files)))
Density(args.files)
plt.tight_layout()
if args.out:
    plt.savefig(args.out)
    print(args.out, "was created.")
else:
    plt.show()
