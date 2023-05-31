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
# par.add_argument(
#     '-m', '--mode', choices=['density', 'displacement'], default='density')
args = par.parse_args()


def Density(f):
    if len(f) == 1:
        body = open(f[0]).read()
        data = re.findall(
            # "Step Dt.+?c_2 \n(.+?)\nLoop time", body, re.DOTALL)
            # "Step           Time.+?Volume    \n(.+?)\nLoop time", body, re.DOTALL)
            "Step           Time.+?c_2    \n(.+?)\nLoop time", body, re.DOTALL)
        n = len(data)
        print(data)
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
                # "Step Dt.+?c_2 \n(.+?)\nLoop time", body, re.DOTALL)
                "Step           Time.+?Volume    \n(.+?)\nLoop time", body, re.DOTALL)
            n = len(data)
            data = np.array(" ".join(data).split())
            print(np.array(" ".join(data).split()).shape)
            print(args.files)
            # for i in data:
            #     print(data)
            # exit()

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


conv = 1.60217e-19 * 6.023e23 / 4186.05  # J/ev * mol^-1 / kcal/J -> kcal/mol
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
