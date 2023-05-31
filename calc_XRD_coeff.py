#!/usr/bin/env python
import numpy as np
import argparse
par = argparse.ArgumentParser(description="bar")
par.add_argument('deg', help='2theta, degree', type=int)
par.add_argument('atoms', help='number of atoms', type=int)
args = par.parse_args()

lam = 1.5406
Q = 4 * np.pi * np.sin(args.deg / 360 * np.pi) / lam
Ar = (Q / 4 / np.pi)**2  # Ar = (q / 4pi) ** 2
fsi = (5.275329*np.exp(-2.631338*Ar) + 3.191038*np.exp(-33.73073*Ar) +
       1.511514*np.exp(-0.081119*Ar) + 1.356849*np.exp(-86.28864*Ar) +
       2.519114*np.exp(-1.170087*Ar) + 0.145073)
fo = (2.960427*np.exp(-14.18226*Ar) + 2.508818*np.exp(-5.936858*Ar) +
      0.637853*np.exp(-0.112726*Ar) + 0.722838*np.exp(-34.95848*Ar) +
      1.142756*np.exp(-0.390240*Ar) + 0.027014)
h = "XRD: CUSTOM ARG=ss,so,oo \
FUNC=(x*{:.4f}+y*{:.4f}+z*{:.4f}+{:.4f}) PERIODIC=NO"
N = args.atoms/2
fss, fso, foo = fsi**2/N, fsi*fo/N, fo**2/N
base = fsi**2/3 + fo**2/3*2
print(args.atoms, "atoms")
print("Q = {:.4f}".format(Q))
print(h.format(fss, fso, foo, base))
