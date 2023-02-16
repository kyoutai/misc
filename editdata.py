#!/usr/bin/env python
import re
import numpy as np
import argparse

# atom_style full で計算、出力したlammps data file を, atom_style charge に対応
# する書式に加工するプログラムです。

par = argparse.ArgumentParser(description="bar")
par.add_argument('infile', help='input lammpsdata file')
par.add_argument('outfile', help='output file')
args = par.parse_args()

with open(args.infile) as f:
    lines = f.readlines()
with open(args.infile) as f:
    body = f.read()

natoms = int(lines[2].split()[0])
ntypes = int(lines[3].split()[0])
head = lines[:12+ntypes] + ['Atoms # charge\n', '\n']
tail = lines[-1*natoms-3:]
# mass = re.findall("Masses\n\n(.+?)\n\nPair", body, re.DOTALL)[0].split()
atom = re.findall("Atoms # full\n\n(.+?)\n\nVeloc", body, re.DOTALL)[0].split()
atom = np.array(atom, dtype=object).reshape(-1, 10)[:, 1:]
atom[:, :2] = atom[:, :2].astype(int)
atom[:, 2:] = atom[:, 2:].astype(float)
atom[:, -3:] = atom[:, -3:].astype(int)


fmt = "{0[0]:d} {0[1]:d} {0[2]:5f} {0[3]:20.14e} {0[4]:20.14e} {0[5]:20.14e}\
 {0[6]:d} {0[7]:d} {0[8]:d}\n"

with open(args.outfile, "w") as w:
    for i in head:
        w.write(i)
    for i in atom:
        w.write(fmt.format(list(i)))
    for i in tail:
        w.write(i)
