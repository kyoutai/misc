#!/usr/bin/env python
import re
import numpy as np
import argparse

# 230215
# this program can be used to SIO2 trajectory file
# 出力したlammps trajectory data file を atom_style charge lammps data に変換
# 最終ステップを自動で読み込み data file に書き換える。
# ついでに、原子データをid 順番に列べる。

par = argparse.ArgumentParser(description="bar")
par.add_argument('body', help='input file body. EX: hoge -> hoge.lammpstrj')
args = par.parse_args()

charge = {"O": "-1.200", "Si": "+2.400"}
# readhead
headitems, bodyitems = 0, 0
with open(args.body+".lammpstrj") as f:
    for i in range(4):
        head = f.readline()
        headitems += len(head.split())
    atoms = int(head)           # NUMBER OF ATOMS_VALUE
    step1 = atoms + 9
    for i in range(5):
        ITEM = f.readline()     # ITEM: ATOMS id type ...
        headitems += len(ITEM.split())
    ITEM = ITEM.split()[2:]
    bodyitems = len(ITEM)
    xidx, yidx, zidx = 0, 0, 0
    for idx, item in enumerate(ITEM):
        if item == "x":
            xidx = idx
        if item == "y":
            yidx = idx
        if item == "z":
            zidx = idx
    if not (xidx and yidx and zidx):
        exit("wrong ITEM: ATOMS!! check your lammps input file.")

with open(args.body+".lammpstrj") as f:
    lines = f.readlines()

lines = np.array(" ".join(lines).split(), dtype=object).reshape(
    -1, headitems + bodyitems*atoms)[-1, :]
L = list(lines[14:20].astype(float))  # xmin, xmax, ymin, ymax, zmin, zmax

fmt = "  {0:>4d} {1[1]:>4d}   {2:s}   "
fmt += "{{{0}}}   ".format("1[{0}]:10.6e".format(xidx))
fmt += "{{{0}}}   ".format("1[{0}]:10.6e".format(yidx))
fmt += "{{{0}}}   ".format("1[{0}]:10.6e".format(zidx))
fmt += "# {{{0}}}\n".format("1[2]:s")

body = lines[headitems:].reshape(-1, bodyitems)
body[:, :2] = body[:, :2].astype(int)
body[:, 3:] = body[:, 3:].astype(float)
sortindex = np.argsort(body[:, 2])[::-1]  # O->Si にソートされるので [::-1]逆転

with open(args.body+".data", "w") as w:
    w.write("from {}.lammpstrj data file\n\n{} atoms\n2 atom types\n\n"
            .format(args.body, atoms))
    w.write("{}  {}  xlo xhi\n{}  {}  ylo yhi\n{}  {}  zlo zhi\n\n".format(*L))
    w.write("Masses\n\n   1   28.085500 # Si  Si\n   2   15.999400 # O  O\n\n")
    w.write("Atoms\n\n")
    for idx, i in enumerate(body[sortindex]):
        w.write(fmt.format(idx+1, i, charge[i[2]]))
