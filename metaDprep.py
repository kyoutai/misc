#!/usr/bin/env python
import numpy as np
import argparse

bar = """
メタダイナミクスに使用するインプットファイルを作成するプログラム
すぐに qsub run.sh で実行できる状態までもっていく。
作成するファイル：
  first.in : メタダイナミクス実行前の平衡構造作成
  init.in : メタダイナミクスの初期化
  init.dat : init.in のメタダイナミクス条件
  meta_morse.in : メタダイナミクス実行ファイル、ループ実行可能
  plumed.dat : meta_morse.in のメタダイナミクス条件
  run.sh : 計算を計算機に流すスクリプト
前もって準備するもの
  ???.mod 力場を記入したファイル
  ???.data MDの初期構造、照射構造が好ましい(溶融ガラス構造は別途スクリプト作成)
初期構造は、原子IDが 1-108:Si, 109-324:O で合計原子数324個とする。
"""
par = argparse.ArgumentParser(description=bar)
par.add_argument('name', help='directory name')
par.add_argument('temp', help='system temperature')
par.add_argument('potential', help='potential .mod file')
par.add_argument('data', help='lammps data file')
par.add_argument("-n", '--nodes', help='number of nodes', default=1, type=int)
par.add_argument(
    '--choice', choices=["silicon", "germanium", "tin"], default="germanium")
par.add_argument('-b', "--bias", help='bias factor', default=100, type=int)
par.add_argument('-l', "--loop", help='run loop (50 ns/loop)',
                 default=5, type=int)
args = par.parse_args()

IN = """log            {log}
boundary       p p p
units          real
atom_style     charge
{read}
replicate      1 1 1

include        ../{potential}
{reset_timestep}
timestep       1.0
thermo         10000
thermo_style   custom step time temp enthalpy etotal press vol
{plumed_file}
dump           1 all custom 100 {phase}.lammpstrj id type element x y z
dump_modify    1 element Si O sort id
"""

# write first.in
w = open("first.in", "w")
head = IN.format(log="first.log",
                 phase="first",
                 read="read_data      ../{}".format(args.data),
                 potential=args.potential,
                 reset_timestep="",
                 plumed_file="")
body = """velocity       all create 300 300
fix            2 all nvt temp 300  {T} 100

run            100000
fix            2 all nvt temp  {T}  {T} 100

run            30000
write_data     init.data
write_restart  init.rst
unfix          2
fix            2 all npt temp  {T}  {T} 100 iso 1 1 10000
run            30000
"""
w.write(head)
w.write(body.format(T=args.temp))
w.close()

# init.in
w = open("init.in", "w")
plumed = "fix            10 all plumed plumedfile init.dat outfile plumed.log"
head = IN.format(log="cv2.log",
                 phase="init",
                 read="read_restart   init.rst",
                 potential=args.potential,
                 reset_timestep="reset_timestep 0",
                 plumed_file=plumed)
body = """fix            2 all npt temp  {T}  {T} 100 iso 1 1 10000

run            10000
write_restart  cv2.rst
""".format(T=args.temp)
w.write(head)
w.write(body)
w.close()

# meta_morse.in
w = open("meta_morse.in", "w")
plumed = "fix            10 all plumed plumedfile plumed.dat outfile plumed.log"
head = IN.format(log="cv2.log append",
                 phase="metaD",
                 read="read_restart   cv2.rst",
                 potential=args.potential,
                 reset_timestep="",
                 plumed_file=plumed)
body = "fix            2 all npt temp  {T}  {T} 100 iso 1 1 10000".format(
    T=args.temp)
w.write(head)
w.write(body)
loop = "\nrun            50000000\nwrite_restart  cv2.rst"
w.write("# {} ns calculation\n".format(args.loop*50))
for i in range(args.loop):
    w.write(loop)
w.close()

DAT_0 = """#atoms=324, theta=27
UNITS ENERGY=kcal/mol LENGTH=A TIME=fs
{RESTART}
SI: GROUP ATOMS=1-108
O: GROUP ATOMS=109-324
"""
DAT_1 = """ss: COORDINATION GROUPA=SI SWITCH={CUSTOM FUNC=sin(1.9041678309861318*x)/1.9041678309861318/x*sin(0.1257*x)/0.1257/x R_0=1}
so: COORDINATION GROUPA=SI GROUPB=O SWITCH={CUSTOM FUNC=sin(1.9041678309861318*x)/1.9041678309861318/x*sin(0.1257*x)/0.1257/x R_0=1}
oo: COORDINATION GROUPA=O SWITCH={CUSTOM FUNC=sin(1.9041678309861318*x)/1.9041678309861318/x*sin(0.1257*x)/0.1257/x R_0=1}
RD: CUSTOM ARG=ss,so,oo FUNC=(x*0.7107+y*0.4270+z*0.2565+66.0829) PERIODIC=NO
"""
DAT_2 = """metad: METAD ARG=RD ...
    PACE=200 HEIGHT=9.56 BIASFACTOR={bias} SIGMA=5 FILE=HILLS TEMP={T}
    GRID_MIN=-100 GRID_MAX=+300
...
PRINT ARG=RD FILE=COLVAR STRIDE=200"""
# init.dat
w = open("init.dat", "w")
w.write(DAT_0.format(RESTART=""))
w.write(DAT_1)
w.write(DAT_2.format(T=args.temp, bias=args.bias))
w.close()
# plumed.dat
w = open("plumed.dat", "w")
w.write(DAT_0.format(RESTART="RESTART"))
w.write(DAT_1)
w.write(DAT_2.format(T=args.temp, bias=args.bias))
w.close()

SH = """#!/bin/sh
#PBS -V
#PBS -l nodes={node}:ppn={core}
#PBS -N {name}
#PBS -M gawagawa888@outlook.jp
#PBS -m abe
#PBS -j oe
#PBS -k oed
#PBS -q {typ}

cd $PBS_O_WORKDIR
OMP_NUM_THREADS=1
NPROCS=`wc -l < $PBS_NODEFILE`

lammps=/home/gawa/bin/lmp_{lmp}
mpirun -machinefile $PBS_NODEFILE -np {np} $lammps -in first.in
mpirun -machinefile $PBS_NODEFILE -np {np} $lammps -in init.in
mpirun -machinefile $PBS_NODEFILE -np {np} $lammps -in meta_morse.in
"""
w = open("run.sh", "w")
if args.choice == "silicon":
    w.write(SH.format(np=4, core=4, node=args.nodes, name=args.name,
                      lmp="plumed_silicon", typ="silicon"))
elif args.choice == "germanium":
    w.write(SH.format(np=10, core=10, node=args.nodes, name=args.name,
                      lmp="eam", typ="germanium"))
else:  # tin
    w.write(SH.format(np=10, core=24, node=args.nodes, name=args.name,
                      lmp="plumed_tin", typ="tin"))
w.close()
