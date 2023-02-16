#!/usr/bin/env python
import numpy as np
import argparse

# xyz ファイルを cif に変換するプログラム。
# 自前の処理を追加しているので、使用前に改変してください。

test = "usage: xyz2cif.py foo.xyz -l 10 10 9 -a 0 0 90"
par = argparse.ArgumentParser(description=test)
par.add_argument('xyz', help='input xyz file')
par.add_argument('lattice', help='[a, b, c] 格子定数を指定、**引数3つ必要**',
                 nargs=3, type=float)
par.add_argument('--ab', help='ab平面で指定角度(deg)だけ回転', type=float)
par.add_argument('--bc', help='bc平面で指定角度(deg)だけ回転', type=float)
par.add_argument('--ca', help='ca平面で指定角度(deg)だけ回転', type=float)
par.add_argument('--alpha', default=90, type=float,
                 help='lattice parameter, alpha(deg)')
par.add_argument('--beta', default=90, type=float,
                 help='lattice parameter, beta(deg)')
par.add_argument('--gamma', default=90, type=float,
                 help='lattice parameter, gamma(deg)')
par.add_argument('-o', '--outfile', required=True, help='output xyz file')
args = par.parse_args()

CIF = """#======================================================================
# CRYSTAL DATA
#----------------------------------------------------------------------
data_VESTA_phase_1

_chemical_name_common                  'test'
_cell_length_a                         %f
_cell_length_b                         %f
_cell_length_c                         %f
_cell_angle_alpha                      %f
_cell_angle_beta                       %f
_cell_angle_gamma                      %f
_space_group_name_H-M_alt              'P 1'
_space_group_IT_number                 1

loop_
_space_group_symop_operation_xyz
   'x, y, z'

loop_
   _atom_site_label
   _atom_site_occupancy
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_adp_type
   _atom_site_B_iso_or_equiv
   _atom_site_type_symbol
"""

a, b, c = args.lattice[0], args.lattice[1], args.lattice[2]
alpha = args.alpha*np.pi/180
beta = args.beta*np.pi/180
gamma = args.gamma*np.pi/180
# 規格化行列
v1 = [a, 0, 0]
v2 = [b*np.cos(gamma),  b*np.sin(gamma), 0]
v3 = [c*np.cos(beta),
      c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma),
      c*np.sqrt(1+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)
                - np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2 /
                np.sin(gamma))]
M = np.array([v1, v2, v3])
M_ = np.linalg.inv(M)

data = np.loadtxt(args.xyz, skiprows=2, dtype=object)
data[:, 1:] = data[:, 1:].astype(float)
# data[:, 1:4] = data[:, 1:4] + 0.5
# 回転行列
if args.ab:
    theta = args.ab*np.pi/180
    turnab = np.array([[np.cos(theta), -np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
    data[:, 1:3] = np.dot(data[:, 1:3], turnab)
    pass
if args.bc:
    pass
if args.ca:
    pass
# 格子定数で規格化
data[:, 1:4] = np.dot(data[:, 1:4], M_)
data[:, 1:4] = data[:, 1:4] - np.min(data[:, 1:4], axis=0)

w = open(args.outfile, "w")
w.write(CIF % (a, b, c, args.alpha, args.beta, args.gamma))
for d in data:
    w.write("{:s} 1.0 {:f} {:f} {:f} Biso 1.0 {:s}\n".format(
        d[0], d[1], d[2], d[3], d[0]))
print("done")
