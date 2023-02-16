#!/usr/bin/env python
import numpy as np
import argparse
import os
from scipy.optimize import leastsq
from scipy.optimize import differential_evolution
from scipy.optimize import basinhopping

par = argparse.ArgumentParser(description="test")
par.add_argument('trjfile',
                 help="lammpstrj file")
par.add_argument('-m', '--method', default="leastsq",
                 choices=['leastsq', 'basinhopping', 'differential_evolution'],
                 help="leastsq")
par.add_argument('-s', '--step', default=-1, type=int,
                 help="Target time step in lammpstrj [-1]")
par.add_argument('-o', '--outlog', default=None,
                 help="output logfile name [basename]")
par.add_argument('--p0', '--p0', nargs=3, default=[0, 0, 0], type=float,
                 help="initial guess of origin")
par.add_argument('-v', '--verbose', default=False, action='store_true',
                 help="verbosity")
args = par.parse_args()

SUPER_CELL = [5, 4, 5]
UNIT = """
Si  0.472000     0.000000     0.666667
Si  0.972000     0.500001     0.666667
Si  0.764000     0.236001     0.333333
Si  0.264000     0.736002     0.333333
Si  0.264000     0.264001     0.000000
Si  0.764000     0.764002     0.000000
O   0.275500     0.132500     0.786000
O   0.775500     0.632502     0.786000
O   0.663500     0.071500     0.452667
O   0.163500     0.571501     0.452667
O   0.561000     0.296001     0.119333
O   0.061000     0.796002     0.119333
O   0.061000     0.204001     0.214000
O   0.561000     0.704002     0.214000
O   0.775500     0.367501     0.547333
O   0.275500     0.867502     0.547333
O   0.163500     0.428501     0.880667
O   0.663500     0.928502     0.880667
"""
# alpha quarts(trigonal)
# Si  0.47200    0.00000    0.66667
# Si  0.00000    0.47200    0.33333
# Si  0.52800    0.52800    0.00000
# O   0.40800    0.26500    0.78600
# O   0.73500    0.14300    0.45267
# O   0.85700    0.59200    0.11933
# O   0.26500    0.40800    0.21400
# O   0.14300    0.73500    0.54733
# O   0.59200    0.85700    0.88067

CENTER_ATOM = "Si"
MAX_DISPLACEMNET = {}
MAX_DISPLACEMNET['Si'] = 0.6
MAX_DISPLACEMNET['O'] = 0.9
BONDS = {}
BONDS['O-Si'] = 1.75
BONDS['Si-O'] = 1.75


class XtalData():
    def __init__(self):
        self.udata = {}
        self.MakeUnit()

    def MakeUnit(self):
        # make SUPER_CELL
        unit = np.array(UNIT.split(), dtype=object).reshape((-1, 4))
        unit[:, 1:4] = unit[:, 1:4].astype(float)
        data = np.empty((0, 4))
        for i in range(SUPER_CELL[0]):
            for j in range(SUPER_CELL[1]):
                for k in range(SUPER_CELL[2]):
                    u = np.copy(unit)
                    u[:, 1:4] = u[:, 1:4] + [i, j, k]
                    data = np.vstack((data, u))
        data[:, 1:4] = data[:, 1:4]/[10, 10, 9]
        self.data = data
        # return data as fractional coord.


class LammpsData():
    def __init__(self, trjfile):
        self.trjfile = trjfile
        self.GetData()

    def GetData(self):
        """
        self.data:  fractional coord.
        axis=0: step
        axis=1: atoms
        axis=2: id type elem x y z (6で固定)
        self.M  座標変換の行列
        axis=0: step
        axis=1: 3
        axis=2: 3
        """
        lines = open(self.trjfile).readlines()
        atoms = int(lines[3])
        steps = int(len(lines)/(9 + atoms))
        rows = len(lines[10].split())
        body = np.array(lines).reshape((steps, 9 + atoms))
        L = body[:, 5:8]
        data = body[:, 9:9+atoms]
        L = np.array(" ".join(L.reshape(-1)).split())
        L = L.reshape((steps, 3, 2)).astype(float)
        # xlo_bound, xhi_bound = L[:, 0, 0], L[:, 0, 1]
        # ylo_bound, yhi_bound = L[:, 1, 0], L[:, 1, 1]
        # zlo_bound, zhi_bound = L[:, 2, 0], L[:, 2, 1]
        data = np.array(" ".join(data.reshape(-1)).split(), dtype=object)
        data = data.reshape((steps, atoms, rows))
        data = data[:, :, 0:6]  # remove vx, vy, vz
        data[:, :, 0:2] = data[:, :, 0:2].astype(int)    # id, type
        data[:, :, 3:6] = data[:, :, 3:6].astype(float)  # xyz
        a, b, c = np.empty(steps), np.empty(steps), np.empty(steps)
        # alpha, beta, gamma = np.empty(steps), np.empty(steps), np.empty(steps)
        # for i in range(steps):
        #     xlo = xlo_bound[i] - np.min([0.0, xy[i], xz[i], xy[i]+xz[i]])
        #     xhi = xhi_bound[i] - np.max([0.0, xy[i], xz[i], xy[i]+xz[i]])
        #     ylo = ylo_bound[i] - np.min([0.0, yz[i]])
        #     yhi = yhi_bound[i] - np.max([0.0, yz[i]])
        #     zlo = zlo_bound[i]
        #     zhi = zhi_bound[i]
        #     lx, ly, lz = xhi-xlo, yhi-ylo, zhi-zlo
        #     a[i] = lx
        #     b[i] = np.sqrt(ly**2 + xy[i]**2)
        #     c[i] = np.sqrt(lz**2 + xz[i]**2 + yz[i]**2)
        #     alpha[i] = np.arccos((xy[i] * xz[i] + ly * yz[i])/(b[i]*c[i]))
        #     beta[i] = np.arccos(xz[i]/c[i])
        #     gamma[i] = np.arccos(xy[i]/b[i])
        # M = np.zeros((steps, 3, 3))
        # M[:, 0, 0] = a
        # M[:, 1, 0] = b * np.cos(gamma)
        # M[:, 1, 1] = b * np.sin(gamma)
        # M[:, 2, 0] = c * np.cos(beta)
        # M[:, 2, 1] = c / np.sin(gamma) * (np.cos(alpha) -
        #                                   np.cos(gamma)*np.cos(beta))
        # M[:, 2, 2] = c * np.sqrt(1 - np.cos(beta)**2 -
        #                          1/np.sin(alpha)**2 *
        #                          (np.cos(alpha)-np.cos(gamma)*np.cos(beta))**2)
        # M_ = np.zeros(M.shape)
        for i in range(steps):
            # M_ = np.linalg.inv(M[i])
            # data[i, :, 3:6] = np.dot(data[i, :, 3:6], M_)
            data[i, :, 3:6] = data[i, :, 3:6] - np.floor(data[i, :, 3:6])
        self.data = data
        # self.M = M
        self.a, self.b, self.c = a, b, c
        # self.alpha = alpha * 180/np.pi
        # self.beta = beta * 180/np.pi
        # self.gamma = gamma * 180/np.pi


class GetCrystallinity():
    def __init__(self, xtal, lmp):
        self.xtal = xtal
        self.lmp = lmp
        if self.lmp.data.shape[1] != self.xtal.data.shape[0]:
            print("atom number in lammps and cif is inconsistent!!!")
            exit()

    def NearestAtoms(self, p):
        """
        Parameters:
        p: origin (0, 0)
        Return:
        rdix:  index list of nearest crystal site
        dr: distance list of nearest crystal site
        """
        dr = self.Si_lmp2D - self.Si_xtl2D + p
        dr = dr - np.round(dr)
        # delta = np.dot(dr, self.M)
        # dr = np.sum(delta**2, axis=2)
        dr = np.sum(dr**2, axis=2)
        ridx = np.argmin(dr, axis=0)
        delta = self.Si_lmp[:, :] - self.Si_xtl[ridx, :] + p
        delta = delta - np.round(delta)
        # dr = np.dot(delta, self.M)
        if args.verbose is True:
            print("sum(dr), origin", np.sum(dr**2), p)
        return ridx, dr

    def Error(self, p):
        ridx, dr = self.NearestAtoms(p)
        return(np.sum(dr**2))

    def Error_vec(self, p):
        ridx, dr = self.NearestAtoms(p)
        # (y-ydata)のvectorを返す
        return dr.reshape(-1).astype(float)

    def CallBack(self, x, f, accepted):
        print(x, f, int(accepted))

    def GetOrigin(self, step):
        # 2D data
        # 巨大メモリになるが2次元にして計算をはやくする
        atoms = self.Si_lmp.shape[0]
        self.Si_lmp2D = np.tile(self.Si_lmp, (atoms, 1))
        self.Si_lmp2D = self.Si_lmp2D.reshape((atoms, atoms, 3))
        self.Si_xtl2D = np.repeat(self.Si_xtl, atoms, axis=0)
        self.Si_xtl2D = self.Si_xtl2D.reshape((atoms, atoms, 3))
        # leastsqの時はvectorを返すfuncを指定
        p0 = np.array(args.p0)
        bounds = np.array([[0, 1], [0, 1], [0, 1]])
        if args.method == "leastsq":
            ret = p0
        if args.method == "differential_evolution":
            ret = differential_evolution(self.Error, bounds=bounds,
                                         maxiter=3).x
        # 局所minimumが不安ならbasinhoppingを使う。
        if args.method == "basinhopping":
            minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds}
            ret = basinhopping(self.Error, p0, niter=100,
                               niter_success=None,
                               callback=self.CallBack, disp=True,
                               minimizer_kwargs=minimizer_kwargs)
            ret = ret.x
        if args.verbose is True and args.method != "leastsq":
            print("stochastic:", ret, self.Error(ret))
        ret = leastsq(self.Error_vec, ret, xtol=1e-6)[0]
        print("Origin setting: {} ({})".format(ret, self.Error(ret)))
        self.origin = ret

    def CountCrystal(self, elem, cutoff):
        lmpdata = self.lmpdata[self.lmpdata[:, 2] == elem, :]
        lmpxyz = lmpdata[:, 3:6].astype(float)
        xtlxyz = self.xtal.data[self.xtal.data[:, 0] == elem, 1:4]
        xtlxyz = (xtlxyz - self.origin).astype(float)
        atoms = xtlxyz.shape[0]
        lmpdata2 = np.tile(lmpxyz, (atoms, 1)).reshape((atoms, atoms, 3))
        xtldata2 = np.repeat(xtlxyz, atoms, axis=0).reshape((atoms, atoms, 3))
        dr = lmpdata2 - xtldata2
        dr = dr - np.round(dr)
        print(dr.shape)
        delta = np.linalg.norm(dr, axis=2)
        rlist = np.min(delta, axis=0)
        flag = rlist < cutoff
        # crystal atoms
        # print("crystal:", elem, lmpdata[flag, :])
        # print("non-cyrstal:", elem, lmpdata[np.logical_not(flag), :])
        return lmpdata[flag, :]

    def CalcCrystal(self, step):
        self.lmpdata = self.lmp.data[step, :, :]
        # self.M = self.lmp.M[step]
        # only Si
        self.Si_lmp = self.lmpdata[self.lmpdata[:, 2] == CENTER_ATOM, 3:6]
        self.Si_xtl = self.xtal.data[self.xtal.data[:, 0] == CENTER_ATOM, 1:4]
        self.Si_lmp = self.Si_lmp.astype(float)
        self.Si_xtl = self.Si_xtl.astype(float)
        # debug
        # self.origin = np.array([-5.381e-08, -4.6977e-08, -3.70370e-02])
        print("Setting origin to be the highest crystalline atoms...")
        self.GetOrigin(step)
        crystal = np.empty((0, 6))  # id, type, elem, x, y, z
        elems = np.unique((self.lmpdata[:, 2]))
        # define crystal atoms from 結晶サイトの許容距離
        for e in elems:
            data = self.CountCrystal(e, MAX_DISPLACEMNET[e])
            crystal = np.vstack((crystal, data))

        print("Detecting isolated atoms...")
        cn = np.zeros(crystal.shape[0], dtype=int)
        for bond, cutoff in BONDS.items():
            elem_a, elem_b = bond.strip().split("-")
            idx_a = np.where(crystal[:, 2] == elem_a)[0]
            idx_b = np.where(crystal[:, 2] == elem_b)[0]
            for a in idx_a:
                dr = (crystal[a, 3:6] - crystal[idx_b, 3:6]).astype(float)
                dr = dr - np.round(dr)
                # dr = np.dot(dr, self.M)
                dr = np.linalg.norm(dr, axis=1)
                cn[a] += np.sum(dr < cutoff)
        cryslta_flag = cn != 0
        isolate_flag = cn == 0
        self.crystal = crystal[cryslta_flag, :]
        self.isolate = crystal[isolate_flag, :]
        self.crystallinity = self.crystal.shape[0]/self.lmp.data.shape[1]*100

    def OutputLog(self, step):
        if args.outlog is None:
            base, ext = os.path.splitext(args.trjfile)
            outfile = base + ".xtalog"
        else:
            outfile = args.outlog
        print("crystal atoms: {}".format(self.crystal.shape[0]))
        print("isolate atoms: {}".format(self.isolate.shape[0]))
        o = open(outfile, "w")
        o.write("# Origin setting\n")
        o.write("# Initial guess {} {} {}\n".format(*args.p0))
        o.write("# Optimized set {} {} {}\n".format(*self.origin))
        o.write("\n")
        o.write("# Atomic information\n")
        o.write("Total atoms {}\n".format(self.lmp.data.shape[1]))
        o.write("Crystalline atoms {}\n".format(self.crystal.shape[0]))
        o.write("Isolated atoms {}\n".format(self.isolate.shape[0]))
        o.write("\n")
        o.write("# Crystallinity\n")
        o.write("{:.4f} %\n\n".format(self.crystallinity))
        o.write("# Crystal system\n")
        o.write("a: {}\n".format(self.lmp.a[step]))
        o.write("b: {}\n".format(self.lmp.b[step]))
        o.write("c: {}\n".format(self.lmp.c[step]))
        # o.write("alpha: {}\n".format(self.lmp.alpha[step]))
        # o.write("beta:  {}\n".format(self.lmp.beta[step]))
        # o.write("gamma: {}\n".format(self.lmp.gamma[step]))
        fmt = "{:d} {:d} {:s} {:g} {:g} {:g}\n"
        s = "\n# Crystal without isolated atom index Total {} atoms\n"
        o.write(s.format(self.crystal.shape[0]))
        for a in self.crystal:
            o.write(fmt.format(*a))
        o.write("\n# Crystal but Isolated atom index Total {} atoms\n".format
                (self.isolate.shape[0]))
        for a in self.isolate:
            o.write(fmt.format(*a))
        print(outfile, "was created.")


xtal = XtalData()
lmp = LammpsData(args.trjfile)
crt = GetCrystallinity(xtal, lmp)
crt.CalcCrystal(args.step)
crt.OutputLog(args.step)
