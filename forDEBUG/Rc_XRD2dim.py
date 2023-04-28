#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

# http://www.pnas.org/cgi/doi/10.1073/pnas.1803919115
# SIO2 全原子でXRD 強度を求めるプログラム。
# 前作XRD2dim.pyと異なり、原子種類別のピークは出さない。
# 5種類のLorch関数使用したXRDピークをまとめて出すことができます。
par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-l', '--lattice', help='lattice parameter', default="14.3")
par.add_argument('-o', '--outfile', help='output name')
# par.add_argument('--onestep', action="store_true")
args = par.parse_args()

lattice = float(args.lattice)

# 図の体裁
plt_dic = {}
plt_dic['legend.fancybox'] = True
plt_dic['legend.labelspacing'] = 0.3
plt_dic['legend.numpoints'] = 3
plt_dic['figure.figsize'] = [8, 6]
plt_dic['axes.grid'] = True
plt_dic['font.size'] = 12
plt_dic['legend.fontsize'] = 12
plt_dic['axes.labelsize'] = 14
plt_dic['xtick.major.size'] = 5
plt_dic['xtick.minor.size'] = 3
plt_dic['xtick.direction'] = 'in'
plt_dic['savefig.bbox'] = 'tight'
plt_dic['savefig.dpi'] = 150
plt_dic['savefig.transparent'] = True
plt.rcParams.update(plt_dic)
fig, ax = plt.subplots()

# 原子散乱因子と、散乱定数Q を設定
lam = 1.5406
Ttheta = np.linspace(15, 80, 651)
Ttheta = np.linspace(15, 50, 351)
Q = 4 * np.pi * np.sin(Ttheta / 360 * np.pi) / lam

fsi = np.zeros_like(Ttheta)
fo = np.zeros_like(Ttheta)
for idx, q in enumerate(Q):  # atom scattering factors
    Ar = (q / 4 / np.pi)**2  # Ar = (q / 4pi) ** 2
    fsi[idx] = (5.275329*np.exp(-2.631338*Ar) + 3.191038*np.exp(-33.73073*Ar) +
                1.511514*np.exp(-0.081119*Ar) + 1.356849*np.exp(-86.28864*Ar) +
                2.519114*np.exp(-1.170087*Ar) + 0.145073)
    fo[idx] = (2.960427*np.exp(-14.18226*Ar) + 2.508818*np.exp(-5.936858*Ar) +
               0.637853*np.exp(-0.112726*Ar) + 0.722838*np.exp(-34.95848*Ar) +
               1.142756*np.exp(-0.390240*Ar) + 0.027014)


def Lorch(r, Rc):  # Window Func.
    return np.sin(np.pi*r/Rc) / (np.pi*r/Rc)


class XRD():
    def __init__(self, frname):
        with open(frname) as f:
            lines = np.array(f.readlines(), dtype=object)
        self.atoms = int(lines[3])
        self.siatoms = int(self.atoms/3)
        self.oatoms = int(self.siatoms*2)
        lines = lines.reshape(-1, self.atoms+9)
        self.L = lines[:, 5:8]
        self.steps = int(self.L.shape[0])
        data = lines[:, 9:]
        self.L = np.array(" ".join(list(self.L.reshape(-1))).split(),
                          dtype=float).reshape(self.steps, 3, 2)
        data = np.array(" ".join(list(data.reshape(-1))).split(),
                        dtype=object).reshape(self.steps, self.atoms, -1)
        data[:, :, :2] = data[:, :, :2].astype(int)
        data[:, :, 3:] = data[:, :, 3:].astype(float)

        si_idx = (data[0, :, 1] == 1)  # o_idx = ~si_idx
        ss_idx = si_idx * si_idx.reshape(-1, 1)
        so_idx = si_idx * ~si_idx.reshape(-1, 1)
        oo_idx = ~si_idx * ~si_idx.reshape(-1, 1)
        # 三角行列にする
        ss_idx *= np.tri(self.atoms, k=-1, dtype=bool)
        # so_idx *= si_idx
        oo_idx *= np.tri(self.atoms, k=-1, dtype=bool)

        fss = np.empty((Ttheta.shape[0]))
        fso = np.empty((Ttheta.shape[0]))
        foo = np.empty((Ttheta.shape[0]))
        # 原子散乱因子を fss, fso, foo にそれぞれ設定
        for idx, i in enumerate(Ttheta):
            fss[idx] = fsi[idx] ** 2
            fso[idx] = fsi[idx] * fo[idx]
            foo[idx] = fo[idx] ** 2

        start, end = 0, self.steps
        start, end = 0, 1       # only onesterp
        loops = end - start
        self.Arc = np.zeros((loops, int(Ttheta.shape[0])))
        self.Brc = np.zeros((loops, int(Ttheta.shape[0])))
        self.Crc = np.zeros((loops, int(Ttheta.shape[0])))
        self.Drc = np.zeros((loops, int(Ttheta.shape[0])))
        self.nolorch = np.zeros((loops, int(Ttheta.shape[0])))
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        # SS = np.zeros((int(Ttheta.shape[0])))
        # SO = np.zeros((int(Ttheta.shape[0])))
        # OO = np.zeros((int(Ttheta.shape[0])))

        for IDX, i in enumerate(range(start, end)):
            if i % 10 == 0:
                print("{} progressing...".format(i))
            lx = self.L[i, 0, 1] - self.L[i, 0, 0]
            dx = data[i, :, -3] - data[i, :, -3].reshape(-1, 1)
            dy = data[i, :, -2] - data[i, :, -2].reshape(-1, 1)
            dz = data[i, :, -1] - data[i, :, -1].reshape(-1, 1)
            dx = dx - lx * np.round((dx/lx).astype(float))
            dy = dy - lx * np.round((dy/lx).astype(float))
            dz = dz - lx * np.round((dz/lx).astype(float))
            dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
            drss = dr[ss_idx]
            drso = dr[so_idx]
            droo = dr[oo_idx]
            # 体格成分 0 消去
            drss = drss[np.nonzero(drss)]
            drso = drso[np.nonzero(drso)]
            droo = droo[np.nonzero(droo)]

            for idx, q in enumerate(Q):
                # si-si XRD
                x = q * drss
                Ac = Lorch(drss, lattice*np.sqrt(3))
                Bc = Lorch(drss, lattice)
                Cc = Lorch(drss, lattice*np.sqrt(3)*0.5)
                Dc = Lorch(drss, lattice*0.5)
                Lc = Lorch(drss, 20)
                self.Arc[IDX, idx] += np.sum(np.sin(x)/x * Ac) * fss[idx] * 2
                self.Brc[IDX, idx] += np.sum(np.sin(x)/x * Bc) * fss[idx] * 2
                self.Crc[IDX, idx] += np.sum(np.sin(x)/x * Cc) * fss[idx] * 2
                self.Drc[IDX, idx] += np.sum(np.sin(x)/x * Dc) * fss[idx] * 2
                self.nolorch[IDX, idx] += np.sum(np.sin(x)/x * fss[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * fss[idx] * 2)
                # si-o XRD
                x = q * drso
                Ac = Lorch(drso, lattice*np.sqrt(3))
                Bc = Lorch(drso, lattice)
                Cc = Lorch(drso, lattice*np.sqrt(3)*0.5)
                Dc = Lorch(drso, lattice*0.5)
                Lc = Lorch(drso, 20)
                self.Arc[IDX, idx] += np.sum(np.sin(x)/x * Ac) * fso[idx] * 2
                self.Brc[IDX, idx] += np.sum(np.sin(x)/x * Bc) * fso[idx] * 2
                self.Crc[IDX, idx] += np.sum(np.sin(x)/x * Cc) * fso[idx] * 2
                self.Drc[IDX, idx] += np.sum(np.sin(x)/x * Dc) * fso[idx] * 2
                self.nolorch[IDX, idx] += np.sum(np.sin(x)/x * fso[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * fso[idx] * 2)
                # o-o XRD
                x = q * droo
                Ac = Lorch(droo, lattice*np.sqrt(3))
                Bc = Lorch(droo, lattice)
                Cc = Lorch(droo, lattice*np.sqrt(3)*0.5)
                Dc = Lorch(droo, lattice*0.5)
                Lc = Lorch(droo, 20)
                self.Arc[IDX, idx] += np.sum(np.sin(x)/x * Ac) * foo[idx] * 2
                self.Brc[IDX, idx] += np.sum(np.sin(x)/x * Bc) * foo[idx] * 2
                self.Crc[IDX, idx] += np.sum(np.sin(x)/x * Cc) * foo[idx] * 2
                self.Drc[IDX, idx] += np.sum(np.sin(x)/x * Dc) * foo[idx] * 2
                self.nolorch[IDX, idx] += np.sum(np.sin(x)/x * foo[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * foo[idx] * 2)
                # sum coefficient
                self.Arc[IDX, idx] += fss[idx] * self.siatoms
                self.Arc[IDX, idx] += foo[idx] * self.oatoms
                self.Brc[IDX, idx] += fss[idx] * self.siatoms
                self.Brc[IDX, idx] += foo[idx] * self.oatoms
                self.Crc[IDX, idx] += fss[idx] * self.siatoms
                self.Crc[IDX, idx] += foo[idx] * self.oatoms
                self.Drc[IDX, idx] += fss[idx] * self.siatoms
                self.Drc[IDX, idx] += foo[idx] * self.oatoms
                self.nolorch[IDX, idx] += fss[idx] * self.siatoms
                self.nolorch[IDX, idx] += foo[idx] * self.oatoms
                self.lorch[IDX, idx] += fss[idx] * self.siatoms
                self.lorch[IDX, idx] += foo[idx] * self.oatoms

        Arc_ave = np.average(self.Arc, axis=0)/self.atoms
        Brc_ave = np.average(self.Brc, axis=0)/self.atoms
        Crc_ave = np.average(self.Crc, axis=0)/self.atoms
        Drc_ave = np.average(self.Drc, axis=0)/self.atoms
        nolorch_ave = np.average(self.nolorch, axis=0)/self.atoms
        lorch_ave = np.average(self.lorch, axis=0)/self.atoms
        ax.plot(Ttheta, Arc_ave, "b-", label="Rc = 1.73L")
        ax.plot(Ttheta, Brc_ave, "k-", label="Rc = 1.00L")
        ax.plot(Ttheta, Crc_ave, "g-", label="Rc = 0.87L")
        ax.plot(Ttheta, Drc_ave, "r-", label="Rc = 0.50L")
        ax.plot(Ttheta, nolorch_ave, "b--", label="nolorch")
        ax.plot(Ttheta, lorch_ave, "k--", label="lorch")
        ax.set_xlabel('2theta')
        ax.set_title('all xrd peak')
        plt.tight_layout()


a = XRD(args.file)
plt.legend()
if args.outfile:
    plt.savefig(args.outfile, dpi=300)
else:
    plt.show()
