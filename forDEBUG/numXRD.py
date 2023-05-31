#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file', nargs="+")
mark = ["-r", "-g", "-b", "-c", "-m", "-y"]
args = par.parse_args()

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

# const
lam = 1.5406
Ttheta = np.linspace(15, 50, 351)
Q = 4 * np.pi * np.sin(Ttheta / 360 * np.pi) / lam
Qmax = 4 * np.pi / lam

fsi = np.zeros_like(Ttheta)
# print("atomic scattering coeff, scattering factor, 2theta")
for idx, q in enumerate(Q):  # atom scattering factors
    Ar = (q / 4 / np.pi)**2  # Ar = (q / 4pi) ** 2
    fsi[idx] = (5.275329*np.exp(-2.631338*Ar) + 3.191038*np.exp(-33.73073*Ar) +
                1.511514*np.exp(-0.081119*Ar) + 1.356849*np.exp(-86.28864*Ar) +
                2.519114*np.exp(-1.170087*Ar) + 0.145073)


def Lorch(r, Rc):  # Window Func.
    return np.sin(np.pi*r/Rc) / (np.pi*r/Rc)


class XRD():
    def __init__(self, number, frname):
        with open(frname) as f:
            lines = np.array(f.readlines(), dtype=object)
        self.base = frname.split(".")[0]
        # ITEM: ATOMS id type element x y z vx vy vz
        # -> Xin = 5 - 2 = 3
        self.Xin = lines[8].split().index("x") - 2
        self.atoms = int(lines[3])
        self.siatoms = int(self.atoms/3)
        lines = lines.reshape(-1, self.atoms+9)
        self.L = lines[:, 5:8]
        self.steps = int(self.L.shape[0])
        data = lines[:, 9:]
        self.L = np.array(" ".join(list(self.L.reshape(-1))).split(),
                          dtype=float).reshape(self.steps, 3, 2)
        data = np.array(" ".join(list(data.reshape(-1))).split(),
                        dtype=object).reshape(self.steps, self.atoms, -1)
        si_mask = data[0, :, 1] == "1"
        data = data[:, si_mask, :]  # Si 原子だけのデータに変える

        data[:, :, :2] = data[:, :, :2].astype(int)
        data[:, :, 3:] = data[:, :, 3:].astype(float)
        start, end = 0, self.steps

        # start, end = 0, 1
        loops = end - start
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        for IDX, i in enumerate(range(start, end)):
            if i % 100 == 0:
                print("{} progressing...".format(i))
            lx = self.L[i, 0, 1] - self.L[i, 0, 0]
            dx = data[i, :, self.Xin] - data[i, :, self.Xin].reshape(-1, 1)
            dy = data[i, :, self.Xin+1] - data[i, :, self.Xin+1].reshape(-1, 1)
            dz = data[i, :, self.Xin+2] - data[i, :, self.Xin+2].reshape(-1, 1)
            dx = dx - lx * np.round((dx/lx).astype(float))
            dy = dy - lx * np.round((dy/lx).astype(float))
            dz = dz - lx * np.round((dz/lx).astype(float))
            dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
            # http://www.beam2d.net/blog/2021/06/27/offdiagonal/
            # 0となっている対角成分を除去
            dr = dr.reshape(self.siatoms ** 2)[:-1].reshape(
                self.siatoms - 1, self.siatoms + 1)[:, 1:].reshape(
                self.siatoms, self.siatoms - 1).reshape(-1)
            # Rc = lx * np.sqrt(3)
            # Rc = lx * np.sqrt(3) / 2
            Rc = 25
            # カットオフ Rc 以上の値を除去
            self.rc = np.copy(dr.reshape(self.siatoms * (self.siatoms-1)))
            np.place(self.rc, self.rc > Rc, 0)  # dr upper limit: Rc = lx / 2
            self.rc = self.rc[np.nonzero(self.rc)]
            Lc = Lorch(self.rc, Rc)
            # Lc = Lorch(dr, Rc)
            for idx, q in enumerate(Q):
                x = q * self.rc
                coeff = fsi[idx] * fsi[idx] / (self.atoms/3)
                self.lorch[IDX, idx] += np.sum(np.sin(x) / x * Lc) * coeff
                self.lorch[IDX, idx] += coeff * self.atoms/3

        lorch_ave = np.average(self.lorch, axis=0)
        # lorch_std = np.std(self.lorch, axis=0)

        # label = frname.split("_")[0]
        #ax.plot(Ttheta, lorch_ave, "-", label="3**0.5 * L/2")
        ax.plot(Ttheta, lorch_ave, mark[number], label=self.base)
        # ax.fill_between(Ttheta, lorch_ave-lorch_std, lorch_ave+lorch_std)
        # ax.plot(Ttheta, fpsi, "-", label="factor")
        ax.set_xlabel('2theta')
        ax.set_ylabel('Intensity')
        ax.legend()


# mark = ["r-", "b-", "g-"]
for idx, name in enumerate(args.file):
    a = XRD(idx, name)
# diff = a.lorch[0, :-1] - a.lorch[0, 1:]  # 増加量
# diff = diff > 0
# print("\nlorch\n")
# for i in range(diff.shape[0] - 1):
#     if diff[i] != diff[i+1]:
#         print(Ttheta[i], "2 theta")
# diff = a.foo[0, :-1] - a.foo[0, 1:]  # 増加量
# diff = diff > 0
# print("\nfoo\n")
# for i in range(diff.shape[0] - 1):
#     if diff[i] != diff[i+1]:
#         print(Ttheta[i], "2 theta")
plt.tight_layout()
# ax.set_ylim(-20, 200)
plt.savefig("numXRD.png")
