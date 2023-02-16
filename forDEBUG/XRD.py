#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
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
# Ttheta = np.linspace(15, 80, 651)
Ttheta = np.linspace(15, 50, 351)
Q = 4 * np.pi * np.sin(Ttheta / 360 * np.pi) / lam

fsi = np.zeros_like(Ttheta)
# print("atomic scattering coeff, scattering factor, 2theta")
for idx, q in enumerate(Q):  # atom scattering factors
    Ar = (q / 4 / np.pi)**2  # Ar = (q / 4pi) ** 2
    fsi[idx] = (5.275329*np.exp(-2.631338*Ar) + 3.191038*np.exp(-33.73073*Ar) +
                1.511514*np.exp(-0.081119*Ar) + 1.356849*np.exp(-86.28864*Ar) +
                2.519114*np.exp(-1.170087*Ar) + 0.145073)
#     if idx % 10 == 0:
#         print(fsi[idx]**2, q, Ttheta[idx])
# exit()


def Lorch(r, Rc):  # Window Func.
    return np.sin(np.pi*r/Rc) / (np.pi*r/Rc)


class XRD():
    def __init__(self, frname):
        with open(frname) as f:
            lines = np.array(f.readlines(), dtype=object)
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
        sidata = data[:, si_mask, :]
        sidata[:, :, :2] = sidata[:, :, :2].astype(int)
        sidata[:, :, 3:] = sidata[:, :, 3:].astype(float)
        start, end = 0, self.steps

        # start, end = 0, 1
        # start, end = 101, 201
        loops = end - start
        self.foo = np.zeros((loops, int(Ttheta.shape[0])))
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        self.morch = np.zeros((loops, int(Ttheta.shape[0])))
        # RANGE = np.linspace(900, 999, 100, dtype=int)
        # for IDX, i in enumerate(RANGE):
        for IDX, i in enumerate(range(start, end)):
            if i % 100 == 0:
                print("{} progressing...".format(i))
            lx = self.L[i, 0, 1] - self.L[i, 0, 0]
            dx = sidata[i, :, -3] - sidata[i, :, -3].reshape(-1, 1)
            dy = sidata[i, :, -2] - sidata[i, :, -2].reshape(-1, 1)
            dz = sidata[i, :, -1] - sidata[i, :, -1].reshape(-1, 1)
            dx = dx - lx * np.round((dx/lx).astype(float))
            dy = dy - lx * np.round((dy/lx).astype(float))
            dz = dz - lx * np.round((dz/lx).astype(float))
            dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
            # http://www.beam2d.net/blog/2021/06/27/offdiagonal/
            # 0となっている対角成分を除去
            dr = dr.reshape(self.siatoms ** 2)[:-1].reshape(
                self.siatoms - 1, self.siatoms + 1)[:, 1:].reshape(
                self.siatoms, self.siatoms - 1).reshape(-1)
            Rc = lx * np.sqrt(3)
            Rc = 25  # 2022: 25, 2018: 18
            # カットオフ Rc 以上の値を除去
            self.rc = np.copy(dr.reshape(self.siatoms * (self.siatoms-1)))
            np.place(self.rc, self.rc > Rc, 0)  # dr upper limit: Rc = lx / 2
            self.rc = self.rc[np.nonzero(self.rc)]
            Lc = Lorch(self.rc, Rc)
            # Lc = Lorch(dr, Rc)
            for idx, q in enumerate(Q):
                # X = q * dr
                x = q * self.rc
                coeff = fsi[idx] * fsi[idx] / (self.atoms/3)
                # Lorch half of calc cell
                self.morch[IDX, idx] += np.sum(np.sin(x) / x * Lc * coeff)
                self.morch[IDX, idx] += fsi[idx]**2
                # self.morch[IDX, idx] += np.sum(np.sin(X) / X *
                #                                Lorch(dr, lx*np.sqrt(3))) * coeff
                # non_Lorch
                self.foo[IDX, idx] += np.sum(np.sin(x) / x) * coeff
                self.foo[IDX, idx] += fsi[idx]**2

        foo_ave = np.average(self.foo, axis=0)
        morch_ave = np.average(self.morch, axis=0)
        foo_std = np.std(self.foo, axis=0)
        # lorch_std = np.std(self.lorch, axis=0)

        # label = frname.split("_")[0]
        ax.plot(Ttheta, foo_ave, "b--", label="hoge")
        ax.plot(Ttheta, morch_ave, "g-", label="glass")
        #ax.fill_between(Ttheta, foo_ave-foo_std, foo_ave+foo_std)
        # ax.fill_between(Ttheta, lorch_ave-lorch_std, lorch_ave+lorch_std)
        # ax.plot(Ttheta, fpsi, "-", label="factor")
        ax.set_xlabel('2theta')
        ax.set_ylabel('Intensity')
        ax.legend()


# a = XRD("crystal_2ns.lammpstrj")
a = XRD(args.file)
diff = a.lorch[0, :-1] - a.lorch[0, 1:]  # 増加量
diff = diff > 0
print("\nlorch\n")
for i in range(diff.shape[0] - 1):
    if diff[i] != diff[i+1]:
        print(Ttheta[i], "2 theta")
diff = a.foo[0, :-1] - a.foo[0, 1:]  # 増加量
diff = diff > 0
print("\nfoo\n")
for i in range(diff.shape[0] - 1):
    if diff[i] != diff[i+1]:
        print(Ttheta[i], "2 theta")
plt.tight_layout()
plt.text(40, 0, "(3/2)**0.5 * L")
# ax.set_ylim(-50, 300)
plt.show()
