#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('--onestep', action="store_true")
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
fig, ax = plt.subplots(nrows=2)

# 原子散乱因子と、散乱定数Q を設定
lam = 1.5406
Ttheta = np.linspace(15, 80, 651)
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
        self.ss = np.zeros_like(Ttheta)
        self.so = np.zeros_like(Ttheta)
        self.oo = np.zeros_like(Ttheta)
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
        
        # make coefficient 3dim array
        si_idx = (data[0, :, 1] == 1)  # o_idx = ~si_idx
        ss_idx = si_idx * si_idx.reshape(-1, 1)
        so_idx = si_idx * ~si_idx.reshape(-1, 1)
        oo_idx = ~si_idx * ~si_idx.reshape(-1, 1)
        # print(np.sum(ss_idx),np.sum(so_idx),np.sum(oo_idx))
        # 三角行列にする
        ss_idx *= np.tri(self.atoms, k=-1, dtype=bool)
        so_idx *= si_idx
        oo_idx *= np.tri(self.atoms, k=-1, dtype=bool)
        # print(np.sum(ss_idx),np.sum(so_idx),np.sum(oo_idx))
        # exit()

        fss = np.empty((Ttheta.shape[0]))
        fso = np.empty((Ttheta.shape[0]))
        foo = np.empty((Ttheta.shape[0]))
        # 原子散乱因子を fss, fso, foo にそれぞれ設定
        for idx, i in enumerate(Ttheta):
            fss[idx] = fsi[idx] ** 2
            fso[idx] = fsi[idx] * fo[idx]
            foo[idx] = fo[idx] ** 2

        start, end = 0, self.steps
        start, end = 0, 1
        loops = end - start
        # self.foo = np.zeros((loops, int(Ttheta.shape[0])))
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        SS = np.zeros((int(Ttheta.shape[0])))
        SO = np.zeros((int(Ttheta.shape[0])))
        OO = np.zeros((int(Ttheta.shape[0])))
        # print(np.sum(ss_idx), np.sum(so_idx), np.sum(oo_idx))

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
            Rc = 15
            # カットオフ Rc 以上の値を除去
            np.place(drss, drss > Rc, 0)  # dr upper limit: Rc
            np.place(drso, drso > Rc, 0)  # dr upper limit: Rc
            np.place(droo, droo > Rc, 0)  # dr upper limit: Rc
            drss = drss[np.nonzero(drss)]
            drso = drso[np.nonzero(drso)]
            droo = droo[np.nonzero(droo)]

            for idx, q in enumerate(Q):
                x = q * drss        # Qr(Si-Si)
                Lc = Lorch(drss, Rc)
                SS[idx] += np.sum(np.sin(x)/x * Lc * fss[idx] * 2)/self.siatoms
                # self.foo[IDX, idx] += np.sum(np.sin(x)/x * fss[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * fss[idx] * 2)
                x = q * droo        # Qr(O-O)
                Lc = Lorch(droo, Rc)
                OO[idx] += np.sum(np.sin(x)/x * Lc * foo[idx] * 2)/self.siatoms
                # self.foo[IDX, idx] += np.sum(np.sin(x)/x * foo[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * foo[idx] * 2)
                x = q * drso        # Qr(Si-O)
                Lc = Lorch(drso, Rc)
                SO[idx] += np.sum(np.sin(x)/x * Lc * fso[idx] * 2)/self.siatoms
                # self.foo[IDX, idx] += np.sum(np.sin(x)/x * fso[idx] * 2)
                self.lorch[IDX, idx] += np.sum(np.sin(x)/x * Lc * fso[idx] * 2)
                # self.foo[IDX, idx] += fss[idx] * self.siatoms
                # self.foo[IDX, idx] += foo[idx] * self.oatoms
                self.lorch[IDX, idx] += fss[idx] * self.siatoms
                self.lorch[IDX, idx] += foo[idx] * self.oatoms
            if args.onestep:
                ax.plot(Ttheta, SS, "b-", label="S-S")
                ax.plot(Ttheta, SO, "g-", label="S-O")
                ax.plot(Ttheta, OO, "-", color="black", label="O-O")
                ax.legend()
                plt.show()
                exit("onestep end")

        # foo_ave = np.average(self.foo, axis=0)/self.atoms
        lorch_ave = np.average(self.lorch, axis=0)/self.atoms
        lorch_std = np.std(self.lorch, axis=0)/self.atoms
        # print(lorch_std)
        # label = frname.split("_")[0]
        SS = SS/loops + fss[idx] / 3
        OO = OO/loops + foo[idx] / 3 * 2
        ax[0].plot(Ttheta, SS, "b-", label="S-S")
        ax[0].plot(Ttheta, SO, "g-", label="S-O")
        ax[0].plot(Ttheta, OO, "-", color="black", label="O-O")
        # ax.plot(Ttheta, foo_ave, "--", label=label)
        ax[1].plot(Ttheta, lorch_ave, "r-", label="3**0.5 * L/2")
        ax[1].fill_between(Ttheta,
                           lorch_ave+lorch_std, lorch_ave-lorch_std, alpha=0.1)
        ax[0].set_xlabel('2theta')
        ax[1].set_xlabel('2theta')
        ax[0].set_title('each atom pairs')
        ax[1].set_title('all xrd peak')
        ax[0].set_ylabel('Intensity')
        ax[0].legend()


a = XRD(args.file)
# b = XRD("../glass/"+args.file)
# data = np.loadtxt("COLVARcd")[:, 1:]
# x = np.linspace(15, 50, data.shape[1])
# ax.plot(x, np.mean(data, axis=0), "ro-")
# data = np.loadtxt("../glass/COLVARcd")[:, 1:]
# x = np.linspace(15, 50, data.shape[1])
# ax.plot(x, np.mean(data, axis=0), "bo-")

plt.show()
