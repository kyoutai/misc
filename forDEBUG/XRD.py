#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument('-e', '--elem', choices=["Si", "Na"], default="Si")
par.add_argument('-r', '--rc', help='Lorch cutoff, Rc', default=25, type=int)
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
Ttheta = np.linspace(15, 80, 651)
# Ttheta = np.linspace(15, 50, 351)
Q = 4 * np.pi * np.sin(Ttheta / 360 * np.pi) / lam
hoge = np.sin(Ttheta/360*np.pi)/lam

fsi = np.zeros_like(Ttheta)
# print("atomic scattering coeff, scattering factor, 2theta")
if args.elem == "Si":
    for idx, q in enumerate(Q):  # atom scattering factors
        A = (q / 4 / np.pi)**2  # Ar = (q / 4pi) ** 2
        fsi[idx] = (5.275329*np.exp(-2.631338*A)+3.191038*np.exp(-33.73073*A) +
                    1.511514*np.exp(-0.081119*A)+1.356849*np.exp(-86.28864*A) +
                    2.519114*np.exp(-1.170087*A)+0.145073)
elif args.elem == "Na":
    for idx, q in enumerate(Q):
        A = (q / 4 / np.pi)**2
        # fsi[idx]=(4.910127*np.exp(-3.281434*A)+3.081783*np.exp(-9.119178*A) +
        #           1.262067*np.exp(-0.102763*A)+1.098938*np.exp(-132.0139*A) +
        #           0.560991*np.exp(-0.405878*A)+0.0797120)
        fsi[idx] = (4.7626*np.exp(-3.2550*A) + 3.1736*np.exp(-8.84220*A) +
                    1.2674*np.exp(-0.3136*A) + 1.1128*np.exp(-129.424*A) +
                    0.676)
        # International Tables for Crystallography: Mathematical, Physical and
        # Chemical Tables, ed Prince E (Springer, Berlin), Vol C, pp 554–595.
else:
    exit("wrong element! check your argument: -e or --elem")


def Lorch(r, Rc):  # Window Func.
    x = np.pi * r / Rc
    return np.sin(x)/x


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
        if args.elem == "Si":
            self.si_XRD(data)
        elif args.elem == "Na":
            self.na_XRD(data)

        # plot XRD peak
        morch_ave = np.average(self.morch, axis=0)
        # foo_ave = np.average(self.foo, axis=0)
        # foo_std = np.std(self.foo, axis=0)
        # ax.plot(Ttheta, foo_ave, "b--", label="hoge")
        ax.plot(Ttheta, morch_ave, "r-", label="crystal")
        # ax.fill_between(Ttheta, foo_ave-foo_std, foo_ave+foo_std)
        ax.set_xlabel('2theta')
        ax.set_ylabel('Intensity')
        # ax.legend(loc="lower right")
        ax.legend(loc="best")

    def si_XRD(self, data):
        si_mask = data[0, :, 1] == "1"
        sidata = data[:, si_mask, :]
        sidata[:, :, :2] = sidata[:, :, :2].astype(int)
        sidata[:, :, 3:] = sidata[:, :, 3:].astype(float)
        start, end = 0, self.steps
        # start, end = 101, 201
        loops = end - start
        self.foo = np.zeros((loops, int(Ttheta.shape[0])))
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        self.morch = np.zeros((loops, int(Ttheta.shape[0])))
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
            Rc = 25  # 2022: 25, 2018: 18
            # カットオフ Rc 以上の値を除去
            self.rc = np.copy(dr.reshape(self.siatoms * (self.siatoms-1)))
            np.place(self.rc, self.rc > Rc, 0)  # dr upper limit: Rc = lx / 2
            self.rc = self.rc[np.nonzero(self.rc)]
            Lc = Lorch(self.rc, Rc)
            for idx, q in enumerate(Q):
                x = q * self.rc
                coeff = fsi[idx] * fsi[idx] / (self.atoms/3)
                # Lorch half of calc cell
                self.morch[IDX, idx] += np.sum(np.sin(x) / x * Lc * coeff)
                self.morch[IDX, idx] += fsi[idx]**2
                self.foo[IDX, idx] += np.sum(np.sin(x) / x) * coeff
                self.foo[IDX, idx] += fsi[idx]**2

    def na_XRD(self, data):
        data[:, :, :2] = data[:, :, :2].astype(int)
        data[:, :, 3:] = data[:, :, 3:].astype(float)
        start, end = 0, self.steps
        start, end = 0, 1
        loops = end - start
        self.foo = np.zeros((loops, int(Ttheta.shape[0])))
        self.lorch = np.zeros((loops, int(Ttheta.shape[0])))
        self.morch = np.zeros((loops, int(Ttheta.shape[0])))
        for IDX, i in enumerate(range(start, end)):
            if i % 100 == 0:
                print("{} progressing...".format(i))
            lx = self.L[i, 0, 1] - self.L[i, 0, 0]
            dx = data[i, :, -3] - data[i, :, -3].reshape(-1, 1)
            dy = data[i, :, -2] - data[i, :, -2].reshape(-1, 1)
            dz = data[i, :, -1] - data[i, :, -1].reshape(-1, 1)
            dx = dx - lx * np.round((dx/lx).astype(float))
            dy = dy - lx * np.round((dy/lx).astype(float))
            dz = dz - lx * np.round((dz/lx).astype(float))
            dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
            # http://www.beam2d.net/blog/2021/06/27/offdiagonal/
            # 0となっている対角成分を除去
            dr = dr.reshape(self.atoms ** 2)[:-1].reshape(
                self.atoms - 1, self.atoms + 1)[:, 1:].reshape(
                self.atoms, self.atoms - 1).reshape(-1)
            # Lc = Lorch(dr, 11)
            # ax.plot(dr*np.pi/11, Lc, ".")
            # plt.show()
            # exit()
            # カットオフ Rc 以上の値を除去
            self.rc = np.copy(dr.reshape(self.atoms * (self.atoms-1)))
            np.place(self.rc, self.rc > args.rc, 0)  # dr upper limit: Rc
            self.rc = self.rc[np.nonzero(self.rc)]
            Lc = Lorch(self.rc, args.rc)

            for idx, q in enumerate(Q):
                x, coeff = q * self.rc, fsi[idx] ** 2 / self.atoms
                # Lorch half of calc cell
                # different id pairs -> sin(pi*r/Rc)/(pi*r/Rc)*W(r)
                self.morch[IDX, idx] += np.sum(np.sin(x) / x * Lc) * coeff
                # same id pairs -> r=0, sin(pi*r/Rc)/(pi*r/Rc)*W(r) = 1
                self.morch[IDX, idx] += fsi[idx]**2
                # self.foo[IDX, idx] += np.sum(np.sin(x) / x) * coeff
                # self.foo[IDX, idx] += fsi[idx]**2


a = XRD(args.file)
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
# plt.text(40, 0, "(3/2)**0.5 * L")
ax.set_ylim(-50, 300)
ax.set_xlim(20, 80)
plt.show()
