#!/usr/bin/env python
import numpy as np
# import matplotlib.pyplot as plt
import argparse

par = argparse.ArgumentParser(description="bar")
par.add_argument('file', help='input file')
par.add_argument("-o", '--outfile', help='output file')
# par.add_argument('--onestep', action="store_true")
args = par.parse_args()

# ========================================
cutoff = 2.0
# ========================================

# # 図の体裁
# plt_dic = {}
# plt_dic['legend.fancybox'] = True
# plt_dic['legend.labelspacing'] = 0.3
# plt_dic['legend.numpoints'] = 3
# plt_dic['figure.figsize'] = [8, 6]
# plt_dic['axes.grid'] = True
# plt_dic['font.size'] = 12
# plt_dic['legend.fontsize'] = 12
# plt_dic['axes.labelsize'] = 14
# plt_dic['xtick.major.size'] = 5
# plt_dic['xtick.minor.size'] = 3
# plt_dic['xtick.direction'] = 'in'
# plt_dic['savefig.bbox'] = 'tight'
# plt_dic['savefig.dpi'] = 150
# plt_dic['savefig.transparent'] = True
# plt.rcParams.update(plt_dic)
# fig, ax = plt.subplots()

with open(args.file) as f:
    lines = np.array(f.readlines(), dtype=object)

atoms = int(lines[3])
siatoms, oatoms = int(atoms/3), int(atoms*2/3)
lines = lines.reshape(-1, atoms+9)
L = lines[:, 5:8]
steps = int(L.shape[0])
data = lines[:, 9:]
L = np.array(" ".join(list(L.reshape(-1))).split(),
             dtype=float).reshape(steps, 3, 2)
data = np.array(" ".join(list(data.reshape(-1))).split(),
                dtype=object).reshape(steps, atoms, -1)
data[:, :, :2] = data[:, :, :2].astype(int)
data[:, :, 3:] = data[:, :, 3:].astype(float)

si_idx = (data[0, :, 1] == 1)  # o_idx = ~si_idx
sidata = data[:, si_idx]
odata = data[:, ~si_idx]

# Si-O distance
Sihistogram = np.zeros(7)
Ohistogram = np.zeros(7)

for i in range(steps):
    if i % 10 == 0:
        print("{} progressing...".format(i))
    lx = L[i, 0, 1] - L[i, 0, 0]
    ly = L[i, 1, 1] - L[i, 1, 0]
    lz = L[i, 2, 1] - L[i, 2, 0]
    # Si 中心に O は何配位している？
    dx = odata[i, :, 3] - sidata[i, :, 3].reshape(-1, 1)
    dy = odata[i, :, 4] - sidata[i, :, 4].reshape(-1, 1)
    dz = odata[i, :, 5] - sidata[i, :, 5].reshape(-1, 1)
    dx = dx - lx * np.round((dx/lx).astype(float))
    dy = dy - lx * np.round((dy/ly).astype(float))
    dz = dz - lx * np.round((dz/lz).astype(float))
    dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
    flag = dr < cutoff
    # flag = dr < 3.0
    flag = np.sum(flag, axis=1)
    hist, bins = np.histogram(flag, bins=7, range=(0.0, 7.0))
    Sihistogram += hist
    # O 中心に Si は何配位している？
    dx = sidata[i, :, 3] - odata[i, :, 3].reshape(-1, 1)
    dy = sidata[i, :, 4] - odata[i, :, 4].reshape(-1, 1)
    dz = sidata[i, :, 5] - odata[i, :, 5].reshape(-1, 1)
    dx = dx - lx * np.round((dx/lx).astype(float))
    dy = dy - lx * np.round((dy/ly).astype(float))
    dz = dz - lx * np.round((dz/lz).astype(float))
    dr = np.sqrt((dx**2 + dy**2 + dz**2).astype(float))
    flag = dr < cutoff
    # flag = dr < 3.0
    flag = np.sum(flag, axis=1)
    hist, bins = np.histogram(flag, bins=7, range=(0.0, 7.0))
    Ohistogram += hist

SH = Sihistogram/np.sum(Sihistogram)*100
OH = Ohistogram/np.sum(Ohistogram)*100
print("|      |{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|".format(
    *list(bins.astype(int))))
print("|Si中心|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|".format(
    *list(SH)))
print("| O中心|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|".format(
    *list(OH)))
if args.outfile:
    w = open(args.outfile+".txt", "w")
    w.write("|      |{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|{:d}配位|\n".format(*list(bins.astype(int))))
    w.write("|Si中心|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|\n".format(*list(SH)))
    w.write("| O中心|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|{:5.2f}|\n".format(*list(OH)))
