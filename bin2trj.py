#!/usr/bin/env python
import argparse
import numpy as np
import os.path
import re
import struct
from collections import OrderedDict

"""
albiteなどでbinary2txtがうまく動かないときがあるので
少し遅いがpythonでtrjファイルに書き出すプログラム
"""

description = """describe"""
par = argparse.ArgumentParser()
par = argparse.ArgumentParser(description=description)
par.add_argument('infiles', nargs="+")
par.add_argument('-d', '--data')
par.add_argument('-o', '--out', help="lammpstrj file BASE name")
par.add_argument('--bin_new', action='store_true')
par.add_argument('--bertrand', action='store_true')
args = par.parse_args()

# label: [mass charge charge(bertrand) elem]
FF = OrderedDict({})
FF["Si"] = [28.0855, 2.4, 1.89, "Si"]
FF["O"] = [15.9994, -1.2, -0.945, "O"]
FF["Na"] = [22.9897, 0.6, 0.4725, "Na"]
FF["K"] = [39.0983, 0.6, 0.4725, "K"]
FF["Al"] = [26.9815, 1.8, 1.4175, "Al"]
FF["Ca"] = [40.078, 1.2, 0.945, "Ca"]
FF["Mg"] = [24.306, 1.2, 0.945, "Mg"]
FF["Fe"] = [55.845, 1.2, 0.945, "Fe"]
FF["Ti"] = [47.867, 2.4, 1.89, "Ti"]
FF["C"] = [12.0107, 2.4, 1.89, "C"]
FF["Si_irr"] = [28.0855, 2.4, 1.89, "Si"]
FF["O_irr"] = [15.9994,  -1.2, -0.945, "O"]
FF["Na_irr"] = [22.9897, 0.6, 0.4725, "Na"]
FF["K_irr"] = [39.0983, 0.6, 0.4725, "K"]
FF["Al_irr"] = [26.9815, 1.8, 1.4175, "Al"]
FF["Ca_irr"] = [40.078, 1.2, 0.945, "Ca"]
FF["Mg_irr"] = [24.306, 1.2, 0.945, "Mg"]
FF["Fe_irr"] = [55.845, 1.2, 0.945, "Fe"]
FF["Ti_irr"] = [47.867, 2.4, 1.89, "Ti"]
FF["C_irr"] = [12.0107, 2.4, 1.89, "C"]

class LammpsData():
    def __init__(self):
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 0
        self.beta = 0
        self.gamma = 0
        self.box = np.empty((3, 2))  # xlo xhi ylo yhi zlo zhi
        self.box_tri = None         # only triclinic
        self.masses = np.empty((0, 3))  # 0:number 1:mass 2:elem 3:label
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        self.atoms = np.empty((0, 11))
        self.knock_idx = None  # knock-on atomのindex
        self.system = "tri"
        self.step = 0

    def Read_init(self, infile):
        self.Bin2Data(args.data, infile, 0)
        self.WriteLammpsTrj(infile)
        print(self.step)
        for i in range(self.steps):
            if i != 0:
                self.Bin2Data(args.data, infile, i)
                self.WriteLammpsTrj(infile)
                print(self.step)

    def CalcMatrix(self, a, b, c, alpha, beta, gamma):
        v1 = [a, 0, 0]
        v2 = [b * np.cos(gamma),  b * np.sin(gamma), 0]
        v3 = [c * np.cos(beta),
              c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
              c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
                          - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2) /
              np.sin(gamma)]
        M = np.array([v1, v2, v3])
        return M

    def GetBinData(self, offset, f):
        f.seek(offset)
        head_size = 16 + 4 + 4*6 + 8*6 + 4 + 4
        #step, atoms = struct.unpack("ll", f.read(16)) # step, atoms
        self.step = struct.unpack("l", f.read(8))[0] # step, atoms
        atoms = struct.unpack("l", f.read(8))[0] # step, atoms
        triclnic = struct.unpack("i", f.read(4))[0]   # triclinic -> 1
        v = struct.unpack("6i", f.read(4*6))  # boundary
        box = struct.unpack("6d", f.read(8*6)) # xlo, xhi, ylo, yhi, zlo, zhi
        M = np.array(box).reshape((3, 2))
        if triclnic == 1:
            box_tri = struct.unpack("3d", f.read(8*3))
            M = np.hstack((M, np.array(box_tri).reshape((-1, 1))))
            head_size += 8*3
        size_one = struct.unpack("i", f.read(4))[0]  # 1行のデータ数
        nchunk = struct.unpack("i", f.read(4))[0]
        p = 0
        data = np.empty((atoms, size_one))
        data_size = 0
        for i in range(nchunk):
            n = struct.unpack("i", f.read(4))[0]
            m = int(n/size_one)
            d = struct.unpack("{:.0f}d".format(n), f.read(n * 8))
            data_size += n*8 + 4
            d = np.array(d).reshape((m, size_one))
            data[p: p+m, :] = d
            p = p + m
        return head_size + data_size, M, data

    def GetBinData2(self, offset, f):
        f.seek(offset)
        # head_size = 16 + 4 + 4*6 + 8*6 + 4 + 4
        #step, atoms = struct.unpack("ll", f.read(16)) # step, atoms
        ntimestep = struct.unpack("q", f.read(8))[0] # step, atoms
        if ntimestep < 0:
            magic_string_len = -ntimestep
            fmt = "<{}s".format(magic_string_len)
            size = magic_string_len
            magic_string = struct.unpack(fmt, f.read(magic_string_len))[0]
            magic_string = magic_string.decode()
            endian = struct.unpack("<i", f.read(4))[0]
            revision = struct.unpack("<i", f.read(4))[0]
            self.step = struct.unpack("<q", f.read(8))[0] # step
        atoms = struct.unpack("<q", f.read(8))[0] # step, atoms
        triclinic = struct.unpack("<i", f.read(4))[0]   # triclinic -> 1
        bound = struct.unpack("<6i", f.read(4*6))
        box = struct.unpack("<6d", f.read(8*6)) # xlo, xhi, ylo, yhi, zlo, zhi
        M = np.array(box).reshape((3, 2))
        if triclinic == 1:
            box_tri = struct.unpack("<3d", f.read(8*3))
            M = np.hstack((M, np.array(box_tri).reshape((-1, 1))))
        size_one = struct.unpack("<i", f.read(4))[0]  # 1行のデータ数
        if magic_string and revision > 1:
            len_ = struct.unpack("<i", f.read(4))[0]
            if len_ > 0:
                fmt = "<{}s".format(len_)
                unit_style = struct.unpack(fmt, f.read(len_))
                unit_style = unit_style.decode()
            flag = struct.unpack("<s", f.read(1))[0].decode()
            if flag == "":
                time = struct.unpack("<d", f.read(8))
            len_ = struct.unpack("<i", f.read(4))[0]
            fmt = "<{}s".format(len_)
            colmuns = struct.unpack(fmt, f.read(len_))[0].decode()
        nchunk = struct.unpack("<i", f.read(4))[0]
        p = 0
        data = np.empty((atoms, size_one))
        data_size = 0
        for i in range(nchunk):
            n = struct.unpack("<i", f.read(4))[0]
            m = int(n/size_one)
            d = struct.unpack("<{:.0f}d".format(n), f.read(n * 8))
            data_size += n*8 + 4
            d = np.array(d).reshape((m, size_one))
            data[p: p+m, :] = d
            p = p + m
        data_size = f.tell()
        return data_size, M, data

    def Bin2Data(self, data, infile, step):
        """
        binaryで保存したデータの読み込み
        dataのMassesからデータを拾ってatoms[:, 10]のlabelに入れておく
        atomsのxyzは分率座標で保存
        """
        body = open(data).read()  # data
        masses = re.search("Masses.*?\n\n(.*?)\n\n", body, re.DOTALL).group(1)
        masses = re.sub("#", "", masses)
        masses = np.array(masses.strip().split(),
                          dtype=object).reshape((-1, 4))
        masses[:, 0] = masses[:, 0].astype(int)
        masses[:, 1] = masses[:, 1].astype(float)
        # bin dataの読み込み
        filesize = os.path.getsize(infile)
        f = open(infile, "rb")
        if step == 0:
            if args.bin_new is False:
                self.size, M, data = self.GetBinData(0, f)
            else:
                self.size, M, data = self.GetBinData2(0, f)
            self.steps = int(filesize / self.size)
        else:
            # 指定ステップのデータ
            if args.bin_new is False:
                size, M, data = self.GetBinData(self.size*step, f)
            else:
                size, M, data = self.GetBinData2(self.size*step, f)
        self.item_box = M
        atoms = data.shape[0]
        if M.shape[1] == 3:  # triclinic
            xlo = M[0, 0] - np.min([0.0, M[0, 2], M[1, 2], M[0, 2] + M[1, 2]])
            xhi = M[0, 1] - np.max([0.0, M[0, 2], M[1, 2], M[0, 2] + M[1, 2]])
            ylo = M[1, 0] - np.min([0.0, M[2, 2]])
            yhi = M[1, 1] - np.max([0.0, M[2, 2]])
            zlo = M[2, 0]
            zhi = M[2, 1]
            lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
            a = lx
            b = np.sqrt(ly**2 + M[0, 2]**2)
            c = np.sqrt(lz**2 + M[1, 2]**2 + M[2, 2]**2)
            alpha = np.arccos((M[2, 2] * ly + M[0, 2] * M[1, 2]) / (b * c))
            beta = np.arccos(M[1, 2] / c)
            gamma = np.arccos(M[0, 2] / b)
            self.box = [[0.0, lx], [0.0, ly], [0.0, lz]]
            self.box_tri = M[:, 2]
            self.system = "tri"
        else:  # orthogonal
            alpha, beta, gamma = 0.5 * np.pi, 0.5 * np.pi, 0.5 * np.pi
            a = M[0, 1] - M[0, 0]
            b = M[1, 1] - M[1, 0]
            c = M[2, 1] - M[2, 0]
            self.box = M
            self.system = "aniso"
        #self.M = self.CalcMatrix(a, b, c, alpha, beta, gamma)
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        self.atoms = np.zeros((atoms, 11), dtype=object)
        self.atoms[:, 0] = data[:, 0].astype(int) # id
        self.atoms[:, 1] = data[:, 1].astype(int) # type
        self.atoms[:, 3:6] = data[:, 3:6].astype(float)  # x y z
        self.atoms[:, 6:9] = data[:, 6:9].astype(float)  # vx vy vz
        # 分率座標にしない
        #M_ = np.linalg.inv(self.M)
        #self.atoms[:, 3:6] = np.dot(self.atoms[:, 3:6].astype(float), M_)
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        # 元素名とlabelをセットする
        # masses 0:type 1:mass 2:element 3:type
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        for m in masses:
            self.atoms[self.atoms[:, 1] == m[0], 9] = m[2]
            self.atoms[self.atoms[:, 1] == m[0], 10] = m[2]

    def UpdateTypeCharge(self):
        """
        labelに基づいて、atomsのtype, charge, elemをupdateする
        """
        labels = np.unique(self.atoms[:, 10])
        self.masses = np.empty((0, 4), dtype=object)
        for k in FF.keys():
            if k in labels:
                self.masses = np.vstack(
                    (self.masses, [[0, FF[k][0], FF[k][3], k]]))
        self.masses[:, 0] = np.linspace(1, labels.shape[0], labels.shape[0])
        self.masses[:, 0] = self.masses[:, 0].astype(int)
        self.masses[:, 1] = self.masses[:, 1].astype(float)
        for m in self.masses:
            self.atoms[self.atoms[:, 10] == m[3], 1] = m[0]        # type
            if args.bertrand is True:
                self.atoms[self.atoms[:, 10] == m[3], 2] = FF[m[3]][2]  # charge bertrandのとき
            else:
                self.atoms[self.atoms[:, 10] == m[3], 2] = FF[m[3]][1]  # charge
            self.atoms[self.atoms[:, 10] == m[3], 9] = FF[m[3]][3]  # elem
        self.atoms[:, 0:2] = self.atoms[:, 0:2].astype(int)
        self.atoms[:, 2:9] = self.atoms[:, 2:9].astype(float)

    def WriteLammpsTrj(self, infile):
        """
        lammpstrj fileの生成
        """
        self.UpdateTypeCharge()
        basename = os.path.splitext(os.path.basename(infile))[0]
        if args.out:
            outfile = args.out
        else:
            outfile = basename + ".lammpstrj"
        # outfile = "tameshi.lammpstrj"
        o = open(outfile, "a")
        o.write("ITEM: TIMESTEP\n")
        o.write("{}\n".format(self.step))
        o.write("ITEM: NUMBER OF ATOMS\n")
        o.write("{}\n".format(self.atoms.shape[0]))
        o.write("ITEM: BOX BOUNDS pp pp pp\n")
        if self.item_box.shape[1] == 2:
            for i in (self.item_box):
                o.write("{} {}\n".format(*i))
        else:
            for i in (self.item_box):
                o.write("{} {} {}\n".format(*i))                
        o.write("ITEM: ATOMS id type element x y z vx vy vz\n")
        for i in (self.atoms):
            i[2] = i[-1]
            o.write("{} {} {} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\n".format(*i))

l = LammpsData()
for i in args.infiles:
    l.Read_init(i)
