#!/usr/bin/env python
import argparse
import numpy as np
import os.path
import re
import struct
from collections import OrderedDict
par = argparse.ArgumentParser()
group = par.add_mutually_exclusive_group()
group.add_argument('--ciffile', metavar=('*.cif'),
                   default=None, help="input cif file (prmitive cell)")
par.add_argument('-s', '--replicate', nargs=3, default=[1, 1, 1], metavar=1,
                 type=int,  help="supercell cell; x y z")
par.add_argument('-o', '--basename', required=True,
                 help="output basename lammps input; *.data and *.in")
args = par.parse_args()

cifbase = """
#======================================================================
# CRYSTAL DATA
#----------------------------------------------------------------------
data_VESTA_phase_1

_chemical_name_common                  'orthogonal tridymite'
_cell_length_a                         {a}
_cell_length_b                         {b}
_cell_length_c                         {c}
_cell_angle_alpha                      90.00
_cell_angle_beta                       90.00
_cell_angle_gamma                      90.00
_cell_volume                           182.256248
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


# label: [mass charge elem]
FF = OrderedDict({})
FF["Si"] = [28.0855, 2.4, "Si"]
FF["O"] = [15.9994,  -1.2, "O"]
FF["Si_irr"] = [28.0855, 2.4, "Si"]
FF["O_irr"] = [15.9994,  -1.2, "O"]

# neutron cross section
# https://www.ncnr.nist.gov/resources/n-lengths/
FF_XS = {}
FF_XS["Si"] = 2.167
FF_XS["O"] = 4.232


def CalcMatrix(a, b, c, alpha, beta, gamma):
    v1 = [a, 0, 0]
    v2 = [b * np.cos(gamma),  b * np.sin(gamma), 0]
    v3 = [c * np.cos(beta),
          c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
          c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
                      - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2) /
          np.sin(gamma)]
    M = np.array([v1, v2, v3])
    return M


class CIFDATA():
    """
    attribute: 
    lattice constant: a b c alpha beta gamma
    data: element x y z (fractional corrd)
    """

    def __init__(self, ciffile):
        self.ciffile = ciffile
        self.GetData()

    def GetPara(self, s, body):
        obj = re.search(s, body)
        v = float(obj.group(1))
        return v

    def GetData(self):
        body = open(self.ciffile).read()
        o = re.search("""_space_group_name_H-M_alt\s+'P 1'""", body)
        if o is None:
            print("{} is not primitive cell.".format(self.ciffile))
            print("Pleae use VESTA to convert P1 and export cif file.")
            exit()

        self.a = self.GetPara("_cell_length_a +([0-9\.]+)", body)
        self.b = self.GetPara("_cell_length_b +([0-9\.]+)", body)
        self.c = self.GetPara("_cell_length_c +([0-9\.]+)", body)
        alpha = self.GetPara("_cell_angle_alpha +([0-9\.]+)", body)
        beta = self.GetPara("_cell_angle_beta +([0-9\.]+)", body)
        gamma = self.GetPara("_cell_angle_gamma +([0-9\.]+)", body)
        self.alpha = alpha * np.pi / 180
        self.beta = beta * np.pi / 180
        self.gamma = gamma * np.pi / 180
        obj = re.search("_atom_site_type_symbol(.+)\nloop_\n", body, re.DOTALL)
        data = obj.group(1)
        data = re.sub("\(\d+?\)", "", data)
        data = np.array(data.split(), dtype=object).reshape(-1, 8)
        data = data[:, [7, 2, 3, 4]]  # element, x, y, z
        data[:, 1:4] = data[:, 1:4].astype(float)
        data[:, 1:4] = data[:, 1:4] - np.floor(data[:, 1:4].astype(float))
        self.data = data
        # trigonal -> orthogonal
        data[:, 1] = (data[:, 1] - 0.5 * data[:, 2]) * self.a
        data[:, 2] = data[:, 2] * self.b * np.sqrt(3) / 2
        data[:, 3] = data[:, 3] * self.c  # no change gamma
        GAWA = data
        # data[:, 1] = (data[:, 1] - 0.5 * data[:, 2]) * self.a
        # data[:, 2] = data[:, 2] * self.b * np.sqrt(3) / 2
        # data[:, 3] = data[:, 3] * self.c  # no change gamma


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

    def Cif2Data(self, ciffile):
        """
        primitiveなcifファイルを読んでlammpsdataの生成
        labelは元素記号とする。 atomsのxyzは分率座標で保存
        """
        cif = CIFDATA(ciffile)
        a, b, c = cif.a, cif.b, cif.c
        alpha, beta, gamma = cif.alpha, cif.beta, cif.gamma
        b = b * np.sin(gamma)
        self.box[:, 1] = np.array([a, b, c])
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        atoms = np.zeros((cif.data.shape[0], 11), dtype=object)
        atoms[:, 10] = None
        for i, d in enumerate(cif.data):
            atoms[i, 0] = i + 1
            atoms[i, 3:6] = d[1:4]
            atoms[i, 10] = d[0]  # label
        self.M = CalcMatrix(a, b, c, alpha, beta, gamma)
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self.atoms = atoms

    def Trj2Data(self, restart):
        """
        dataのMassesからデータを拾ってatoms[:, 10]のlabelに入れておく
        atomsのxyzは分率座標で保存
        """
        body = open(restart[0]).read()  # data
        masses = re.search("Masses.*?\n\n(.*?)\n\n", body, re.DOTALL).group(1)
        masses = re.sub("#", "", masses)
        masses = np.array(masses.strip().split(),
                          dtype=object).reshape((-1, 4))
        masses[:, 0] = masses[:, 0].astype(int)
        masses[:, 1] = masses[:, 1].astype(float)
        lines = open(restart[1]).readlines()  # trj
        atoms = int(lines[3])
        lines = lines[-(atoms + 9):]
        M = np.array([d.split() for d in lines[5:8]], dtype=float)
        if M.shape[1]:  # triclinic
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
        else:  # orthogonal
            alpha, beta, gamma = 0.5 * np.pi, 0.5 * np.pi, 0.5 * np.pi
            a = M[0, 1] - M[0, 0]
            b = M[1, 1] - M[1, 0]
            c = M[2, 1] - M[2, 0]
            self.box = M
        self.M = CalcMatrix(a, b, c, alpha, beta, gamma)
        data = np.array(" ".join(lines[9:9 + atoms]).split(), dtype=object)
        data = data.reshape((atoms, -1))
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        self.atoms = np.zeros((atoms, 11), dtype=object)
        self.atoms[:, 0] = data[:, 0].astype(int)
        self.atoms[:, 1] = data[:, 1].astype(int)
        self.atoms[:, 3:6] = data[:, 3:6].astype(float)  # x y z
        self.atoms[:, 6:9] = data[:, 6:9].astype(float)  # vx vy vz
        self.atoms[:, 10] = data[:, 2]                   # label
        M_ = np.linalg.inv(self.M)
        self.atoms[:, 3:6] = np.dot(self.atoms[:, 3:6].astype(float), M_)
        self.a, self.b, self.c = a, b, c
        self.alpha, self.beta, self.gamma = 90, 90, 90
        # 元素名とlabelをセットする
        # masses 0:type 1:mass 2:element 3:type
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        for m in masses:
            self.atoms[self.atoms[:, 1] == m[0], 9] = m[2]
            self.atoms[self.atoms[:, 1] == m[0], 10] = m[2]

    def GetBinData(self, offset, f):
        f.seek(offset)
        head_size = 16 + 4 + 4*6 + 8*6 + 4 + 4
        step, atoms = struct.unpack("ll", f.read(16))  # step, atoms
        triclnic = struct.unpack("i", f.read(4))[0]   # triclinic -> 1
        v = struct.unpack("6i", f.read(4*6))  # boundary
        box = struct.unpack("6d", f.read(8*6))  # xlo, xhi, ylo, yhi, zlo, zhi
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

    # def Bin2Data(self, restart):
    #     """
    #     binaryで保存したデータの読み込み
    #     dataのMassesからデータを拾ってatoms[:, 10]のlabelに入れておく
    #     atomsのxyzは分率座標で保存
    #     """
    #     body = open(restart[0]).read()  # data
    #     masses = re.search("Masses.*?\n\n(.*?)\n\n", body, re.DOTALL).group(1)
    #     masses = re.sub("#", "", masses)
    #     masses = np.array(masses.strip().split(),
    #                       dtype=object).reshape((-1, 4))
    #     masses[:, 0] = masses[:, 0].astype(int)
    #     masses[:, 1] = masses[:, 1].astype(float)
    #     # bin dataの読み込み
    #     filesize = os.path.getsize(restart[1])
    #     f = open(restart[1], "rb")
    #     size, M, data = self.GetBinData(0, f)
    #     steps = int(filesize / size)
    #     # 最終ステップのデータ
    #     size, M, data = self.GetBinData(size*(steps - 1), f)
    #     atoms = data.shape[0]
    #     if M.shape[1] == 3:  # triclinic
    #         xlo = M[0, 0] - np.min([0.0, M[0, 2], M[1, 2], M[0, 2] + M[1, 2]])
    #         xhi = M[0, 1] - np.max([0.0, M[0, 2], M[1, 2], M[0, 2] + M[1, 2]])
    #         ylo = M[1, 0] - np.min([0.0, M[2, 2]])
    #         yhi = M[1, 1] - np.max([0.0, M[2, 2]])
    #         zlo = M[2, 0]
    #         zhi = M[2, 1]
    #         lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
    #         a = lx
    #         b = np.sqrt(ly**2 + M[0, 2]**2)
    #         c = np.sqrt(lz**2 + M[1, 2]**2 + M[2, 2]**2)
    #         alpha = np.arccos((M[2, 2] * ly + M[0, 2] * M[1, 2]) / (b * c))
    #         beta = np.arccos(M[1, 2] / c)
    #         gamma = np.arccos(M[0, 2] / b)
    #         self.box = [[0.0, lx], [0.0, ly], [0.0, lz]]
    #         self.box_tri = M[:, 2]
    #     else:  # orthogonal
    #         alpha, beta, gamma = 0.5 * np.pi, 0.5 * np.pi, 0.5 * np.pi
    #         a = M[0, 1] - M[0, 0]
    #         b = M[1, 1] - M[1, 0]
    #         c = M[2, 1] - M[2, 0]
    #         self.box = M
    #     self.M = CalcMatrix(a, b, c, alpha, beta, gamma)
    #     # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
    #     self.atoms = np.zeros((atoms, 11), dtype=object)
    #     self.atoms[:, 0] = data[:, 0].astype(int)  # id
    #     self.atoms[:, 1] = data[:, 1].astype(int)  # type
    #     self.atoms[:, 3:6] = data[:, 3:6].astype(float)  # x y z
    #     self.atoms[:, 6:9] = data[:, 6:9].astype(float)  # vx vy vz
    #     M_ = np.linalg.inv(self.M)
    #     self.atoms[:, 3:6] = np.dot(self.atoms[:, 3:6].astype(float), M_)
    #     self.a, self.b, self.c = a, b, c
    #     self.alpha, self.beta, self.gamma = alpha, beta, gamma
    #     # 元素名とlabelをセットする
    #     # masses 0:type 1:mass 2:element 3:type
    #     # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
    #     for m in masses:
    #         self.atoms[self.atoms[:, 1] == m[0], 9] = m[2]
    #         self.atoms[self.atoms[:, 1] == m[0], 10] = m[2]

    def UpdateTypeCharge(self):
        """
        labelに基づいて、atomsのtype, charge, elemをupdateする
        """
        labels = np.unique(self.atoms[:, 10])
        self.masses = np.empty((0, 4), dtype=object)
        for k in FF.keys():
            if k in labels:
                self.masses = np.vstack(
                    (self.masses, [[0, FF[k][0], FF[k][2], k]]))
        self.masses[:, 0] = np.linspace(1, labels.shape[0], labels.shape[0])
        self.masses[:, 0] = self.masses[:, 0].astype(int)
        self.masses[:, 1] = self.masses[:, 1].astype(float)
        for m in self.masses:
            self.atoms[self.atoms[:, 10] == m[3], 1] = m[0]        # type
            self.atoms[self.atoms[:, 10] == m[3], 2] = FF[m[2]][1]  # charge
            self.atoms[self.atoms[:, 10] == m[3], 9] = FF[m[2]][2]  # elem
        self.atoms[:, 0:2] = self.atoms[:, 0:2].astype(int)
        self.atoms[:, 2:9] = self.atoms[:, 2:9].astype(float)

    def WriteLammpsData(self):
        """
        lammps data fileの生成
        """
        self.UpdateTypeCharge()
        datafile = args.basename + ".data"
        o = open(datafile, "w")
        o.write("TEST\n\n")
        # sec1
        o.write("{} atoms\n".format(self.atoms.shape[0]))
        # sec2
        o.write("{} atom types\n".format(self.masses.shape[0]))
        o.write("\n")
        # sec3
        o.write("{:<24.16e} ".format(self.box[0][0]))
        o.write("{:<24.16e} ".format(self.box[0][1]))
        o.write("xlo xhi\n")
        o.write("{:<24.16e} ".format(self.box[1][0]))
        o.write("{:<24.16e} ".format(self.box[1][1]))
        o.write("ylo yhi\n")
        o.write("{:<24.16e} ".format(self.box[2][0]))
        o.write("{:<24.16e} ".format(self.box[2][1]))
        o.write("zlo zhi\n")
        # if self.box_tri is not None:
        #     o.write("{:<24.16e} ".format(self.box_tri[0]))
        #     o.write("{:<24.16e} ".format(self.box_tri[1]))
        #     o.write("{:<24.16e} ".format(self.box_tri[2]))
        #     o.write(" xy xz yz\n".format(self.box_tri))
        o.write("\n")
        # sec4
        o.write("Masses\n\n")
        for m in self.masses:
            o.write("%4d" % (m[0]))   # id
            o.write("%12f" % (m[1]))  # mass
            o.write(" # %s " % (m[2]))  # elem
            o.write(" %s" % (m[3]))  # label
            o.write("\n")
        o.write("\n")
        # sec5
        o.write("Atoms\n\n")
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:elem 7:vx 8:vy 9:vz
        # xyz = np.dot(self.atoms[:, 3:6], self.M)
        xyz = self.atoms[:, 3:6]
        for i, a in enumerate(self.atoms):
            o.write("{:6d} ".format(a[0]))  # id
            o.write("{:4d} ".format(a[1]))  # type_id
            o.write("{:8.4f} ".format(a[2]))  # charge
            o.write("{:24.16e} ".format(xyz[i, 0]))  # x
            o.write("{:24.16e} ".format(xyz[i, 1]))  # y
            o.write("{:24.16e} ".format(xyz[i, 2]))  # z
            o.write("# {}".format(a[10]))   # label
            o.write("\n")
        o.write("\n")
        # sec5
        o.write("Velocities\n\n")
        for i, a in enumerate(self.atoms):
            o.write("{:6d} ".format(a[0]))  # id
            o.write("{:24.16e} ".format(a[6]))  # vx
            o.write("{:24.16e} ".format(a[7]))  # vy
            o.write("{:24.16e} ".format(a[8]))  # vz
            o.write("# {}".format(a[10]))   # label
            if i == self.knock_idx:
                o.write(", PKA")
            o.write("\n")
        o.write("\n")
        o.close()
        print(datafile, "was created.")


def Output(l):
    l.gamma = np.pi/2
    print("-" * 60)
    print("input cif:", args.ciffile)
    if args.ciffile is not None:
        restart = args.ciffile
    print("input restart:", restart)
    print("super cell:", args.replicate)
    print("output basename:", args.basename)
    print("-" * 60)
    print("-" * 60)
    print("total atoms:", l.atoms.shape[0])
    print("atomic type list:")
    for k in FF.keys():
        n = np.sum(l.atoms[:, 10] == k)
        print("{:>8s}: {}".format(k, n))
    print("-" * 60)
    s = "(a, b, c) = ({:.2f}, {:.2f}, {:.2f})"
    print(s.format(l.a, l.b, l.c))
    ang = np.array([l.alpha, l.beta, l.gamma]) * 180 / np.pi
    s = "(alpha, beta, gamma) = ({:.2f}, {:.2f}, {:.2f})"
    print(s.format(ang[0], ang[1], ang[2]))
    # V = np.dot(l.M[:, 0], np.cross(l.M[:, 1], l.M[:, 2]))
    V = l.a * l.b * l.c
    W = 0
    for a in l.atoms:
        W = W + FF[a[10]][0]
    W = W  # g/mol -> g
    print("volume: {:.2f} Ang3".format(V))
    print("weight: {:.2f} g/mol".format(W))
    print("density: {:.3f} g/cm3".format(W / 6.02214129e23 / (V * 1e-24)))


l = LammpsData()
if args.ciffile is not None:
    l.Cif2Data(args.ciffile)
l.WriteLammpsData()
Output(l)
