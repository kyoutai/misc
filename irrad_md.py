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
group.add_argument('--cascade', metavar=('*.data', '*.bin', '600'), nargs=3,
                   default=None, help="cascade run *.data *.bin bullet_energy_(eV)")
group.add_argument('--equilib', metavar=('*.data', '*.bin', '20'), nargs=3,
                   default=None, help="equilibrium run *.data *.bin NVE_radius")
par.add_argument('-s', '--replicate', nargs=3, default=[1, 1, 1], metavar=1,
                 type=int,  help="supercell cell; x y z")
par.add_argument('-o', '--basename', required=True,
                 help="output basename lammps input; *.data and *.in")
args = par.parse_args()


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

LAMMPSIN = """
boundary        p p p
units	        real

atom_style      charge
box             tilt large
read_data       {base}.data
replicate       {replicate}
pair_style      reax/c NULL checkqeq no
pair_coeff      * * ffieldCSH.reax {elems}
fix             1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c
neighbor        1 bin
neigh_modify    every 1 check yes
"""

LAMMPSIN_INIT = """
velocity        all create 300 777
fix	        10 all npt temp 300 300 100 tri 0.0 0.0 1000
dump            1 all custom 1000 {base}.bin id type element x y z vx vy vz
dump_modify     1 element {elems} sort id
compute         1 all displace/atom
compute         2 all reduce max c_1[4]
thermo          10
thermo_style    custom step dt time temp press pe ke etotal enthalpy density c_2

timestep        0.5
run             20000
"""

LAMMPSIN_IRR = """
group pka       id {pka_id}      # {pka_elem}
group non-pka   subtract all pka # pka以外の原子

fix             10 all dt/reset 1 0.0001 0.5 0.05 units box   # variable dt
fix             11 all npt temp 300 300 1000000000000 tri 0 0 1000

dump             1 all custom 1000 {base}.bin id type element x y z vx vy vz
dump_modify      1 element {elems} sort id

compute         1 non-pka displace/atom
compute         2 non-pka reduce max c_1[4]

# displacementのmaxが2.0 Ang.以下になるまでdt/resetで計算
thermo          100
print           "CASCADE LOOP START, {base}"
label           loop
  compute         3 all displace/atom
  compute         4 all reduce max c_3[4]
  thermo_style    custom step dt time temp press pe ke etotal enthalpy density c_2 c_4
  run             1000
  variable        dispmax equal c_4
  if              "${{dispmax}} < 2.0" then "jump SELF break_loop"
  uncompute        3
  uncompute        4
  jump            SELF loop
label           break_loop
print           "CASCADE LOOP END, {base}"
"""

LAMMPSIN_EQUILIBRIUM = """
thermo          100
neigh_modify    every 10 check yes
group sink      type 1 2
group damage    type 3 4

fix             10 all npt temp 300 300 100 tri 0 0 1000
fix             12 sink temp/berendsen 300 300 100
fix             13 sink recenter INIT INIT INIT

dump            1 all custom 1000 {base}.bin id type element x y z vx vy vz
dump_modify     1 element {elems} sort id

compute         1 all displace/atom
compute         2 all reduce max c_1[4]
thermo_style    custom step dt time temp press pe ke etotal enthalpy density c_2

timestep        0.5
run             4000
print           "EQUILIBRIUM RUN START, {base}"
run             5000
print           "EQUILIBRIUM RUN END, {base}"
"""


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
        obj = re.search("_atom_site_type_symbol(.+)", body, re.DOTALL)
        data = obj.group(1)
        data = re.sub("\(\d+?\)", "", data)
        data = np.array(data.split(), dtype=object).reshape(-1, 8)
        data = data[:, [7, 2, 3, 4]]  # element, x, y, z
        data[:, 1:4] = data[:, 1:4].astype(float)
        data[:, 1:4] = data[:, 1:4] - np.floor(data[:, 1:4].astype(float))
        self.data = data


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
        if np.round(alpha * 180 / np.pi) == 90 and\
           np.round(beta * 180 / np.pi) == 90 and\
           np.round(gamma * 180 / np.pi) == 90:
            self.box[0, :] = np.array([[0.0, a], [0.0, b], [0.0, c]])
        else:  # triclnic
            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b * c * np.cos(alpha) - xy * xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)
            self.box = np.array([[0.0, lx], [0.0, ly], [0.0, lz]])
            self.box_tri = np.array([xy, xz, yz])
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
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        # 元素名とlabelをセットする
        # masses 0:type 1:mass 2:element 3:type
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        for m in masses:
            self.atoms[self.atoms[:, 1] == m[0], 9] = m[2]
            self.atoms[self.atoms[:, 1] == m[0], 10] = m[2]

    def GetBinData(self, offset, f):
        f.seek(offset)
        head_size = 16 + 4 + 4*6 + 8*6 + 4 + 4
        step, atoms = struct.unpack("ll", f.read(16)) # step, atoms
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
            

    def Bin2Data(self, restart):
        """
        binaryで保存したデータの読み込み
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
        # bin dataの読み込み
        filesize = os.path.getsize(restart[1])
        f = open(restart[1], "rb")
        size, M, data = self.GetBinData(0, f)
        steps = int(filesize / size)
        # 最終ステップのデータ
        size, M, data = self.GetBinData(size*(steps - 1), f)
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
        else:  # orthogonal
            alpha, beta, gamma = 0.5 * np.pi, 0.5 * np.pi, 0.5 * np.pi
            a = M[0, 1] - M[0, 0]
            b = M[1, 1] - M[1, 0]
            c = M[2, 1] - M[2, 0]
            self.box = M
        self.M = CalcMatrix(a, b, c, alpha, beta, gamma)
        # 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
        self.atoms = np.zeros((atoms, 11), dtype=object)
        self.atoms[:, 0] = data[:, 0].astype(int) # id
        self.atoms[:, 1] = data[:, 1].astype(int) # type
        self.atoms[:, 3:6] = data[:, 3:6].astype(float)  # x y z
        self.atoms[:, 6:9] = data[:, 6:9].astype(float)  # vx vy vz
        M_ = np.linalg.inv(self.M)
        self.atoms[:, 3:6] = np.dot(self.atoms[:, 3:6].astype(float), M_)
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
        o.write("xlo xhi\n".format(self.box[0][1]))
        o.write("{:<24.16e} ".format(self.box[1][0]))
        o.write("{:<24.16e} ".format(self.box[1][1]))
        o.write("ylo yhi\n".format(self.box[0][1]))
        o.write("{:<24.16e} ".format(self.box[2][0]))
        o.write("{:<24.16e} ".format(self.box[2][1]))
        o.write("zlo zhi\n".format(self.box[0][1]))
        if self.box_tri is not None:
            o.write("{:<24.16e} ".format(self.box_tri[0]))
            o.write("{:<24.16e} ".format(self.box_tri[1]))
            o.write("{:<24.16e} ".format(self.box_tri[2]))
            o.write(" xy xz yz\n".format(self.box_tri))
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
        xyz = np.dot(self.atoms[:, 3:6], self.M)
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

    def WriteLammpsIn(self):
        """
        lammps input fileの生成
        """
        infile = args.basename + ".in"
        o = open(infile, "w")
        elems = " ".join(self.masses[:, 2])
        supercell = " ".join([str(d) for d in args.replicate])
        lammps_in = LAMMPSIN.format(base=args.basename,
                                    elems=elems,
                                    replicate=supercell)
        if args.ciffile is not None:
            lammps_init = LAMMPSIN_INIT
            lammps_init = lammps_init.format(base=args.basename,
                                             elems=elems)
        if args.cascade is not None:
            lammps_init = LAMMPSIN_IRR
            lammps_init = lammps_init.format(base=args.basename,
                                             elems=elems,
                                             pka_id=self.knock_idx + 1,
                                             pka_elem=l.atoms[self.knock_idx, 9])
        if args.equilib is not None:
            lammps_init = LAMMPSIN_EQUILIBRIUM
            lammps_init = lammps_init.format(base=args.basename,
                                             elems=elems)
        o.write(lammps_in)
        o.write(lammps_init)
        o.close()
        print(infile, "was created.")


def AddIrrDamageCascade(d):
    # 中性子cross sectionを考慮してランダムに原子を選択
    # atoms 0:id 1:type 2:charge 3:x 4:y 5:z 6:vx 7:vy 8:vz 9:elem 10:label
    XS_sum, FF_XS_p = 0, {}
    for x in np.unique(d.atoms[:, 9]):
        XS_sum += FF_XS[x]
        FF_XS_p[x] = FF_XS[x]
    # 全原子種の散乱断面積を規格化
    for k in FF_XS_p.keys():
        FF_XS_p[k] = FF_XS_p[k]/XS_sum
    while True:    
        p = np.floor(np.random.random() * d.atoms.shape[0]).astype(int)
        # 規格化した散乱断面積以下なら選択
        if np.random.random() < FF_XS_p[d.atoms[p, 9]]: 
            break
    # ランダムな角度を生成
    rnd = np.random.random(2)
    theta, phi = rnd[0] * np.pi, 2 * rnd[1] * np.pi
    d.knock_idx = p
    d.theta, d.phi = theta * 180 / np.pi, phi * 180 / np.pi
    E = float(args.cascade[2]) * 1.602176462e-19  # eV -> J 
    m = FF[d.atoms[p, 10]][0]          # g/mol
    m = m / 6.02214129e23 * 1e-3       # g/mol -> kg
    v_norm = np.sqrt(2 * E / m)        # m/s
    v_norm = v_norm * 1e-5             # m/s -> Ang./fs
    vx = v_norm * np.sin(theta) * np.cos(phi)
    vy = v_norm * np.sin(theta) * np.sin(phi)
    vz = v_norm * np.cos(theta)
    d.atoms[p, 6:9] = [vx, vy, vz] # update velocity
    d.irr_v = np.array([vx, vy, vz])  # Ang/fs


def AddIrrDamageEquilibrium(d):
    data = open(args.equilib[0]).read()
    obj = re.search(".*PKA", data)
    if obj is  None:
        print("PKA atom is nothing.")
        print("PKA atom must be defined in {}.".format(args.equilib[0]))
        print("Please assign *.data file generated for cascade run.")
    idx = int(obj.group(0).split()[0])
    p = np.where(d.atoms[:, 0] == idx)[0]
    d.knock_idx = p
    # 運動エネルギーを変えた原子まわりのtypeを変更
    radius = float(args.equilib[2]) # NVEの半径
    xyz = d.atoms[:, 3:6].astype(float)
    r = xyz[p, 0:3] - xyz[:, 0:3]
    r = np.dot(r - np.rint(r), d.M)
    r = np.sqrt(np.sum(r**2, axis=1))
    idx = np.where(r < radius)[0]
    # damage zoneのlabelを*_irrにする
    for i in idx:
        d.atoms[i, 10] = d.atoms[i, 10] + "_irr"
    

def AddIrrDamageDebug(d):
    """
    特定の原子に特定の速度を与える。debug用
    """
    # 最近接原子方向に速度ベクトルを与える time stepを決めるdebug用
    E = 600 # eV
    NEV_radius = 20 # Ang
    args.irr = [E, NEV_radius]
    p = 4091 # O lammps id=4092 Si lammps id=4089
    d.knock_idx = p
    # 再近接原子を探索して再近接原子のベクトルからthetaとphiを決める
    r = d.atoms[p, 3:6] - d.atoms[:, 3:6]
    r = r - np.rint(r.astype(float))
    r = np.dot(r, d.M)
    r = np.sum(r ** 2, axis=1)
    idx = np.argsort(r)[1]
    print("The nearest atom id: {} (VMD={})".format(idx+1, idx))
    vec = d.atoms[idx, 3:6] - d.atoms[p, 3:6]
    vec = np.dot(vec, d.M)
    vec_norm = np.sqrt(np.sum(vec ** 2))
    theta = np.arccos(vec[2]/vec_norm)
    phi = np.sign(vec[1]) * np.arccos(vec[0]/np.sqrt(vec[1]**2 + vec[2]**2))
    # 運動エネルギーからkock-onする速度を計算
    E = E * 1.602176462e-19         # eV -> J 
    m = FF[d.atoms[p, 10]][0]       # g/mol
    m = m / 6.02214129e23 * 1e-3    # g/mol -> kg
    v_norm = np.sqrt(2 * E / m)     # m/s
    v_norm = v_norm * 1e-5          # m/s -> Ang./fs
    vx = v_norm * np.sin(theta) * np.cos(phi)
    vy = v_norm * np.sin(theta) * np.sin(phi)
    vz = v_norm * np.cos(theta)
    d.atoms[p, 6:9] = [vx, vy, vz] # update velocity
    # 運動エネルギーを変えた原子まわりのtypeを変更
    xyz = d.atoms[:, 3:6].astype(float)
    r = xyz[p, 0:3] - xyz[:, 0:3]
    r = np.dot(r - np.rint(r), d.M)
    r = np.sqrt(np.sum(r**2, axis=1))
    idx = np.where(r < NEV_radius)[0]
    # damage zoneのlabelを*_irrにする
    for i in idx:
        d.atoms[i, 10] = d.atoms[i, 10] + "_irr"
    d.irr_v = np.array([vx, vy, vz])  # Ang/fs
    d.theta, d.phi = theta, phi


def Output(l):
    print("-" * 60)
    print("input cif:", args.ciffile)
    if args.ciffile is not None:
        restart = args.ciffile
    if args.cascade is not None:
        restart = args.cascade[0:2]
    if args.equilib is not None:
        restart = args.equilib[0:2]
    print("input restart:", restart)
    print("super cell:", args.replicate)
    print("output basename:", args.basename)
    print("-" * 60)
    if args.cascade is not None:
        print("irradiation energy: {} ev".format(args.cascade[2]))
        s = "(vx, vy, vz) = ({:.6f}, {:.6f}, {:.6f}) Ang./fs"
        print(s.format(l.irr_v[0], l.irr_v[1], l.irr_v[2]))
        s = "|v| = {:.6f} Ang./fs"
        print(s.format(np.sqrt(np.sum(l.irr_v**2))))
        s = "(theta, phi) =  ({:.2f} {:.2f})"
        print(s.format(l.theta, l.phi))
        print("knock-on atom id: {}".format(l.knock_idx + 1))
        print("knock-on atom index in VMD: {}".format(l.knock_idx))
        print("knock-on atom type: {}".format(l.atoms[l.knock_idx, 10]))
    else:
        print("No-irradiation")
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
    V = np.dot(l.M[:, 0], np.cross(l.M[:, 1], l.M[:, 2]))
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
if args.cascade is not None:
    l.Bin2Data(args.cascade[0:2])
    AddIrrDamageCascade(l)
if args.equilib is not None:
    l.Bin2Data(args.equilib[0:2])
    AddIrrDamageEquilibrium(l)

l.WriteLammpsData()
l.WriteLammpsIn()
Output(l)
