#!/usr/bin/env python

# # https://doi.org/10.1016/j.jnoncrysol.2004.08.247
# p = {}
# p["Si-Si"] = ["1 1", 0.007695, 2.0446, 3.7598]
# p["Si-O"] = ["1 2",1.99597, 2.6518, 1.628]
# p["O-O"] = ["2 2",0.023272, 1.3731, 3.791]

# coeff = 1.602176462e-19 / 4.184 * 1e-3 * 6.02214129e23
# # coeff = 23.06000                # ev -> kcal/mol
# fmt = "pair_coeff     {} morse {:12.4f} {:7.4f} {:7.2f} # {}"

# for k, v in p.items():
#     d0 = v[1] * coeff
#     print(fmt.format(v[0], d0, v[2], v[3], k))

# morse ab initio parametized
# https://doi.org/10.1063/1.153312
### atomic unit ###
# D[Eh], gamma[a0**-2], r0[a0]
p = {}
p["Si-Si"] = ["1 1", -0.0020846, 10.45517, 5.75038]
p["Si-O"] = ["1 2", 0.0019033, 11.15230, 4.63710]
p["O-O"] = ["2 2", 0.00024748, 12.07092, 7.17005]

a0 = 0.5291                     # angstrome, distance
Eh = 4.359744e-18               # J, energy
cal = 4.184                     # J, energy
mol = 6.02214129e23             # /mol, number
# coeff = 1.602176462e-19 / 4.184 * 1e-3 * 6.02214129e23
Ecoeff = Eh / cal * 1e-3 * 6.02214129e23  # kcal/mol
# coeff = 23.06000                # ev -> kcal/mol
fmt = "pair_coeff     {} morse {:12.4f} {:7.4f} {:7.2f} # {}"

for k, v in p.items():
    d0 = v[1] * Ecoeff
    d1 = v[2] / v[3] * 0.5 / a0
    d2 = v[3] * a0
    print(fmt.format(v[0], d0, d1, d2, k))
