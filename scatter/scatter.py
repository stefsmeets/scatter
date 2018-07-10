#!/usr/bin/env python

#    Scatter - A python tool to plot and output atomic scattering factors
#    Copyright (C) 2018 Stef Smeets
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from __future__ import print_function
try:
    range = xrange
except NameError:
    pass

from random import random as r
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import argparse

from it_table_4322 import it_table_4322
from it_table_4323 import it_table_4323
from peng1998 import peng1998
from wk1995 import wk1995
from dt1969 import dt1969
from atomic_radii import atomic_radii


XRS = {  # a1 b1 a2 b2 a3 b3 a4 b4 c
    'Al': "FORMGN Al    6.420203.038701.900200.742601.5936031.54721.9646085.08861.11510",
    'C': "FORMGN C     2.3100020.84391.0200010.20751.588600.568700.8650051.65120.21560",
    'Ca': "FORMGN Ca    8.6266010.44217.387300.659901.5899085.74841.02110178.4371.37510",
    'F': "FORMGN F     3.5392010.28252.641204.294401.517000.261501.0243026.14760.27760",
    'Ge': "FORMGN Ge    16.08162.850906.374700.251603.7068011.44683.6830054.76252.13130",
    'H': "FORMGN H     0.4899220.65930.262007.740390.1967749.55190.049882.201590.00131",
    'N': "FORMGN N     12.21260.005703.132209.893302.0125028.99751.166300.58260-11.529",
    'Na': "FORMGN Na    4.762603.285003.173608.842201.267400.313601.11280129.4240.67600",
    'O': "FORMGN O     3.0448513.22772.286805.701101.546300.323900.8670032.90890.25080",
    'P': "FORMGN P     6.434501.906704.1791027.15701.780000.526001.4908068.16451.11490",
    'Si': "FORMGN Si    6.291502.438603.0353032.33371.989100.678501.5410081.69371.14070",
    'Zn': "FORMGN Zn    14.07433.265507.031800.233305.1652010.31632.4100058.70971.30410",
    'Zr': "FORMGN Zr    17.87651.2761810.948011.91605.417320.117623.6572187.66272.06929"
}

elements = ['H',    'He',   'Li',   'Be',   'B',    'C',    'N',    'O',    'F',
            'Ne',   'Na',   'Mg',   'Al',   'Si',   'P',    'S',    'Cl',   'Ar',
            'K',    'Ca',   'Sc',   'Ti',   'V',    'Cr',   'Mn',   'Fe',
            'Co',   'Ni',   'Cu',   'Zn',   'Ga',   'Ge',   'As',   'Se',   'Br',
            'Kr',   'Rb',   'Sr',   'Y',    'Zr',   'Nb',   'Mo',   'Tc',   'Ru',
            'Rh',   'Pd',   'Ag',   'Cd',   'In',   'Sn',   'Sb',   'Te',   'I',
            'Xe',   'Cs',   'Ba',   'La',   'Ce',   'Pr',   'Nd',   'Pm',   'Sm',
            'Eu',   'Gd',   'Tb',   'Dy',   'Ho',   'Er',   'Tm',   'Yb',   'Lu',
            'Hf',   'Ta',   'W',    'Re',   'Os',   'Ir',   'Pt',   'Au',   'Hg',
            'Tl',   'Pb',   'Bi',   'Po',   'At',   'Rn',   'Fr',   'Ra',   'Ac',
            'Th',   'Pa',   'U',    'Np',   'Pu',   'Am',   'Cm',   'Bk',   'Cf']

ions = ['H1-', 'Li1+', 'Be2+', 'Cval', 'O1-', 'O2-', 'F1-', 'Na1+', 'Mg2+',
        'Al3+', 'Sival', 'Si4+', 'Cl1-', 'K1+', 'Ca2+', 'Sc3+', 'Ti2+',
        'Ti3+', 'Ti4+', 'V2+', 'V3+', 'V5+', 'Cr2+', 'Cr3+', 'Mn2+',
        'Mn3+', 'Mn4+', 'Fe2+', 'Fe3+', 'Co2+', 'Co3+', 'Ni2+', 'Ni3+',
        'Cu1+', 'Cu2+', 'Zn2+', 'Ga3+', 'Ge4+', 'Br1-', 'Rb1+', 'Sr2+',
        'Y3+', 'Zr4+', 'Nb3+', 'Nb5+', 'Mo3+', 'Mo5+', 'Mo6+', 'Ru3+',
        'Ru4+', 'Rh3+', 'Rh4+', 'Pd2+', 'Pd4+', 'Ag1+', 'Ag2+', 'Cd2+',
        'In3+', 'Sn2+', 'Sn4+', 'Sb3+', 'Sb5+', 'I1-', 'Cs1+', 'Ba2+',
        'La3+', 'Ce3+', 'Ce4+', 'Pr3+', 'Pr4+', 'Nd3+', 'Pm3+', 'Sm3+',
        'Eu2+', 'Eu3+', 'Gd3+', 'Tb3+', 'Dy3+', 'Ho3+', 'Er3+', 'Tm3+',
        'Yb2+', 'Yb3+', 'Lu3+', 'Hf4+', 'Ta5+', 'W6+', 'Os4+', 'Ir3+',
        'Ir4+', 'Pt2+', 'Pt4+', 'Au1+', 'Au3+', 'Hg1+', 'Hg2+', 'Tl1+',
        'Tl3+', 'Pb2+', 'Pb4+', 'Bi3+', 'Bi5+', 'Ra2+', 'Ac3+', 'Th4+',
        'U3+', 'U4+', 'U6+', 'Np3+', 'Np4+', 'Np6+', 'Pu3+', 'Pu4+',
        'Pu6+']

other = ['D', 'Es', 'Fm', 'Md', 'No', 'Lr', 'NULL']

# n = r()*len(it_table_4322.keys()) // 1
# key = it_table_4322.keys()[int(n)]

# print it_table_4322[key]


def calc_s(ld, r):
    return 4 * np.pi * (1/ld) * r


def gaussian_fit(a, b, s):
    """General Gaussian"""
    return a * np.exp(-b * s**2)


def plot_sf_atoms(atoms, s, kind="xray"):
    for atom in atoms:
        data = tables[atom]
        y = calc_sf(atom, data, s, kind)
        plt.plot(s, y, label=atom)

    plt.ylabel("f")
    plt.xlabel(r"$\sin(\theta)/\lambda (1/\mathrm{\AA})$")
    plt.legend()
    plt.show()


def calc_sf(atom, data, s, kind="xray"):
    if kind == 'xray':
        return calc_sf_xray(atom, data, s)
    elif kind == "electron":
        return calc_sf_electron(atom, data, s)
    else:
        raise NameError


def calc_sf_electron(atom,  data,   s):
    """scattering factor function for electron/it432x table"""
    total = None
    for i in range(5):
        a, b, dZ = data[2+i], data[7+i], data[1]
        y = gaussian_fit(a, b, s)

        if total is None:
            total = y
        else:
            total += y
    return total


def calc_sf_xray(atom,  data,   s):
    """Scattering factor function for xray/wk1995 table"""
    total = None
    for i in range(5):
        a, b, c = data[0+i], data[5+i], data[10]   # wk95
        y = gaussian_fit(a, b, s)

        if total is None:
            total = y + c
        else:
            total += y
    return total


def print_xy_atoms(atoms, s, kind="xray"):
    ys = []
    for atom in atoms:
        data = tables[atom]
        ys.append(calc_sf(atom, data, s, kind))
    print_xy(atoms, s, ys)


def print_xy(atoms, s,  ys):

    print("\n      ", end=' ')
    for atom in atoms:
        print("{:>10s}".format(atom), end=' ')
    for i, val in enumerate(s):
        print("\n{:6.2f}".format(val), end=' ')
        for j in range(len(atoms)):
            print("{:10.5f}".format(ys[j][i]), end=' ')


def check_consistency(atoms, s, plot=False,   show=False, threshold=0):
    for atom in atoms:
        total1 = None
        data1 = it_table_4322[atom]
        for i in range(5):
            a, b, dZ = data1[2+i], data1[7+i], data1[1]
            y = gaussian_fit(a, b, s)
            if total1 == None:
                total1 = y
            else:
                total1 += y

        data2 = it_table_4323[atom]
        total2 = None
        for i in range(5):
            a, b, dZ = data2[2+i], data2[7+i], data2[1]
            y = gaussian_fit(a, b, s)
            if total2 == None:
                total2 = y
            else:
                total2 += y

        r = sum((total1-total2)**2)
        print("%4s %7.3f" % (atom, r))

        if r > threshold:
            if show:
                print_table([atom], data1)
                print_table([atom], data2)
            if plot:
                plt.plot(s, total1)
                plt.plot(s, total2)
                plt.show()


def print_combine_tables(atoms):
    for atom in atoms:
        if atom in wk1995:
            data = wk1995[atom]
            print("{ \"%s\",\n" % atom, end=' ')
            print("    { %f, %f, %f, %f, %f },   \n" % (data[0],    data[1],    data[2],    data[3],    data[4]), end=' ')
            print("    { %f, %f, %f, %f, %f },   \n" % (data[5],    data[6],    data[7],    data[8],    data[9]), end=' ')
            print("      %f,                     \n" % data[10], end=' ')
        else:
            print('atom not found ?', atom)
            exit()
        if atom in it_table_4322:
            data = it_table_4322[atom]
            print("    { %f, %f, %f, %f, %f },   \n" % (data[2],    data[3],    data[4],    data[5],    data[6]), end=' ')
            print("    { %f, %f, %f, %f, %f },   \n" % (data[7],    data[8],    data[9],    data[10],   data[11]), end=' ')
        else:
            print("    { 0.0, 0.0, 0.0, 0.0, 0.0 }, \n", end=' ')
            print("    { 0.0, 0.0, 0.0, 0.0, 0.0 }, \n", end=' ')
        if atom in it_table_4323:
            data = it_table_4323[atom]
            print("    { %f, %f, %f, %f, %f },   \n" % (data[2],    data[3],    data[4],    data[5],    data[6]), end=' ')
            print("    { %f, %f, %f, %f, %f } }, \n" % (data[7],    data[8],    data[9],    data[10],   data[11]), end=' ')
        else:
            print("    { 0.0, 0.0, 0.0, 0.0, 0.0 }, \n", end=' ')
            print("    { 0.0, 0.0, 0.0, 0.0, 0.0 } }, \n", end=' ')


def xrs2table(d):
    e = {}
    a5 = 0.0
    b5 = 0.0
    for key, line in d.items():
        label = line[0:11].strip()
        a1, b1, a2, b2, a3, b3, a4, b4, c = [
            float(line[13+7*i:13+7+7*i]) for i in range(9)]
        e['xrs'+key] = [a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c]
    return e


def combine_scfacts(atom1, atom2, ratio):
    ratio = float(ratio)

    data1 = tables[atom1]
    data2 = tables[atom2]

    sfact = []

    for n in range(12):
        item = data1[n]*ratio + data2[n]*(1-ratio)
        sfact.append(item)

    atoms = [atom1, atom2, 'custom']
    return atoms, sfact


def print_table(atoms, kind, custom_data=None):
    if kind == "xray":
        print("""name 
    a1 a1 a3 a4 a5
    b1 b2 b3 b4 b5
    c""")
    if kind == "electron":
        print("""name 
    a1 a1 a3 a4 a5
    b1 b2 b3 b4 b5""")
    for atom in atoms:
        if custom_data:
            data = custom_data
        else:
            data = tables[atom]
        if kind == "xray":
            print("{ \"%s\",\n" % atom, end=' ')
            print("    { %f, %f, %f, %f, %f },   \n" % (data[0],    data[1],    data[2],    data[3],    data[4]), end=' ')
            print("    { %f, %f, %f, %f, %f },   \n" % (data[5],    data[6],    data[7],    data[8],    data[9]), end=' ')
            print("      %f,                     \n" % data[10], end=' ')
        if kind == "electron":
            print("{ \"%s\",\n" % atom, end=' ')
            print("    { %f,    %f, %f, %f, %f },   \n" % (data[2], data[3],    data[4],    data[5],    data[6]), end=' ')
            print("    { %f,    %f, %f, %f, %f } },\n" % (data[7],  data[8],    data[9],    data[10],   data[11]), end=' ')


def add_us():
    atoms, tables['custom'] = combine_scfacts(*atoms)
    check_consistency(atoms, s, True, True, 0.01)


def gaussian(x, height, width, offset):
    return height*np.exp(-(x)**2*width) + offset

def four_gaussians(x, h1, w1, h2, w2, h3, w3, h4, w4, offset=0):
    return (gaussian(x, h1, w1, offset=0) +
            gaussian(x, h2, w2, offset=0) +
            gaussian(x, h3, w3, offset=0) + 
            gaussian(x, h4, w4, offset=0) + offset)

def five_gaussians(x, h1, w1, h2, w2, h3, w3, h4, w4, h5, w5, offset=0):
    return (gaussian(x, h1, w1, offset=0) +
        gaussian(x, h2,  w2, offset=0) +
        gaussian(x, h3,  w3, offset=0) + 
        gaussian(x, h4,  w4, offset=0) + 
        gaussian(x, h5,  w5, offset=0) + offset)

def print_table_shelx(atoms, s, kind):
    # from IPython import embed
    # embed()

    sfacs = []

    from scipy import optimize

    for atom in atoms:
        if atom not in tables:
            continue
        data = tables[atom]
        if kind == "electron":
            values = data[2], data[7], data[3], data[8], data[4], data[9], data[5], data[10], data[6], data[11]
            
            args = s, five_gaussians(s, *values, offset=0)
            
            def errfunc(p, x, y):
                # print(p.shape, x.shape, y.shape)
                ret = (four_gaussians(x, *p) - y)**2
                # print(ret.shape, np.sum(ret))
                return ret
            
            x_init = [0.3352, 0.4001, 0.8630, 3.0532, 2.5725, 17.1148, 2.0490, 59.0519]

            result = optimize.least_squares(errfunc, x_init[:], args=args)
            x_values = result.x

            if result.success:
                print("{name} minimized nfev: {nfev}, njev: {njev}, optimality: {optimality}, cost: {cost}".format(name=atom, **result))
            else:
                raise ValueError(result.message)

            element = ''.join([i for i in atom if i.isalpha()])
            covrad = atomic_radii[element][1]
            wght = atomic_radii[element][2]

            h1, w1, h2, w2, h3, w3, h4, w4 = x_values
            sfac = "SFAC {atom:4s} {h1:7.4f} {w1:7.4f} {h2:7.4f} {w2:7.4f} {h3:7.4f} {w3:7.4f} {h4:7.4f} {w4:7.4f} =\n          {no_data:7.4f} {no_data:7.4f} {no_data:7.4f} {covrad:7.4f} {wght:7.4f}".format(
                atom=atom, h1=h1, w1=w1, h2=h2, w2=w2, h3=h3, w3=w3, h4=h4, w4=w4, no_data=0.0, covrad=covrad, wght=wght)
            sfacs.append(sfac)

            plt.plot(s, four_gaussians(s, *x_values), label="{} fitted (4 params)".format(atom))
            plt.plot(s, five_gaussians(s, *values), label="{} tabulated (5 params)".format(atom))

        else:
            raise NotImplementedError

    plt.ylabel("f")
    plt.xlabel(r"$\sin(\theta)/\lambda (1/\mathrm{\AA})$")
    plt.legend()
    plt.show()

    print()
    for sfac in sfacs:
        print(sfac)


def print_table_topas(atoms):
    print("x-ray {")
    for atom in atoms:
        if atom not in tables:
            continue
        data = tables[atom]
        if len(data) == 12:
            print('%8s' % atom, end=' ')
            print('%f  %f  %f  %f  %f' % (data[2],  data[3],    data[4],    data[5],    data[6]), end=' ')
            print('%f' % 0.0, end=' ')
            print('%f  %f  %f  %f  %f' % (data[7],  data[8],    data[9],    data[10],   data[11]))
        else:
            print('%8s' % atom, end=' ')
            print('%f  %f  %f  %f  %f' % (data[0],  data[1],    data[2],    data[3],    data[4]), end=' ')
            print('%f' % 0.0, end=' ')
            print('%f  %f  %f  %f  %f' % (data[5],  data[6],    data[7],    data[8],    data[9]))
    print('}')


def main():
    description = """Notes:
    - Passing 'all' as an argument adds all atoms.
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("args",
                        type=str, metavar="atoms", nargs='*',
                        help="List of atoms")

    parser.add_argument("-t", "--table",
                        action="store", type=str, dest="table",
                        help="Scattering factor table to use (xray,electron,wk1995,it4322,it4323,peng1998,dt1969). Defaults: xray")

    parser.add_argument("-r", "--range", metavar="val",
                        action="store", type=float, nargs=3, dest="s_range",
                        help="sin(th)/lambda range. Requires 3 values: start finish step")

    parser.add_argument("-w", "--raw",
                        action="store_true", dest="print_raw_data",
                        help="Outputs raw data used for plotting graphs")

    parser.add_argument("-n", "--noplot",
                        action="store_false", dest="plot",
                        help="Skips plotting procedure")

    parser.add_argument("-m", "--merged",
                        action="store_true", dest="merged",
                        help="Combines all scattering factor tables. Used for updating the tables in focus/atoms.h")

    parser.add_argument("-c", "--xrs",
                        action="store_true", dest="xrs",
                        help="Plots scattering factors for XRS, including FORMGN lines for drcard.dat")

    parser.add_argument("--topas",
                        action="store_true", dest="topas",
                        help="Print table compatible with topas (save as atmscat.cpp in topas4 root dir)")
    
    parser.add_argument("--shelx",
                        action="store_true", dest="shelx",
                        help="Print SFAC cards for SHELX")

    parser.set_defaults(table="xray",
                        plot=True,
                        print_raw_data=False,
                        merged=False,
                        combine=False,
                        s_range=[0, 2, 0.01],
                        xrs=False,
                        print_table=True,
                        topas=False)

    options = parser.parse_args()
    args = options.args

    global tables

    s = np.arange(*options.s_range)

    atoms = []
    if 'elements' in args:
        args.remove('elements')
        atoms += elements
    if 'ions' in args:
        args.remove('ions')
        atoms += ions
    if 'other' in args:
        args.remove('other')
        atoms += other
    if 'all' in args:
        args.remove('all')
        atoms += elements + ions + other
    if 'xrsall' in args:
        args.remove('xrsall')
        atoms += XRS.keys()
    atoms += args

    if options.table in ('xray', 'wk1995'):
        kind = "xray"
        tables = wk1995
    elif options.table == 'electron':
        kind = "electron"
        tables = dict(peng1998.items() + it_table_4322.items())
    elif options.table == 'it_table_4322':
        kind = "electron"
        tables = it_table_4322
    elif options.table == 'it_table_4323':
        kind = "electron"
        tables = it_table_4323
    elif options.table == 'peng1998':
        kind = "electron"
        tables = peng1998
    else:
        raise NameError('Unknown scattering factor table: {}'.format(options.table))

    if options.xrs:
        options.print_table = False
        tables = dict(tables.items() + xrs2table(XRS).items())
        print('Add these lines to drcard.dat and run datrdn\n')
        for atom in atoms:
            print(XRS.get(atom, 'Atom {} not in the table!'.format(atom)))

        atoms += ['xrs'+atom for atom in atoms if atom in XRS.keys()]

    if options.merged:
        print_combine_tables(atoms)
    elif options.topas:
        print_table_topas(atoms)
    elif options.shelx:
        print_table_shelx(atoms, s, kind)
    else:
        if options.print_table:
            print_table(atoms, kind)

        if options.plot:
            plot_sf_atoms(atoms, s, kind)

        if options.print_raw_data:
            print_xy_atoms(atoms, s, kind)


if __name__ == '__main__':
    main()
