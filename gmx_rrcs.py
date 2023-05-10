# -*- coding: utf-8 -*-

'''
A script can calculate residues-residues contact scores (rrcs) of a trajectory file created by GROMACS.
Au : Chen Zhenghan
'''

import MDAnalysis as mda
import sys
import os
import subprocess
import re
import itertools
import numpy as np
import timeit
import matplotlib.pyplot as plt
import timeit

vanderwalls_radius = {'C': 1.9080, 'N': 1.8240, 'O': 1.6612, 'F': 1.7500,
                      'S': 2.0000, 'P': 2.1000, 'Cl': 1.9480, 'Br': 2.2200,
                      'I': 2.3500}

basic = {'plot': 'matplotlib',
         'r_min': 3.23,
         'r_max': 4.63,
         'rrcs': 'rrcs',
         'bt': 0.0,
         'et': 9999999.0,
         'dt': 0.1,
         'res': None,
         'top': None,
         'traj': None,
         'inter': False,
         'res_num': 0}

class InputFileError(FileNotFoundError):
    pass


class ParameterWrongError(Exception):
    pass


contact_score = {}
res_num_dict = {}
res_list = []
res_pair = []

def input_para(input_file):
    "Read in the parameter file."
    if not os.path.exists(input_file):
        raise InputFileError
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            para, si = line.strip().split(';', 1)[0].split('=')
            para = para.strip()
            si = si.strip()
            if not para:
                continue
            if para not in basic.keys():
                raise ParameterWrongError(para)
            if para in ('r_min', 'r_max', 'et', 'bt', 'dt'):
                basic[para] = float(si)
            else:
                basic[para] = si


def check():
    "Verify the correct format of the parameter file and initialize the Universe."
    # plot
    if basic['plot'] not in ('matplotlib', 'gnuplot', ''):
        raise ParameterWrongError('plot')
    # r_min,r_max
    basic['r_min'] = max(basic['r_min'], 0)
    basic['r_max'] = max(basic['r_min'], basic['r_max'])
    # rrcs
    if basic['rrcs'] not in ('rrcs', 'rc'):
        raise ParameterWrongError('rrcs')
    # file_exists
    if not os.path.exists(basic['top']):
        raise InputFileError(basic['top'])
    if not os.path.exists(basic['traj']):
        raise InputFileError(basic['traj'])
    else:
        if not os.path.exists(basic['res']):
            raise InputFileError(basic['res'])
        set_res(basic['res'])
    global u, dt_min
    u = mda.Universe(basic['top'], basic['traj'])
    basic['res_num'] = len(u.residues)
    # ti
    dt_min = (u.trajectory.time[-1] - u.trajectory.time[0]) / (len(u.trajectory) - 1)
    basic['bt'] = max(basic['bt'], 0.0)
    basic['et'] = min(basic['et'], u.trajectory.time[-1])
    basic['dt'] = max(dt_min, basic['dt'])
    if basic['dt'] > basic['et']:
        raise ParameterWrongError('dt')


def set_res(file_index):
    "Read in the residue index file."
    global res_pair
    res_pair_set = set()
    f = open(file_index, 'r')
    lines = f.readlines()
    for line in lines:
        if '$' in line:
            res_former = int(line.split("$")[0])
            res_latter_selc = []
            if 'all' in line.split("$")[0]:
                res_latter_selc = list(range(1,basic['res_num']+1))
            res_latter = line.split("$")[1].strip().split()
            for res in res_latter:
                if '-' in res:
                    start, end = res.split('-')
                    res_range = range(int(start), int(end) + 1)
                    res_latter_selc.extend(res_range)
                else:
                    res_latter_selc.append(int(res))
            res_latter_selc = list(set(res_latter_selc))
            for res in res_latter_selc:
                res_pair_set.add(tuple(sorted([res_former, res])))

        else:
            res_selc = []
            res_list = line.split()
            for res in res_list:
                if '-' in res:
                    start, end = res.split('-')
                    res_range = range(int(start), int(end) + 1)
                    res_selc.extend(res_range)
                else:
                    res_selc.append(int(res))
            res_selc = sorted(res_selc)
            combinations = itertools.combinations(res_selc, 2)
            res_pair_set.update(combinations)
    res_pair = list(res_pair_set)


def res_build():
    "Construct and save the three-dimensional adjacency matrix of the results."
    global dt_min
    t1 = timeit.default_timer()
    bt = int(basic['bt']/dt_min)
    et = int(10*basic['et']/dt_min)
    dt = int(10*basic['dt']/dt_min)
    pdbbase = basic['top']
    op = basic['res']
    fi = open(pdbbase, 'r')
    all_lines = fi.readlines()
    fi.close()
    atom_lines = [l for l in all_lines if l[0:6] == 'ATOM  ']
    for line in atom_lines:
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].replace(' ', '_')
        res_name = line[17:20]
        res_num = int(line[22:26].strip())
        chain_id = line[21:22]
        if chain_id == ' ':
            chain_id = '#'  # default:$
        res = chain_id + '+' + str(res_num) + '_' + res_name
        if res not in contact_score:
            res_num_dict[res_num] = res
            contact_score[res] = {}
            res_list.append(res)
    if op:
        for pair in res_pair:
            init_dict(res_num_dict[pair[0]], res_num_dict[pair[1]], bt, et, dt)
    else:
        for item in contact_score.keys():
            for index in res_list:
                if not (item == index) and item not in contact_score[index]:
                    init_dict(item, index, bt, et, dt)

    t2 = timeit.default_timer()
    print("build time:%10.5f:" % (t2-t1))


def init_dict(a_res, b_res, bt, et, dt):
    contact_score[a_res][b_res] = {}
    for t in range(bt, et + 1, dt):
        contact_score[a_res][b_res][t] = 0.0


def res_con2():
    "Calculate whether there is contact between residue pairs."
    global u,dt_min
    t3 = timeit.default_timer()
    bt = basic['bt']
    et = basic['et']
    dt = basic['dt']
    d_min = basic['r_min']
    d_max = basic['r_max']
    count = 1
    for ts in u.trajectory[int(bt/dt_min):int(et/dt_min)+1:int(dt/dt_min)]:
        count +=1
        time = u.trajectory.time
        if time % 50000 == 0:
            t4 = timeit.default_timer()
            print(time,f"这{count}次计算用时",t4-t3,'s')
            count = 0
            t3=t4
        res_contact = {}
        for ires in res_list:
            ires_id = int(ires.split('+')[1].split('_')[0].strip())
            ires_contact = []
            ires_atom_index = list(u.select_atoms(f'resid {ires_id} and not name H*').ids - 1)
            '''
            ires_atom_coord = u.atoms[ires_atom_index].positions
            '''
            ires_atom_name = list(u.atoms[ires_atom_index].groupby(['ids', 'names', 'occupancies']))
            ires_atom = {i[1]: [u.atoms[i[0] - 1].position, i[2]] for i in
                         ires_atom_name}  # {atom_name:[array[x,y,z],occ]}
            for kres in contact_score[ires].keys():
                kres_id = int(kres.split('+')[1].split('_')[0].strip())
                kres_atom_index = list(u.select_atoms(f'resid {kres_id} and not name H*').ids - 1)
                kres_atom_name = list(u.atoms[kres_atom_index].groupby(['ids', 'names', 'occupancies']))
                kres_atom = {i[1]: [u.atoms[i[0] - 1].position, i[2]] for i in kres_atom_name}
                kres_flag = 0
                for iatom in ires_atom.keys():
                    #print(ires_atom[iatom])
                    (ix, iy, iz), iocc = ires_atom[iatom]
                    for katom in kres_atom.keys():
                        (kx, ky, kz), kocc = kres_atom[katom]
                        dx = abs(ix - kx)
                        dy = abs(iy - ky)
                        dz = abs(iz - kz)
                        if dx < 4.14 and dy < 4.14 and dz < 4.14:
                            kres_flag = 1
                            break
                    if kres_flag:
                        ires_contact.append(kres)
                        break
            for jres in ires_contact:
                jres_id = int(jres.split('+')[1].split('_')[0].strip())
                jres_atom_index = list(u.select_atoms(f'resid {jres_id} and not name H*').ids - 1)
                jres_atom_name = list(u.atoms[jres_atom_index].groupby(['ids', 'names', 'occupancies']))
                jres_atom = {i[1]: [u.atoms[i[0] - 1].position, i[2]] for i in jres_atom_name}
                total_score = 0
                ires2np = []
                jres2np = []
                sub_heavy_atom = (abs(ires_id - jres_id) < 5)
                for iatom in ires_atom.keys():
                    if sub_heavy_atom and (iatom in ['N', 'CA', 'C', 'O']):
                        continue
                    ires2np.append(ires_atom[iatom][0] * ires_atom[iatom][1])  # position * occ
                for jatom in jres_atom.keys():
                    if sub_heavy_atom and (jatom in ['N', 'CA', 'C', 'O']):
                        continue
                    jres2np.append(jres_atom[jatom][0] * jres_atom[jatom][1])
                if not jres2np:
                    contact_score[ires][jres][time] = total_score
                    continue
                ires_np = np.array(ires2np)
                jres_np = np.array(jres2np)
                dis_np = compute_distances(ires_np, jres_np)
                total_score = rrcs(dis_np, d_max, d_min)
                contact_score[ires][jres][time] = total_score


def rrcs(mat, d_max, d_min):
    "Calculate rrcs."
    rows, cols = mat.shape
    total_score = 0
    for i in range(rows):
        for j in range(cols):
            dis = mat[i, j]
            if dis >= d_max:  # 4.63*4.63 = 21.4369
                score = 0
            elif dis <= d_min:  # 3.23*3.23 = 10.4329
                score = 1.0
            else:
                score = (1 - (dis - d_min) / 1.4)
            total_score += score
    return total_score

def compute_distances(A, B):
    "Calculate the Euclidean distance between each pair of vectors between two matrices."
    m = np.shape(A)[0]
    n = np.shape(B)[0]
    M = np.dot(A, B.T)
    H = np.tile(np.matrix(np.square(A).sum(axis=1)).T,(1,n))
    K = np.tile(np.matrix(np.square(B).sum(axis=1)),(m,1))
    return np.sqrt(-2*M+H+K)

def save_rrcs():
    "Save and output the results."
    global u
    for pair in res_pair:
        ares = res_num_dict[pair[0]]
        bres = res_num_dict[pair[1]]
        x = np.array(list(contact_score[ares][bres].keys()))
        y = np.array(list(contact_score[ares][bres].values()))
        title = f"{ares} and {bres} RRCS vs Time"
        "create xvg"
        with open(f'{ares}&{bres}.xvg', 'w') as f:
            # 写入文件头
            f.write('# Created by Python script rrcs3_1.py\n')
            f.write(f'@    title "{title}"\n')
            f.write('@    xaxis  label "Time (ps)"\n')
            f.write('@    yaxis  label "RRCS"\n')
            f.write('@TYPE xy\n')
            for i in range(len(x)):
                f.write('{:.3f}\t{:.3f}\n'.format(x[i], y[i]))
        "create png"
        plt.plot(x, y)
        plt.title(title)
        plt.xlabel("Time (ps)")
        plt.ylabel("RRCS")
        #plt.show()
        wd = os.getcwd()
        plt.savefig(f"{ares}&{bres}.png")

def chge(titl,count):
    "If there is an existing folder with the same name in the current directory, change the name of the existing folder."
    if os.path.exists(titl+str(count)):
        count+=1
        chge(titl,count)
        count-=1
    print('The dictionary '+titl+str(count-1)+' is changed to '+titl+str(count-1))
    os.rename(titl+str(count-1), titl+str(count))

try:
    if len(sys.argv) != 2:
        raise InputFileError
    input_file = sys.argv[1]
    input_para(input_file)
    check()
except InputFileError as e:
    print(f"{e.args} doesn't exist.")
except ParameterWrongError as e:
    print(f"The input parameters {e.args} do not meet the specification.")
else:
    print("Everything goes well.")
    res_build()
    tile = basic['traj'][:-4]
    res_con2()
    if not os.path.exists(tile):
        os.mkdir(tile)
    else:
        os.rename(tile, tile + '_0')
        chge(tile + '_', 1)
        os.mkdir(tile)
    os.chdir(tile)
    save_rrcs()