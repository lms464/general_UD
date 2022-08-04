import sys

import pandas as pd
import numpy as np
# from itertools import permutations

print(sys.argv)

####################################
## Parse through charmm36-lipid files
###################

def parse_ff(inputfile):
    f = open(inputfile,'r')
    lines = f.readlines()
    f.close()
    
    ff_lines = []
    
    for l in lines:
        
        if l.startswith(';'):
            continue
        if l.startswith('\n'):
            continue
        l_tmp = l.split(';')[0]
        ff_lines.append(l_tmp)
    
    tmp_lines = []
    for flines in ff_lines:
        tmp_split = flines.split()
        tmp = []
        for ti in tmp_split:
            if ti == ' ' or ti == '[' or ti == ']':
                continue
            else:
                tmp.append(ti)
        tmp_lines.append(tmp)
    return tmp_lines

def parse_ff_bonds(bonds_ff):
    bond_ = []
    for fi,ff in enumerate(bonds_ff):
        if ff[0] == 'angletypes':
            break
        elif ff[0] == 'bondtypes':
            continue
        bond_.append(ff)
    bond_map =  bond_map = np.array(bond_)#pd.DataFrame(bond_,columns=(['atom1','atom2','bondtype','b0','kb']))
    return bond_map

def parse_ff_angles(bonds_ff):
    bond_ = []
    checker = 0
    for ff in bonds_ff:
        
        if ff[0] == 'angletypes':
            checker = 1
            continue
        elif ff[0] == 'bondtypes':
            pass
        elif ff[0] == 'dihedraltypes':
            break
              
        if checker == 1:
            bond_.append(ff)
        elif checker == 0:
            pass
    bond_map =  bond_map = np.array(bond_)#pd.DataFrame(bond_,columns=(['atom1','atom2','atom3','bondtype','theta0','ktheta','ubo','ukb']))
    return bond_map

def parse_ff_dihedral(bonds_ff):
    bond_ = []
    checker = 0
    for ff in bonds_ff:
        
        if ff[0] == 'angletypes':
            pass
        elif ff[0] == 'bondtypes':
            pass
        elif ff[0] == 'dihedraltypes':
            checker = 1
            continue
              
        if checker == 1:
            bond_.append(ff)
        elif checker == 0:
            pass
    bond_map = np.array(bond_)#pd.DataFrame(bond_,columns=(['atom1','atom2','atom3','atom4','bondtype','phi0','kphi','mult']))
    return bond_map


######################################
#### Parse info from molecule itp's
######################################
def get_atom_map(inputfile):
    f = open(inputfile,'r')
    lines = f.readlines()
    f.close()
    
    atoms = []
    for l in lines:
        line_tmp= []
        if l.startswith("[ bonds ]"):
            break
        if l.startswith(';'):
            continue
        if l.startswith('\n'):
            continue
        tmp = l.split(';')[0].split()
        for t in tmp:
            if t == ' ': # or t == '[' or t == ']':
                continue
            else:
                line_tmp.append(t)
        atoms.append(line_tmp[1:])
    atoms = atoms[3:]
    atom_map = pd.DataFrame(atoms, columns=(["AtomType","residue","resname","atom",'ind','q','m']))
    return atom_map

def get_bonds_itp(inputfile):
    f = open(inputfile,'r')
    lines = f.readlines()
    f.close()
    
    check = 0
    bonds = []
    for l in lines:
        if l.startswith("[ pairs ]"):
            break
        if l.startswith('[ bonds ]') == False and check==0:
            continue
        elif l.startswith('[ bonds ]'):
            check = 1
        tmp = l.split()
        bond_tmp = []
        for t in tmp:
            if t == ' ' or t.startswith(';') or t.startswith(']') or t.startswith('[') :
                continue
            else:
                bond_tmp.append(t)
        bonds.append(bond_tmp)
    bonds = bonds[2:-1]
    bonds_map = pd.DataFrame(bonds,columns=(["ind1","ind2","bond"]))
    return bonds_map    

def get_angles_itp(inputfile):
    f = open(inputfile,'r')
    lines = f.readlines()
    f.close()
    
    check = 0
    bonds = []
    for l in lines:
        
        if l.startswith("[ angles ]"):
            check = 1
            continue
        elif l.startswith("[ pairs ]"):
            check = 0
            pass
        elif l.startswith('[ bonds ]'):
            check = 0
        elif l.startswith('[ dihedrals ]'):
            break
        tmp = l.split()
        bond_tmp = []
        if check == 1:
            for t in tmp:
                if t == ' ' or t.startswith(';') or t.startswith(']') or t.startswith('[') :
                    continue
                else:
                    bond_tmp.append(t)
            bonds.append(bond_tmp)
    bonds = bonds[1:-1]
    bonds_map = pd.DataFrame(bonds,columns=(["ind1","ind2",'ind3',"bond"]))
    return bonds_map    

#######################################
#### Merge bond/angles/... with atomtypes
#########################################

def rewrite_bonds(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["AtomType"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["AtomType"].values[0]
    return update_bonds

def rewrite_angle(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj,bk) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'],inptBonds['ind3'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["AtomType"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["AtomType"].values[0]
        update_bonds["ind3"][k] = inputitp.query("ind == '%s'"%bk)["AtomType"].values[0]
    return update_bonds


#####################################
#### Compare bonds/angles/... to FF
#####################################

# for now ff_bonds1/2 are a list of atomtypes you want to compare
def check_force_fields_bonds(ff_bond1, ff_bond2, ff_bonds):
    
    def bondcompareB (inpt1,inpt_ref,bond_list):
        if (inpt1[0] == inpt_ref[0] and inpt1[1] == inpt_ref[1]) or (inpt1[1] == inpt_ref[0] and inpt1[0] == inpt_ref[1]):
            bond_list.append(inpt_ref)
        return bond_list
        
    bond = []
    for b1 in[ff_bond1,ff_bond2]:
        for fi in ff_bond:
            bond = bondcompareB(b1, fi,bond)
    return bond

def check_force_fields_angles(ff_ang1, ff_ang2, ff_angs):
    
    def bondcompareA (inpt1,inpt_ref,_list):
        # for a in permutations(inpt1):
        if (inpt1[0] == inpt_ref[0] and inpt1[1] == inpt_ref[1] and inpt1[2] == inpt_ref[2]) or (inpt1[2] == inpt_ref[0] and inpt1[1] == inpt_ref[1] and inpt1[0] == inpt_ref[2]):
            _list.append(inpt_ref)
        return _list
        
    angs = []
    for a1 in[ff_ang1,ff_ang2]:
        for fi in ff_angs:
            angs = bondcompareA(a1, fi,angs)
    return angs

def compare_angles(angles_in):
    ang_dpop = np.loadtxt("dpop.ang",dtype=str)
    ang_lano = np.loadtxt('lano.ang',dtype=str)
    angs = []
    for ad, al in zip(ang_dpop,ang_lano):
        angs.append(check_force_fields_angles(ad.tolist(),al.tolist(),angles_in))
    for a in angs:
        if len(a)<2:
            continue
        else:
            print(a[0][5:],a[1][5:])
            
            
###############
### RUN
##############
bonded_ff = parse_ff("/usr/local/gromacs/share/gromacs/top/charmm36-jul2021.ff/ffbonded.itp")
    #"/home/sharplm/Programs/gromacs-2021.4/share/top/charmm36-jul2021.ff/ffbonded.itp")#(

# Lanterol =  get_atom_map('/home/liam/Downloads/charmm-gui-5955053520/gromacs/toppar/LANO.itp')
# lanterol_bonds = get_bonds_itp('/home/liam/Downloads/charmm-gui-5955053520/gromacs/toppar/LANO.itp')
# lanterol_bonds = rewrite_bonds(lanterol_bonds,Lanterol)
ff_dih = parse_ff_dihedral(bonded_ff)
ff_ang = parse_ff_angles(bonded_ff)
ff_bond = parse_ff_bonds(bonded_ff)
compare_angles(ff_ang)
# diplop = get_atom_map('/home/sharplm/resources/Hapnoid/test_Diplop/gromacs/toppar/DPOP.itp')
# diplop_bonds = get_bonds_itp('/home/sharplm/resources/Hapnoid/test_Diplop/gromacs/toppar/DPOP.itp')
# diplop_bonds = rewrite_bonds(diplop_bonds,diplop)
# diplop_ang = get_angles_itp('/home/sharplm/resources/Hapnoid/test_Diplop/gromacs/toppar/DPOP.itp')
# diplop_ang = rewrite_angle(diplop_ang,diplop)



# if __name__ == "__main__":
# in1,in2= sys.argv[0],sys.argv[1]
# tmp1,tmp2 = check_force_fields_angles(in1,in2,ff_ang)
# print(tmp1)
# print(tmp2)

# for bi in diplop_bonds.values:
#     for fi in ff_bond:
#         if (bi[0] == fi[0] and bi[1] == fi[1]) or (bi[1] == fi[0] and bi[0] == fi[1]):
#             print(fi)
# import itertools
         
# for bi in diplop_ang.values:
#     for fi in ff_ang:
#         for a in itertools.permutations(bi[:3]):
#             if (a[0] == fi[0] and a[1] == fi[1] and a[2] == fi[2]):
#                 print(fi)
#                 continue