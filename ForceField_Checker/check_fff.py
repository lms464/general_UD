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
    for fi, ff in enumerate(bonds_ff):
        if len(ff)==0:
            continue
        # if fi == 7431:
        #     print(ff)
        
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
    atom_map = atom_map.set_index('atom')
    atom_map['residue'] = atom_map.index
    atom_map = atom_map.rename(columns={'residue':'atom'})
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

def get_dihyd_itp(inputfile):
    f = open(inputfile,'r')
    lines = f.readlines()
    f.close()
    
    check = 0
    bonds = []
    for l in lines:
        
        if l.startswith("[ angles ]"):
            check = 0
            pass
        elif l.startswith("[ pairs ]"):
            check = 0
            pass
        elif l.startswith('[ bonds ]'):
            check = 0
        elif l.startswith('[ dihedrals ]'):
            check = 1
            continue
        tmp = l.split()
        bond_tmp = []
        if check == 1:
            for t in tmp:
                if t.startswith("position") or t.startswith("#"):
                    check = 0
                    break
                elif t == ' ' or t.startswith(';') or t.startswith(']') or t.startswith('[') :
                    continue
                else:
                    bond_tmp.append(t)
            bonds.append(bond_tmp)
    bonds = bonds[1:-2]
    bonds_map = pd.DataFrame(bonds,columns=(["ind1","ind2",'ind3','ind4',"bond"]))
    return bonds_map    

#######################################
#### Merge bond/angles/... with atomtypes
#########################################

def rewrite_bonds(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["atom"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["atom"].values[0]
    return update_bonds

def rewrite_angle(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj,bk) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'],inptBonds['ind3'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["atom"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["atom"].values[0]
        update_bonds["ind3"][k] = inputitp.query("ind == '%s'"%bk)["atom"].values[0]
    return update_bonds

def rewrite_dihedral(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj,bk,bl) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'],inptBonds['ind3'],inptBonds['ind4'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["atom"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["atom"].values[0]
        update_bonds["ind3"][k] = inputitp.query("ind == '%s'"%bk)["atom"].values[0]
        update_bonds["ind4"][k] = inputitp.query("ind == '%s'"%bl)["atom"].values[0]
    return update_bonds

#####################################
#### Compare bonds/angles/... to FF
#####################################

# for now ff_bonds1/2 are a list of atomtypes you want to compare
def check_force_fields_bonds(ff_bond1, ff_bonds):
    
    def bondcompareB (inpt1,inpt_ref,bond_list):
        if (inpt1[0] == inpt_ref[0] and inpt1[1] == inpt_ref[1]) or (inpt1[1] == inpt_ref[0] and inpt1[0] == inpt_ref[1]):
            bond_list.append(inpt_ref)
        return bond_list
        
    bond = []
    for fi in ff_bonds:
        bond = bondcompareB(ff_bond1, fi,bond)
    return bond


def compare_bonds(bond_in,bond_1,mol1):
    # ang_dpop = np.loadtxt("dpop.ang",dtype=str)
    # ang_lano = np.loadtxt('lano.ang',dtype=str)
    bonds_ = []
    for b in bond_1.values:
        bonds_.append(check_force_fields_bonds(mol1.T[b[:-1]].T['AtomType'].tolist(),bond_in))
    for ai, a in enumerate(bonds_):
        if len(a)<2:
            continue
        else:
            print(a[0][5:],a[1][5:])
    bond_final = []    
    for d,d1 in zip(bonds_,bond_1.values):
        if len(d)>1:
            bond_final.append([d1[:-1].tolist(),d[0][-3:]])
        elif len(d)<1:
            continue
        else:
            bond_final.append([d1[:-1].tolist(),d[0][-3:]])
    return bond_final


def check_force_fields_angles(ff_ang1, ff_angs):
    
    def bondcompareA (inpt1,inpt_ref,_list):
        # for a in permutations(inpt1):
        if (inpt1[0] == inpt_ref[0] and inpt1[1] == inpt_ref[1] and inpt1[2] == inpt_ref[2]) or (inpt1[2] == inpt_ref[0] and inpt1[1] == inpt_ref[1] and inpt1[0] == inpt_ref[2]):
            _list.append(inpt_ref)
        return _list
        
    angs = []

    for fi in ff_angs:
        angs = bondcompareA(ff_ang1, fi,angs)
    return angs

def check_force_fields_dihedral(ff_ang1, ff_angs):
    
    def bondcompareD (inpt1,inpt_ref,D_list):
        # for a in permutations(inpt1):
        if (inpt1[0] == inpt_ref[0] and inpt1[1] == inpt_ref[1] and inpt1[2] == inpt_ref[2] and inpt1[3] == inpt_ref[3]) or (inpt1[3] == inpt_ref[0] and inpt1[2] == inpt_ref[1] and inpt1[1] == inpt_ref[2] and inpt1[0] == inpt_ref[3]):
           D_list.append(inpt_ref)
        return D_list
        
    angs = []
    for fi in ff_angs:
        angs = bondcompareD(ff_ang1, fi,angs)
    return angs

def compare_angles(angles_in,ang_1,mol1):
    # ang_dpop = np.loadtxt("dpop.ang",dtype=str)
    # ang_lano = np.loadtxt('lano.ang',dtype=str)
    angs = []
    for ad in ang_1.values:
        angs.append(check_force_fields_angles(mol1.T[ad[:-1]].T['AtomType'].tolist(),angles_in))
    for ai, a in enumerate(angs):
        if len(a)<2:
            continue
        else:
            print(a[0][5:],a[1][5:])
    ang_final = []    
    for d,d1 in zip(angs,ang_1.values):
        if len(d)>1:
            ang_final.append([d1[:-1].tolist(),d[0][-3:]])
        elif len(d)<1:
            continue
        else:
            ang_final.append([d1[:-1].tolist(),d[0][-3:]])
    return ang_final
            
def compare_dihydrals(dyhedral_in,dih_1, mol1): #, dih_2, mol2):
    # dih_dpop = np.loadtxt("dpop.dih",dtype=str)
    # dih_lano = np.loadtxt('lano.dih',dtype=str)
    dihds = []
    for ai, ad in enumerate(dih_1.values):
        dihds.append(check_force_fields_dihedral(mol1.T[ad[:-1]].T['AtomType'].tolist(),dyhedral_in)) #,mol2.T[al[:-1]].T['AtomType'].tolist(),dyhedral_in))
    dihds_final = []
    for d,d1 in zip(dihds,dih_1.values):
        if len(d)>1:
            dihds_final.append([d1[:-1].tolist(),d[0][-3:]])
        elif len(d)<1:
            continue
        else:
            dihds_final.append([d1[:-1].tolist(),d[0][-3:]])
    return dihds_final
            
###############
### RUN
##############
bonded_ff = parse_ff("ffbonded.itp")#/home/liam/Downloads/charmm-gui-5955053520/gromacs/toppar/forcefield.itp")
    #"/home/sharplm/Programs/gromacs-2021.4/share/top/charmm36-jul2021.ff/ffbonded.itp")#(

Lanterol =  get_atom_map('LANO.itp')
lanterol_bnds = get_bonds_itp('LANO.itp')
# lanterol_dih = get_dihyd_itp('LANO.itp')
lanterol_bnds = rewrite_bonds(lanterol_bnds,Lanterol)
# lanterol_dih = rewrite_dihedral(lanterol_dih,Lanterol)
# lanterol_ang = get_angles_itp("LANO.itp")
# lanterol_ang = rewrite_angle(lanterol_ang, Lanterol)
# # carbon_l = Lanterol["atom"].str.startswith("C")[Lanterol["atom"].str.startswith("C")].index.values

ff_dih = parse_ff_dihedral(bonded_ff)
ff_ang = parse_ff_angles(bonded_ff)
ff_bond = parse_ff_bonds(bonded_ff)

diplop = get_atom_map('DPOP.itp')
# diplop_dih = get_dihyd_itp('DPOP.itp')
# diplop_dih = rewrite_dihedral(diplop_dih,diplop)


diplop_bonds = get_bonds_itp('DPOP.itp')
diplop_bonds = rewrite_bonds(diplop_bonds,diplop)
dip_bond = compare_bonds(ff_bond,diplop_bonds,diplop)
lan_bond = compare_bonds(ff_bond,lanterol_bnds,Lanterol)

# diplop_ang = get_angles_itp('DPOP.itp')
# diplop_ang = rewrite_angle(diplop_ang,diplop)

# dpl_d = compare_dihydrals(ff_dih,diplop_dih,diplop)
# lan_a = compare_angles(ff_ang,lanterol_ang,Lanterol)
# dip_a = compare_angles(ff_ang,diplop_ang,diplop)

# lan_d = compare_dihydrals(ff_dih,lanterol_dih,Lanterol)

dpl_di = open("Bond_DPOP.bnd",'w')
for d in dip_bond:
    tmp_txt = ''
    for d1 in d[0]:
        tmp_txt = tmp_txt+d1+' '
    for d1 in d[1]:
        tmp_txt = tmp_txt+d1+' '
    tmp_txt = tmp_txt+'\n'
    dpl_di.write(tmp_txt)
dpl_di.close()

lan_di = open("Bond_LANO.bnd",'w')
for d in lan_bond:
    tmp_txt = ''
    for d1 in d[0]:
        tmp_txt = tmp_txt+d1+' '
    for d1 in d[1]:
        tmp_txt = tmp_txt+d1+' '
    tmp_txt = tmp_txt+'\n'
    lan_di.write(tmp_txt)
lan_di.close()

# print("angles")
# compare_angles(ff_ang,diplop_ang,diplop_ang)

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