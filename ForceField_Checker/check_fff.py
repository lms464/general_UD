import pandas as pd

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

def rewrite_bonds(inptBonds,inputitp):
    update_bonds = inptBonds.copy()
    for k,(bi,bj) in enumerate(zip(inptBonds["ind1"],inptBonds['ind2'])):
        update_bonds["ind1"][k] = inputitp.query("ind == '%s'"%bi)["AtomType"].values[0]
        update_bonds["ind2"][k] = inputitp.query("ind == '%s'"%bj)["AtomType"].values[0]
    return update_bonds

bonded_ff = parse_ff("/usr/local/gromacs/share/gromacs/top/charmm36-jul2021.ff/ffbonded.itp")

Lanterol =  get_atom_map('/home/liam/Downloads/charmm-gui-5955053520/gromacs/toppar/LANO.itp')
lanterol_bonds = get_bonds_itp('/home/liam/Downloads/charmm-gui-5955053520/gromacs/toppar/LANO.itp')
lanterol_bonds = rewrite_bonds(lanterol_bonds,Lanterol)

diplop = get_atom_map('/home/liam/lms464/resources/Hapnoid/test_Diplop/gromacs/toppar/DPOP.itp')
diplop_bonds = get_bonds_itp('/home/liam/lms464/resources/Hapnoid/test_Diplop/gromacs/toppar/DPOP.itp')
diplop_bonds = rewrite_bonds(diplop_bonds,diplop)




                                                      