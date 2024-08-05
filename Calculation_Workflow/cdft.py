#!/usr/bin/python
from __future__ import print_function, absolute_import

#Python Libraries
import math, sys, os
from glob import glob
from optparse import OptionParser
import pandas as pd
import numpy as np
from rdkit import Chem
import cclib as cc
import dbstep.Dbstep as db

hartree_to_ev = 27.211396641308

elements = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
            "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni",
            "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo",
            "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
            "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
            "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
            "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
            "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq",
            "Uup","Uuh","Uus","Uuo"]

## covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197 ##
## values for metals decreased by 10 % ##
rcov = [0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
        1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,
        1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
        1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,
        1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,
        1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,
        1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
        1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
        1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,
        1.52, 1.53, 1.54, 1.55]

#Some useful arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
    "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

# Calculation of atomic coordination numbers from dftd3
def ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2):
   cn =[]
   for i in range(0,natom):
      xn = 0.0
      for iat in range(0,natom):
         if iat != i:
            dx, dy, dz = xco[iat] - xco[i], yco[iat] - yco[i], zco[iat] - zco[i]
            r = math.pow(dx*dx+dy*dy+dz*dz,0.5)

            for k in range(0,max_elem):
               if atomtype[i].find(elements[k])>-1:Zi=k
               if atomtype[iat].find(elements[k])>-1:Ziat=k

            rco = k2 * (rcov[Zi]+rcov[Ziat])
            rr=rco/r
            damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
            xn=xn+damp
      cn.append(xn)
   return cn

def parse_cc_data(calc_dir, filename):
    ## Use cclib to parse QM data - unfortunately doesn't parse natural populations correctly for radicals
    file = calc_dir+'/'+filename
    parser = cc.io.ccopen(file)
    try:
        data = parser.parse()
        cc_data = data
    except:
        cc_data = None

    try: # grab Wiberg Bond Orders summed for each atom
        start_wiberg, end_wiberg = None, None
        outfile = open(file,"r")
        lines = outfile.readlines()
        for i,line in enumerate(lines):
            if line.find("Wiberg bond index, Totals by atom:") > -1: start_wiberg = i+4
            if line.find("NBI: Natural Binding Index (NCU strength parameters)") > -1: end_wiberg = i-2

        if start_wiberg != None and end_wiberg != None:
            wiberg_bos = []
            for i in range(start_wiberg, end_wiberg):
                wiberg_bos.append(float(lines[i].split()[2]))
            setattr(cc_data, "bondorders", wiberg_bos)
    except: pass

    v_nbo = 6
    try: # annoying but necessary due to cclib's bug
        start_npop = None
        outfile = open(file,"r")
        lines = outfile.readlines()
        for i,line in enumerate(lines):
            if line.find("-    Spin") > -1: start_npop = i+3
            if line.find("NBO Version 3") >-1: v_nbo =3
        if start_npop != None:
            nat_charges = []
            end_npop = start_npop + len(cc_data.atomnos)
            for i in range(start_npop, end_npop):
                nat_charges.append(float(lines[i].split()[2]))
            cc_data.atomcharges['natural'] = nat_charges
    except: pass

    if v_nbo == 3:
        print("   THIS SOFTWARE IS NOT COMPATIBLE WITH NBO v3!")
        sys.exit(0)
    return cc_data

def dbstep_vbur(calc_dir, filename, id):
    ## use dbstep to compute buried volume from DFT log file
    file = calc_dir+'/'+filename
    #try:
    db_data = db.dbstep(file, atom1=id, volume=True, commandline=True, quiet=True)
    #except:
    #    print('DBSTEP failed for file', file)
    if hasattr(db_data, 'bur_vol'): v_bur = db_data.bur_vol
    else: v_bur = None
    return v_bur

def generate_features(dir, file, oxidized=None, reduced=None, mulliken=False, spin=False, vbur=False, coord=False, ref_oxidized=None, ref_reduced=None):
    file_cc_data = parse_cc_data(dir, file)

    if oxidized != None:
        ox_cc_data = parse_cc_data(dir, oxidized)
    else: ox_cc_data = None

    if reduced != None:
        red_cc_data = parse_cc_data(dir, reduced)
    else: red_cc_data = None

    if ref_oxidized != None:
        ref_ox_data = parse_cc_data(dir, ref_oxidized)
    else: ref_ox_data = None

    if ref_reduced != None:
        ref_red_data = parse_cc_data(dir, ref_reduced)
    else: ref_red_data = None

    atomnos = file_cc_data.atomnos if hasattr(file_cc_data, 'atomnos') else []
    elements = [periodictable[atomno] for atomno in atomnos]
    mulliken_spins = file_cc_data.atomspins['mulliken'] if hasattr(file_cc_data, 'atomspins') else [0] * len(atomnos)
    if not mulliken: npa_charges = file_cc_data.atomcharges['natural'] if hasattr(file_cc_data, 'atomcharges') else [0] * len(atomnos)
    else: npa_charges = file_cc_data.atomcharges['mulliken'] if hasattr(file_cc_data, 'atomcharges') else [0] * len(atomnos)
    coords = file_cc_data.atomcoords[-1] if hasattr(file_cc_data, 'atomcoords') else []
    idx = [i+1 for i in range(0,len(atomnos))]

    # df of atomic descriptors
    CDFT = pd.DataFrame(list(zip(elements, idx, npa_charges)), columns =['element', 'idx', 'q_NPA'])

    if spin: CDFT['spin'] = mulliken_spins
    if hasattr(file_cc_data, 'bondorders'): CDFT['WBOs'] = file_cc_data.bondorders

    # chemical potential and hardness
    e_i, e_a, mu, eta, gei = None, None, None, None, None
    if ox_cc_data != None:
        if hasattr(ox_cc_data, 'scfenergies'): e_i = ox_cc_data.scfenergies[-1]
    if red_cc_data != None:
        if hasattr(red_cc_data, 'scfenergies'): e_a = red_cc_data.scfenergies[-1]
    if e_i != None and e_a != None:
        mu = (e_a - e_i ) / 2
        eta = (e_i - e_a) # Org. Lett. 2007, 9, 14, 2721 suggests that factor of 2 is no longer used
    if mu != None and eta != None:
        gei = (mu ** 2) / (2 * eta)

    # relative nucleophilicity
    ref_e_i, ref_e_a, mu_ref, eta_ref, gni = None, None, None, None, None
    if ref_ox_data != None:
        if hasattr(ref_ox_data, 'scfenergies'): ref_e_i = ref_ox_data.scfenergies[-1]
    if ref_red_data != None:
        if hasattr(ref_red_data, 'scfenergies'): ref_e_a = ref_red_data.scfenergies[-1]
    if ref_e_i != None and ref_e_a != None:
        mu_ref = (ref_e_a - ref_e_i ) / 2
        eta_ref = (ref_e_i - ref_e_a) # Org. Lett. 2007, 9, 14, 2721 suggests that factor of 2 is no longer used
    if None not in [mu, eta, mu_ref, eta_ref]:
        gni = (((mu - mu_ref) ** 2) / ((eta + eta_ref) ** 2)) * (eta/2)

    # natural charges
    if ox_cc_data != None:
        if not spin:
            if not mulliken: ox_npa_charges = ox_cc_data.atomcharges['natural'] if hasattr(ox_cc_data, 'atomcharges') else []
            else: ox_npa_charges = ox_cc_data.atomcharges['mulliken'] if hasattr(ox_cc_data, 'atomcharges') else []
            CDFT['ox_q_NPA'] = ox_npa_charges
        else:
            ox_mulliken_spins = ox_cc_data.atomspins['mulliken'] if hasattr(ox_cc_data, 'atomspins') else []
            CDFT['ox_spin'] = ox_mulliken_spins

    if red_cc_data != None:
        if not spin:
            if not mulliken:
                red_npa_charges = red_cc_data.atomcharges['natural'] if hasattr(red_cc_data, 'atomcharges') else []
                #print(red_cc_data.atomcharges)
            else: red_npa_charges = red_cc_data.atomcharges['mulliken'] if hasattr(red_cc_data, 'atomcharges') else []
            CDFT['red_q_NPA'] = red_npa_charges
        else:
            red_mulliken_spins = red_cc_data.atomspins['mulliken'] if hasattr(red_cc_data, 'atomspins') else []
            CDFT['red_spin'] = red_mulliken_spins

    # fukui functions
    if not spin:
        if ox_cc_data != None: CDFT['f-'] =   CDFT['ox_q_NPA'] - CDFT['q_NPA']
        if red_cc_data != None: CDFT['f+'] =   CDFT['q_NPA'] - CDFT['red_q_NPA']
        if ox_cc_data != None and red_cc_data != None: CDFT['frad'] =   0.5 * CDFT['ox_q_NPA'] - 0.5 * CDFT['red_q_NPA']

    # Parr functions
    else:
        if ox_cc_data != None: CDFT['P-'] =  CDFT['ox_spin'] - CDFT['spin']
        if red_cc_data != None: CDFT['P+'] =  CDFT['red_spin'] - CDFT['spin']
        if ox_cc_data != None and red_cc_data != None: CDFT['Prad'] =  0.5 * CDFT['red_spin'] - 0.5 * CDFT['ox_spin']

    # Local Electrophilicity Index
    if 'f+' in CDFT and gei != None: CDFT['mu+'] = CDFT['f+'] * gei
    if 'f-' in CDFT and gni != None: CDFT['mu-'] = CDFT['f-'] * gni

    #buried volumes computed for all atoms except H
    if vbur: CDFT['V_bur'] = [dbstep_vbur(dir, file, (id+1)) if atomnos[id] == 6 else np.nan for id in range(0,len(atomnos))]

    # D3 coordination numbers for each atom
    if coord:
        xco, yco, zco = [x for [x,y,z] in coords], [y for [x,y,z] in coords], [z for [x,y,z] in coords]
        CDFT['d3_cdn_no'] = ncoord(len(atomnos), rcov, elements, xco, yco, zco, 94, 0.52917726, 16.0, 4.0/3.0)

    # Sum of Wiberg bond indices around each atom

    # Add coordinates
    #CDFT['coords'] = coords

    return CDFT, gei, gni

def main():
    # get command line inputs. Use -h to list all possible arguments and default values
    parser = OptionParser(usage="Usage: %prog [options] <input1>.log <input2>.log ...")
    parser.add_option("--dir", dest="dir", action="store", help="directory containing files", default='.')
    parser.add_option("-f", dest="file", action="store", help="filename with NBO data", default=None)
    parser.add_option("--mulliken", dest="mulliken", action="store_true", help="Request Mulliken charges be used rather than natural (NPA) charges", default=False)
    parser.add_option("--red", dest="reduced", action="store", help="filename of reduced species with NBO data", default=None)
    parser.add_option("--ox", dest="oxidized", action="store", help="filename of oxidized species with NBO data", default=None)
    parser.add_option("--vbur", dest="vbur", action="store_true", help="Compute % buried volumes", default=False)
    parser.add_option("--coord", dest="coord", action="store_true", help="Compute D3 coordination numbers", default=False)
    parser.add_option("--spin", dest="spin", action="store_true", help="Compute Fukui function using Mulliken spin rather than natural population", default=False)
    parser.add_option("--dropH", dest="dropH", action="store_true", help="Drop H atoms from the final dataframe", default=False)
    parser.add_option("--ref_ox", dest="ref_ox", action="store", help="Define a reference oxidized species (a file, typically F+) used to compute relative nucleophilicity", default=None)
    parser.add_option("--ref_red", dest="ref_red", action="store", help="Define a reference reduced species (a file, typically F-) used to compute relative nucleophilicity", default=None)
    (options, args) = parser.parse_args()

    if options.file != None:
        dir = options.dir
        file = options.file
        ox = options.oxidized
        red = options.reduced
        mulliken = options.mulliken
        coord = options.coord
        spin = options.spin
        vbur = options.vbur
        ref_ox = options.ref_ox
        ref_red = options.ref_red

        CDFT, GEI, GNI = generate_features(dir, file, ox, red, mulliken, spin, vbur, coord, ref_ox, ref_red)
        if GEI != None: print('   Global Electrophilicity Index {:.3f} eV'.format(GEI))
        if GNI != None: print('   Relative Nucleophilicity Index {:.3f} eV (defined relative to {}/{})'.format(GNI, ref_ox, ref_red))

        if options.dropH == True: # drop H atoms
            h_indices = CDFT[ CDFT['element'] == 'H'].index
            CDFT.drop(h_indices, inplace = True)

        CDFT.to_csv('fukui.csv')
        print(CDFT)

    else:
        print('  A file needs to be specified!')

if __name__ == "__main__":
    main()
