
import os
import numpy as np 
from aqme.qprep import qprep

destination_fukui = 'fukui_nbo'
destination_redox = 'redox_opt'
nprocs = 16
mem = '64GB'
qm_end = '$nbo NBOSUM BNDIDX $end'
qm_input = 'pbe1pbe/def2tzvp emp=gd3bj scrf=(smd,solvent=water) pop=(nbo6read)' #match other spc calculations
qm_input_redox = 'pbe1pbe/6-31+g(d) emp=gd3bj scrf=(smd,solvent=water) opt freq=noraman' #match other opt calculations
path = os.getcwd()

files = [file for file in sorted(os.listdir(path)) if file.endswith('opt.log')]

for file in files:
    file_name = os.path.basename(file).split(".")[0]
    output = open(file, 'r').readlines()
    for line in output:
        if 'Charge' and 'Multiplicity' in line:
            charge = int(line.strip().split()[2])
            mult = int(line.strip().split()[5])
            if mult == 1:
                new_mult = 2
            if mult == 2:
                new_mult = 1
            #qprep(destination=destination, files=file, program='gaussian', qm_input=qm_input, mem=mem, nprocs=nprocs, suffix= 'nbo', qm_end=qm_end, charge=charge, mult=mult)
            qprep(destination=destination_fukui, files=file, program='gaussian', qm_input=qm_input, mem=mem, nprocs=nprocs, suffix= '+1e_nbo', qm_end=qm_end, charge=charge-1, mult=new_mult) #spc
            qprep(destination=destination_fukui, files=file, program='gaussian', qm_input=qm_input, mem=mem, nprocs=nprocs, suffix= '-1e_nbo', qm_end=qm_end, charge=charge+1, mult=new_mult) #spc
            qprep(destination=destination_redox, files=file, program='gaussian', qm_input=qm_input_redox, mem=mem, nprocs=nprocs, suffix= '+1e_red', charge=charge-1, mult=new_mult) #opt, spc later
            qprep(destination=destination_redox, files=file, program='gaussian', qm_input=qm_input_redox, mem=mem, nprocs=nprocs, suffix= '-1e_ox', charge=charge+1, mult=new_mult) #opt, spc later

            break



