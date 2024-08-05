#!/bin/bash

# prepare nbo single point calculations
python -m aqme --qprep --program gaussian --files "*.log" --destination "sp_nbo" --qm_input 'pbe1pbe/def2tzvp emp=gd3bj scrf=(smd,solvent=water) pop=(nbo6read)' --qm_end '$nbo NBOSUM BNDIDX $end' --mem '64GB' --nprocs 16 --suffix 'spc_nbo'

# prepare Fukui index calculations
python fukui_nbo_input.py
