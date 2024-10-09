#!/bin/bash

##################################################
################    MODULES     ##################
##################################################

source $HOME/.bash_profile

source activate NCBs # activate your conda environment with AQME, OpenBabel, CREST, and CENSO

##################################################
###############    MAIN SCRIPT     ###############
##################################################


# Run CSEARCH conformational search
python -m aqme --csearch --program crest --input 'generated_ncbs.csv' --crest_keywords '-T 1 --ewin 12 -noreftopo' --cregen True --cregen_keywords '-T 1 --ethr 0.1 --rthr 0.125 --bthr 0.1'

# Use obabel to convert sdf files to xyz for censo
cd CSEARCH/
for f in *.sdf; 
	do obabel -isdf $f -oxyz -O ${f%.*}.xyz; 
done

# Running CENSO calculations in a new folder
cd CSEARCH/
mkdir censo_calcs/
mv *_crest.xyz censo_calcs/
cd censo_calcs/
 while read line; do
        filename=$(echo $line | awk '{print($1)}');
        charge=$(echo $line | awk '{print($3)}');
        unpaired=$(echo $line | awk '{print($5)}');
       censo -inp ${filename} -chrg ${charge} -u ${unpaired} -P 48 -O 1 -part0 on -part1 on -part2 off > ${filename%_*}_censo.out;
        mv enso_ensemble_part1.xyz ${filename%_*}_censo_ensemble.xyz;
	rm coord*;
	rm *.json;
	rm *.dat;
 	rm -r CONF*/
done < ../../censo_assist.csv
mv *_censo_ensemble.xyz ../
cd ../../

# Running CREGEN to cluster the conformers
cd CSEARCH/
while read line; do
	crest_file=$(echo $line | awk '{print($6)}');
	cregen_file=$(echo $line | awk '{print($2)}');
	crest censo_calcs/${crest_file} --cregen ${cregen_file} --cluster 10 > ${cregen_file%.*}_clustering.out;
	mv crest_clustered.xyz ${cregen_file%.*}_clustered.xyz;
	rm anmr*;
	rm cre*;
	rm *coord*;
	rm *.xyz.sorted;
	rm struc.xyz;
	rm cluster.order;
done < ../censo_assist.csv
cd ../

# AQME QPREP to prepare gaussian input files
while read line; do
	charge=$(echo $line | awk '{print($3)}');
	mult=$(echo $line | awk '{print($4)}');
	filename=$(echo $line | awk '{print($7)}');
	python -m aqme --qprep --suffix optimization --charge ${charge} --mult ${mult} --files 'CSEARCH/'${filename} --qm_input 'pbe1pbe/6-31+g(d) emp=GD3BJ scrf=(smd,solvent=water) opt freq' --program gaussian --mem 64GB --nprocs 16;
done < censo_assist.csv

# Submitting Gaussian Jobs
for f in QCALC/*.com; do /usr/local/Gaussian/G16C/g16/g16 $f; done

# Checking the Gaussian Jobs
# Run first QCORR analysis
python -m aqme --qcorr --files "QCALC/*log" --mem '64GB' --nprocs 16 --freq_conv 'opt=(calcfc,maxstep=5)' --isom_type 'com' --isom_inputs "QCALC"

# Corrections
if [ -d "QCALC/failed/run_1/fixed_QM_inputs" ]; then
    echo "Starting one rerun of failed calculations"
    for f in QCALC/failed/run_1/fixed_QM_inputs/*.com; do /usr/local/Gaussian/G16C/g16/g16 $f; done
    python -m aqme --qcorr --files "QCALC/failed/run_1/fixed_QM_inputs/*.log" --fullcheck --mem '64GB' --nprocs 16 --isom_type 'com' --isom_inputs "QCALC/failed/run_1/fixed_QM_inputs/*.log"
    else
            echo "No fixable errors found"
    fi

echo "Done with opt-freq calculations"

