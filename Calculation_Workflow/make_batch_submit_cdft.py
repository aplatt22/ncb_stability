### must have cdft.py file in same directory ###
### edit red_suffix and ox_suffix appropriately ###
### code assumes original molecule has same base filename as redox files, but without +/1e. Ends with _nbo
### run in directory with nbo output files: python make_batch_submit ###

import os

original_suffix = '_spc_nbo.log'
red_suffix = '_+1e_nbo.log'
ox_suffix = '_-1e_nbo.log'

path = os.getcwd()

files_nbo = [file for file in sorted(os.listdir(path)) if file.endswith(original_suffix)]
files_red = [file for file in sorted(os.listdir(path)) if file.endswith(red_suffix)]
files_ox = [file for file in sorted(os.listdir(path)) if file.endswith(ox_suffix)]

base = []
for file in files_nbo:
    check = []
    file_base = file.split(original_suffix)[0]
    temp_files = [x for x in sorted(os.listdir(path)) if x.startswith(file_base) and x.endswith('.log')]
    for temp in temp_files:
        c = open(temp,'r')
        lines = c.readlines()
        if 'Normal termination' in lines[-1]:
            check.append('pass')
        else:
            check.append('fail')
            print('Abnormal termination of:',temp, '... Cannot calculate fukui indices.')
    if 'fail' in check:
        continue
    else:
        base.append(file_base)

f = open('batch_submit_cdft.sh', 'w+')

for conf in base:
    line = 'python cdft.py -f {}'.format(conf)+original_suffix+' --red {}'.format(conf)+red_suffix+' --ox {}'.format(conf)+ox_suffix+' > {}_fukui_indices.out'.format(conf)
    f.writelines(line)
    f.writelines('\n')
    f.writelines('echo -e "command line: {}" >> {}_fukui_indices.out'.format(line,conf))
    f.writelines('\n')
    f.writelines('mv fukui.csv {}'.format(conf)+'_fukui_indices.csv')
    f.writelines('\n')
    f.writelines('\n')
f.close()        



