#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Begins fileio.py
#
# File input and output. Read in trajectory from previous simulation.
# =============================================================================

def save(conf_list, filename="traj.txt", interval = 1):
    f = open(filename,'w+')
    for i in range(len(conf_list)):
        f.write("{}\t".format(i+1) + "{}\t".format(conf_list[i]) + "\n")
    f.close()
    
def read(filename):
    from src.energy import energy
    from proteins.chignolin import seq, native_list
    conf_list = []
    acceptance = []
    etotal = []
    conf_list = []
    native = []
    contact = []
    
    f = open(filename, 'r')
    line = f.readline()
    conf_list.append(line.split('\t')[1])
    etotal.append(energy(conf_list[0], seq, native_list)[0])
    contact.append(energy(conf_list[0], seq, native_list)[1])
    native.append(energy(conf_list[0], seq, native_list)[2])
    
    while True:
        if len(etotal)%100000 == 0:
            print(len(etotal)/100000000)
        line = f.readline()
        if line == '':
            break
        if line.split('\t')[1] == conf_list[-1]:
            conf_list.append(conf_list[-1])
            etotal.append(etotal[-1])
            contact.append(contact[-1])
            native.append(native[-1])
        else:
            conf_list.append(line.split('\t')[1])
            temp1, temp2, temp3 = energy(conf_list[-1], seq, native_list)
            etotal.append(temp1)
            contact.append(temp2)
            native.append(temp3)

    
    return acceptance, etotal, conf_list, native, contact