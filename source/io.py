# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 09:05:49 2018

@author: Amy
"""
# =============================================================================
# Begins *** io.py ***
# 
# Writes and reads trajectories from text files 
#
# =============================================================================
def write_to_file (conf_list, name ='traj.txt'):
    """
    writes conformation list to text file
    
    conf_list is list of strings representing the conformation trajectory 
    name is a string for the the file name
    
    """
   
    with open(name, 'w') as f:
        for i in conf_list:
            if type(i) == 'str':
                f.write(i+'\n')
            else: 
                f.write(str(i) + '\n')

def read_from_file(name = 'traj.txt'):
    """
    reads conformation or native contact trajectory from a text file
    returns a list to of integers for native traj or strings for conf traj
    
    name is a string for the the file name to be read 
    """
    conf_traj = []
    with open (name, 'r') as f:
        if 'NATIVES' in name:
            for i in f:
                x = int(i.strip())
                conf_traj.append(x)
        else: 
            for i in f:
                x = i.strip()
                conf_traj.append(x)
            
    return conf_traj