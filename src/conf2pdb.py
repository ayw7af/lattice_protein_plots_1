#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Begins conf2pdb.py
#
# Function that converts a conformation trajectory to its pdb format.
# =============================================================================

from src.AAdict import *

def conf2pdb(seq, conf, frequency = 1):
    ofilename = "latticeprotein.pdb"
    ofile = open(ofilename, "a+")
    x = []
    y = []
    z = []
    
    # Lattice constant in \AA.
    lattice_const = 1.4
    
    for i in range(len(conf) + 1):
        x.append(0.0)
        y.append(0.0)
        z.append(0.0)
        
    # print abcd
    for i in range((len(conf)//2)):
        if conf[i]=="A":
            for j in range(i+1):
                y[j] = y[j] + lattice_const
        if conf[i]=="B":
            for j in range(i+1):
                x[j] = x[j] - lattice_const         
        if conf[i]=="C":
            for j in range(i+1):
                y[j] = y[j] - lattice_const 
        if conf[i]=="D":
            for j in range(i+1):
                x[j] = x[j] + lattice_const 
        if conf[i]=="E":
            for j in range(i+1):
                z[j] = z[j] - lattice_const
        if conf[i]=="F":
            for j in range(i+1):
                z[j] = z[j] + lattice_const

    for i in range((len(conf)//2),len(conf)):
        if conf[i]=="A":
            for j in range(i,len(conf)):
                y[j+1] = y[j+1] - lattice_const  
        if conf[i]=="B":
            for j in range(i,len(conf)):
                x[j+1] = x[j+1] + lattice_const
        if conf[i]=="C":
            for j in range(i,len(conf)):
                y[j+1] = y[j+1] + lattice_const
        if conf[i]=="D":
            for j in range(i,len(conf)):
                x[j+1] = x[j+1] - lattice_const  
        if conf[i]=="E":
            for j in range(i,len(conf)):
                z[j+1] = z[j+1] + lattice_const
        if conf[i]=="F":
            for j in range(i,len(conf)):
                z[j+1] = z[j+1] - lattice_const
    

    for i in range(len(seq)):
        ofile.write("ATOM%7d%4s%5s%6d    %8.3f%8.3f%8.3f\n" % 
                    (i+1,"CA",name[seq[i]],i+1,x[i],y[i],z[i]))
      
    ofile.write("END\n")