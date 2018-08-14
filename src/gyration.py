#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Begins gyration.py

Gyration function computes the radius of gyration of a given conformation.
"""
from src.AAdict import *

def gyration(conformation, seq):

    import numpy as np

    # Lattice constant in \AA.
    lattice_const = 1.4

    # Protein sequence is from N terminal to C terminal.
    residue_mass = []
    for i in range(len(conformation)+1):
        residue_mass.append(wght[seq[i]] - 18.016)
    residue_mass[0] = residue_mass[0] + 1.008
    residue_mass[-1] = residue_mass[-1] + 17.008
    residue_mass = np.asarray(residue_mass)

    """
    Step 1. Construct coordinates from the conformation sequence.
    """
    x = []
    y = []
    z = []
    for i in range(len(conformation)+1):
        x.append(0)
        y.append(0)
        z.append(0)
    for i in range(len(conformation)):
        if conformation[i]=="A":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] - 1
        if conformation[i]=="B":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] + 1
        if conformation[i]=="C":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] + 1
        if conformation[i]=="D":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] - 1
        if conformation[i]=="E":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] + 1
        if conformation[i]=="F":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] - 1

    """
    Step 2. Compute the r_mean.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    mean_x = (residue_mass * x).sum() / residue_mass.sum()
    mean_y = (residue_mass * y).sum() / residue_mass.sum()
    mean_z = (residue_mass * z).sum() / residue_mass.sum()

    """
    Step 3. Compute the R_g.
    """

    # R_g_x = ((residue_mass * (x - mean_x)**2).sum()/residue_mass.sum())**0.5
    # R_g_y = ((residue_mass * (y - mean_y)**2).sum()/residue_mass.sum())**0.5
    # R_g_z = ((residue_mass * (z - mean_z)**2).sum()/residue_mass.sum())**0.5
    R_g = (((residue_mass * (x - mean_x)**2).sum() +\
            (residue_mass * (y - mean_y)**2).sum() +\
            (residue_mass * (z - mean_z)**2).sum() )\
            /residue_mass.sum())**0.5

    # return round(R_g_x * lattice_const, 2),\
    #        round(R_g_y * lattice_const, 2),\
    #        round(R_g_z * lattice_const, 2),\
    #        round(R_g * lattice_const, 2)
    return round(R_g * lattice_const, 2)
