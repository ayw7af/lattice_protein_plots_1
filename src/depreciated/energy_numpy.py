#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# Begins energy.py.
#
# Module to compute the energy of a given conformation.
# Returns one tuple of three values.
# A float value indicates the energy of the input conformation.
# Two integer values containing the number of contacts and native contacts.
#
# In the future OOP, this will become a member function of the ``conf'' class.
# =============================================================================
from src.AAdict import aa as aa, imat as imat
import numpy as np

def energy(conformation, seq, native_list):
    """
    Depreciated due to slowness.
    """

    N = len(conformation) + 1
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)

    for i in range(N - 1):
        if conformation[i]=="A": y[i+1:] = y[i+1:] - 1
        elif conformation[i]=="B": x[i+1:] = x[i+1:] + 1
        elif conformation[i]=="C": y[i+1:] = y[i+1:] + 1
        elif conformation[i]=="D": x[i+1:] = x[i+1:] - 1
        elif conformation[i]=="E": z[i+1:] = z[i+1:] + 1
        else: z[i+1:] = z[i+1:] - 1

    energy = 0
    total_contact = 0
    native_contact = 0

    for (i, j) in native_list:
        if abs(x[i]-x[j])+abs(y[i]-y[j])+abs(z[i]-z[j])==1:
            native_contact = native_contact + 1

    # Self-avoiding checks.
    for i in range(len(conformation)-1):
        for j in range(i+2, len(conformation)+1):
            if (abs(x[i]-x[j])**2 + abs(y[i]-y[j])**2 + abs(z[i]-z[j])**2)==0:
                return float('nan'), 0, 0

    for i in range(len(conformation)-2):
        for j in range(i+3, len(conformation)+1):
            if abs(x[i]-x[j])+abs(y[i]-y[j])+abs(z[i]-z[j])==1:
                energy = energy + imat[aa.index(seq[i])][aa.index(seq[j])]
                total_contact = total_contact + 1

    return energy, total_contact, native_contact