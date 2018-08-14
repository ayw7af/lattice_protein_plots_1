#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Begins pfold.py.
# 
# Module to calculate the folding probability, p_fold, of a given conformation.
#
# In the future OOP, this will become a member function of the ``conf'' class.
# =============================================================================

from changeit import *
from energy import *
import sys, os, math, random
import matplotlib.pyplot as plt

sys.path.append(os.path.realpath("../proteins/"))
from merA import *
# from chignolin import *

def q(conf):
    """
    Computes the fraction of native contacts, q, of a given conformation.
        It returns a float number between 0 and 1, which is the 
        fraction of native_contacts formed.
    """
    temp, contact, native_contact = energy(conf, native_list, seq)
    fraction_of_native_contacts = float(native_contact) / len(native_list)
    assert (fraction_of_native_contacts >= 0) and \
           (fraction_of_native_contacts <= 1)
    return fraction_of_native_contacts
        

def pfold(conf, T, qlo = 0.2, qhi = 0.9, loop = 500):
    """
    Computes the folding probability, p_fold, of a given conformation.
        It returns a float number between 0 and 1, which is the 
        folding probability of that given conformation.
        
    By default, the limit of fraction of native contacts for unfolded state
        is 0.2, and that for folded state is 0.9. The number of loops is 500.
    """ 
    folded = 0
    assert len(conf) == len(seq) - 1
        
    e = 0
    contact = 0
    native_contact = 0
    metro = 0
    
    for i in range(loop):
        fraction = 0.5
        while (fraction > qlo) and (fraction < qhi):
            fraction = q(conf)
            
            newconf = changeit(conf)
            (newe, newcontact, newnative_contact) = \
                energy(newconf, native_list, seq)
            if newe<e:
                conf = newconf
                e = newe  
            else:
                # Metropolis criteria.
                metro = math.exp((e-newe)/T)
                if random.random() < metro:
                    conf = newconf
                    e = newe
        if (fraction >= qhi):
            folded = folded + 1
            
    return float(folded) / loop


def pfold_plot(conf_list, stats = 1):
    """
    Computes the pfold trajectory of a given conformation list.
        It returns a list of pfold trajectory, and a visualized plot.
        
    By default, the sampling interval is 1.
    """
    prob_list = []
    steps = len(conf_list) - 1
    
    # A percentage marker.
    if iteration > 100:
        percentage = int(iteration/100)
    else:
        percentage = 1
        
    for i in np.arange(0, steps, stats):       
        # A percentage output to screen
        if i % percentage == 0:
            print (i/percentage, "% done.")
        prob_list.append(pfold(conf_list[i], 1.0, loop = 100))
        
        
    # =========================================================================
    # Figure 3 plots the folding probability trajectory, p_fold.
    # =========================================================================
    fig3 = plt.figure()
    plt.scatter(np.arange(0.0, steps, stats), 
                prob_list, c='black', s=1, marker='.')
    #plt.scatter(np.arange(5000 * stats, steps, stats), 
    #            running_mean(fold_prob, 501), c='blue', s=0.0007, marker='.')
    plt.title("Folding Probability, $p_\mathrm{fold}$")
    plt.ylabel("$p_\mathrm{fold}$")
    plt.xlabel("Monte Carlo time (steps)")
    plt.ylim(-0.1, 1.1)