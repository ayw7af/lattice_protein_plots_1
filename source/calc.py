# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 08:57:40 2018

@author: Amy Wang
"""
# =============================================================================
# Begins calc.py.
#
# Calculates averages and number of native contacts.
#
# Functions here are used to create the probability trajectory (probTraj.py)
# =============================================================================

def get_native_list (protein_name):
    """
    protein_name is a string representing the lattice protein. 
    Returns list of tuples. Each tuple contains the native contact pair indices
    Currently there are 5 proteins: "mer4", "mer48A", "mer15", "chignolin"
    """
    if protein_name == "mer4":
        from proteins.mer4 import native_list as native_list
    elif protein_name == "mer48A":
        from proteins.mer48A import native_list as native_list
    elif protein_name == "mer27":
        from proteins.mer27 import native_list as native_list
    elif protein_name == "mer15":
        from proteins.mer15 import native_list as native_list
    elif protein_name == "chignolin":
        from proteins.chignolin import native_list as native_list
    else: raise AssertionError("Invalid protein name. Please check the input!")
    
    return native_list 

def cum_avg(mylist):
    """
    Calculates the cumulative average of an integer list (mylist).
    Returns float list of cumulative averages.
    Ex: cum_avg ([0,1,2,3,4,5]) returns [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
    """
    cumsum, cum_aves = [0], []
    
    for i, x in enumerate(mylist, 1):
        cumsum.append(cumsum[i-1] + x)
        cum_ave = (cumsum[i])/(i)
        cum_aves.append(cum_ave)
        
    return cum_aves


def running_avg (mylist, N):
    """
    Calculates the moving average of a list (x) with N integer window size. 
    Returns a numpy array of moving average values
    Ex: running_avg ([0,1,2,3,4,5], 2) returns [ 0.5  1.5  2.5  3.5  4.5]
    """
    import numpy as np
    
    cumsum = np.cumsum(np.insert(mylist, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def isContact (i, j, conformation):
    """
    i and j are  Python indices of contact residues in lattice protein.
    Conformation is latice protein conformation represented by a list of strings
    """       
    from src.energy import getCoord as getCoord
    
    x,y,z = getCoord(conformation) 
   
    if (abs(x[i]-x[j])**2 + abs(y[i]-y[j])**2 + abs(z[i]-z[j])**2) == 1:
        
        return 1
    else: 
        return 0

def num_native (conformation, protein_name = 'mer15'):
    """
    Conformation is latice protein conformation represented by a list of strings.
    protein_name is a string representing the lattice protein. 
    Currently there are 5 proteins: "mer4", "mer48A", "mer15", "chignolin"
    Returns the integer number of native contacts in the given conformation
    
    Default protein_name is set to 'mer15' , a lattice protein with length 15
    Conformation is represented by capital letters each representing one of the 
    six possible directions in a Cartesian coord system.
    """
    from src.energy import getCoord as getCoord
    
    x,y,z = getCoord(conformation) 
    native_list = get_native_list (protein_name)
    
    count = 0 
    for e in native_list:
        i = e[0]
        j = e[1]
        count += isContact(i,j, conformation)
    return count 


def square_func(i, T, amp, p = 10000):
    """
    Function called to create square function for temp
    i is the element in the x axis 
    T is the temperature to oscillate around
    amp is the amplitude
    p is the period
    p is set to 10,000 MC steps by default
    
    returns integer T + amp and T-amp for each period
    """
    if (i//p)%2 == 0:
        return T + amp
    else:
        return T - amp