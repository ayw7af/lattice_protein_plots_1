# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:11:31 2018

@author: Amy
"""
# =============================================================================
# Begins probTraj.py
#
# Calculates the probability trajectory of given protein.
#
# =============================================================================
import numpy as np
import random
from source.calc import isContact as isContact

def ContactProbTraj (i,j,conf_traj, num_sample, prob_traj_length\
                     , protein_name = 'mer15' ):
    """
    Calculates the probablity trajectory of one contact formation by flagging 
    whether contact exists in each conformation of conf_traj and selects from a
    list of valid indices
    
    i and j are Python indices of the contact pair chosen. 
    conf_traj is a list of lattice protein conformations represented by strings 
    num_sample is the integer sample size
    prob_traj_length is the length of the probability trajectory to be calc
    
    Returns list p_traj of length prob_traj_length.
    p_traj is the probability traj of contact pair formation 
    """
    p_traj = []
    #probability list of contact formation to be returned
    
    has_contact = []
    #list of 0s and 1s of the conformation trajectory 
    #1 if contact exists in reference frame, 0 otherwise
    
    valid_sample_ind = []
    #list of indices of element 0 from the has_contact trajectory
    #valid initial reference frame indices 
    
    for e in range (len(conf_traj)):
        #appends boolean flags to has_contact
        exists = isContact(i,j,conf_traj[e])
        has_contact.append (exists)
        
        #if contact exists, index of conformation in the traj is appended 
        if exists == 0 and (e < len(conf_traj) - prob_traj_length):
            valid_sample_ind.append(e)
    
    for i in range (prob_traj_length):
        if i == int(prob_traj_length / 4):
            print("25.0 % done.")
        elif i == int(prob_traj_length / 2):
            print("50.0 % done.")
        elif i == int(0.75 * prob_traj_length):
            print("75.0 % done.")

        occurs = 0
        
        for p in range(num_sample):
            
            # random conf selected from list of valid indices
            # n is the index of nth MC step
            n = valid_sample_ind [random.randint (0,len(valid_sample_ind)-1)]
            
            
            #records num times contact occurs at chosen MC step (either 0 or 1)
            occurs += has_contact [n + i]
            
        p_traj.append(occurs/num_sample)
            
    return p_traj

def NativeProbTraj (conf_list,num_sample,prob_traj_length, native\
              , num_native = 0, protein_name = 'mer15'): 
    """
    Calculates the probablity trajectory of the native state in a folding traj
    by flagging whether contact exists in each conformation of conf_traj and selects from a
    list of valid indices
 
    conf_traj is a list of lattice protein conformations represented by strings 
    num_sample is the integer sample size
    prob_traj_length is the length of the probability trajectory to be calc
    native is a list of integers indicating the number of native contacts 
    at each step of  conf_list
    
    Returns list p_traj of length prob_traj_length.
    p_traj is the probability traj of native conformation 
    """
    from source.calc import get_native_list as get_native_list
    
    p_traj = []
    #list of 0s and 1s of the conformation trajectory 
    has_contact = []
    
    #pool of indices to choose from 
    valid_sample_ind = []
    
    
    for e in native:
        if e == num_native:
            has_contact.append(1)
        else: 
            has_contact.append(0)
         
        if e == 0 and (e < len(conf_list) - prob_traj_length):
            valid_sample_ind.append(e)
       
    for i in range (prob_traj_length):
        if i == int(prob_traj_length / 4):
            print("25.0 % done.")
        elif i == int(prob_traj_length / 2):
            print("50.0 % done.")
        elif i == int(0.75 * prob_traj_length):
            print("75.0 % done.")
        
        occurs = 0
        
        for p in range(num_sample):
            
            # random conformation selected from list of valid indices
            # n is the index of nth MC step
            n = valid_sample_ind [np.random.randint (0,len(valid_sample_ind)-1)]
            
            #records num times contact occurs at chosen MC step (either 0 or 1)
            occurs += has_contact [n + i]
            
        p_traj.append(occurs/num_sample)
            
    return p_traj