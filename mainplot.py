# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Fri Aug  3 15:04:14 2018

@author: Amy
"""
# =============================================================================
# Begins *** mainplot.py ***
# 
# Main plot function for the program. This function will call the appropriate
#   functions from source.plots dependng on the plotting argument
#
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import seaborn as sns

#Dependencies
from source.io import *
from source.calc import *
from source.probTraj import*
from source.plots import *
from random import randint

   
def plot (T, plotting = 'native', r= [10000], iteration = 1000000\
          , num_trial = 1, period = 10000, num_sample = 1000\
          , prob_traj_length = 100000, a = 1, b = 8\
          , amplitude = 100, protein_name = 'chignolin', data= 'read'):
    """
    Function displays plot(s) based on arguments given. 
    
    T is a list of temperatures in Kelvin.
    plotting is a string, either 'native' ,'contact pair', or 'num native'.
    'native' plots the native state probability trajectory as well as a histogram
    'contact pair' plots the contain pair prob traj
    'num native' plots the number of native contact probability
    r is the running average window 
    num_trial is the number of MC runs. Set to 1, by default.
    period is the period of the temp oscillation, set to 10,000 by default
    num_sample is the sample size 
    prob_traj_length is the length of the prob traj returned
    a is the integer value of the first contact in contact pair
    b is the integer value of the second contact in contact pair 
    amplitude is the temperature amplitude in Kelvin
    protein_name is the name of the protein expressed as a string
    data can be 'read' or 'write' to read the probability traj
    from file or write to file
    

    If type set to 'native',function plots the native conformation probability.
    If set to 'contact pair', function plots the contact pair probability traj.
    If set to 'num native', function plots the prob of number of native, else
    error message. 
    
    Period is 10,000 by default
    num_sample is 1,000 by default.
    protein_name is 'chignolin' by default. 
    data is set to 'read' by default
    
    prints out the mean and standard deviation according to a normal distrib
    returns "Done"
    
    Note: if data == 'read' but file does not exist, function will write to file
    first before plotting
    RA means "running average" on graph labels. 
    """
    from src.main import main
    
    #used for labeling graphs
    if amplitude == 0: 
        environ = 'Static'
    else: 
        environ = 'Oscillating'     
    
    #Creates a new graph for every Temp and trial or run        
    for z in range (len(T)):
        for n in range(num_trial):
            
            #generic file name for conformation trajectory
            conf_filename = protein_name + '_' + str(n + 1) + ' ' +str(T[z]) + 'K'\
    + ', Steps = ' + str(iteration) + ', Period = ' \
    + str(period) +  ', Amplitude = ' + str(amplitude)   
            
            if plotting == 'num native' or 'native': 
                
                # if number of native contacts need to be considered, create a 
                #file name that records the num of natives trajectory
                native_filename = 'NATIVES_' + str(n+1) + ' ' + conf_filename\
            +'.txt'
                
                conf_filename = conf_filename + '.txt'
                
                #if data needs to be read, and file exits, read file, else
                #run simulation and write new files
                if data == 'read' and os.path.exists(native_filename) \
                and os.path.exists(conf_filename) :
                    
                    native = read_from_file(name = native_filename)
                    conf_list = read_from_file(name = conf_filename)
            
                else: 
                    print('Running simulation... ')
                    acceptance, etotal, conf_list, native, all_contact \
                    = main(T[z], amplitude, period, iteration, protein_name)
                    write_to_file(conf_list, name = conf_filename)
                    write_to_file (native, name = native_filename) 
                
               
                if plotting == 'num native':
                     #plots number of native contact prob for each temperature, trial
                     # and running average
                    for i in r:
                        
                        #calls function found in source.plots
                        plotNumNative (protein_name, native, T[z], conf_list\
                                       , num_sample, prob_traj_length, i)              
                else:
                     #plots number of native contact prob for each temperature, trial
                     # and running average
                    for i in r:
                        #calls function found in source.plots
                        plotNative (protein_name, native, T[z]\
                                       , iteration, i, environ = environ)
                   
                
            elif plotting == 'contact pair':
                
                #additional info needed for cntact_pair calc
                conf_filename =  conf_filename + ', Samples = ' + str(num_sample)\
                    + ', Prob Traj length = ' + str(prob_traj_length) + '.txt'
                    
                #if data needs to be read, and file exits, read file, else
                #run simulation and write new files
                if data == 'read' and os.path.exists(conf_filename):
                    conf_list = read_from_file(name = conf_filename)
                else: 
                    print('Running simulation... ')
                    acceptance, etotal, conf_list, native, all_contact \
                    = main(T[z], amplitude, period, iteration, protein_name)
                    write_to_file(conf_list, name = conf_filename)
             
                from source.probTraj import ContactProbTraj
                
                #calulates the contact pair probability
                p_traj = ContactProbTraj(a,b,conf_list, num_sample, prob_traj_length\
                     , protein_name = protein_name)
                
                print('Plotting: ' + conf_filename)
                #function to display the contact pair probability plots
                plotContactPair (protein_name,a, b, T[z], p_traj,period\
                                 ,amplitude, r = r)
            
            #If none of the 3 possible types of plots entered for arg plotting
            else: 
                print ('Please set plotting argument to \'native\', \
                       \'num native\', or \'contact pair\'')
                break
            
    return "Done"

"""
Enter simulation parameters below.
"""
T = [450] 
plotting = 'native' #choice of 'contact pair', 'native', and 'num native'
iteration = 100000
r = [10000]
num_trial = 1
period = 1000
num_sample = 1000
prob_traj_length = 10000 #Enter arbitrary number if plotting is set to 'native'
i = 1 #Enter arbitrary num if arg. plotting is set to 'native' or 'num native'
j = 8 #Enter arbitrary num if arg. plotting is set to 'native' or 'num native'
amplitude = 100
protein_name = 'chignolin'
data = 'read'

"""
Enter simulatin parameters above.
"""
print(plot (T, plotting = plotting, r= r, iteration = iteration\
      , num_trial = num_trial, period = period\
      , num_sample =num_sample, prob_traj_length = prob_traj_length\
      , a = i, b = j, amplitude = amplitude, protein_name = protein_name\
      , data= data))