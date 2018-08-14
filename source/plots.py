# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 14:00:29 2018

@author: Amy
"""
# =============================================================================
# Begins *** plots.py ***
# 
# Different plot functions to be used in mainplot.py. 
#   The 3 plot types can diplay a histogram of the native state porbability,
#   prbability of the number of native contact given 0 native contact in the 
#   initial reference frame, and the contact pair probability
#
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

#Dependencies
from source.io import *
from source.calc import *
from source.probTraj import*
from source.plots import *
from random import randint


def plotNumNative (protein_name, native, T, conf_list, num_sample\
                   , prob_traj_length, r = 100):
    """
    Plots the number of native contacts after n steps, given the condition that
    0 contacts appear in the initial chosen step. 
    
    protein_name is the name of the protein expressed as a string. See proteins
    folder for a list of currently available lattice proteins
    native is a list of integer values reresenting the number of natiive 
    contacts present at each MC step in a simulation trajectory
    T is the temperature of the simulation expressed as an int
    conf_list is the conformation trajectory expressed as a list of strings
    num_sample is an int for sample size chosen from the simualtion trajectory
    prob_traj_length is the length of the probability traj as an int
    r is the running average window size expressed as an int
    """
    legend = []
    from source.calc import get_native_list as get_native
    native_list = get_native(protein_name)
    
    fig, ax1 = plt.subplots()
    
    for i in range(len(native_list) + 1):
        c=np.random.rand(3,)
        print('Calculating ' + 'P(' + str(i) + '|0)...')
        p_traj = NativeProbTraj (conf_list,num_sample,prob_traj_length, native\
              , num_native = i, protein_name = protein_name)
        x = np.array(range(len(p_traj)))
        y = np.array(p_traj)
        
        yMA = running_avg(y, r)
        xMA = x[len(x)-len(yMA):]

        ax1.scatter(xMA,yMA,c=c, s=.01, marker='.')
            
        legend_bar = mpatches.Patch(color= c, label= 'P(' + str(i) + '|0)')
        legend.append(legend_bar)
        
    ax1.set_xlabel('MC Step')
    ax1.set_ylabel('Probability')
    plt.title(' Native Contact ' \
              + 'Probability of ' + protein_name + ' in ' + str(T) + 'K'\
              + ' RA=' + str(r) )
    plt.legend (handles = legend)
    plt.show()
    
def plotNative (protein_name, native, T, iteration,  r, environ = ''):
    """
    Functions plots the native state probability over a trajectory. 
    
    protein_name is the name of the protein expressed as a string. See proteins
    native is a list of integer values reresenting the number of natiive 
    contacts present at each MC step in a simulation trajectory
    T is the temperature of the simulation expressed as an int
    iteration is the length of the folding trajectory in MC steps
    r is the running average window size expressed as an int
    environ is a tring either "Static" or "Oscillating" to descript the temp.
    """
    legend = []
    from source.calc import get_native_list as get_native
    native_list = get_native(protein_name)
      
    fig, ax1 = plt.subplots()
    print('Calculating Native State Probability...')
    for i in range(len(native)):
                
        if native[i] == len(native_list):
            
            native[i] = 1
        else:
            native[i] = 0
            
    
    c=np.random.rand(3,)
    
    x = np.array(range(len(native)))
    y = np.array(native)
    
    yMA = running_avg(y, r)
    xMA = x[len(x)-len(yMA):]
    ax1.scatter(xMA,yMA,c= 'blue', s=.01, marker='.')        
    
    ax1.set_xlabel ('MC Step')
    ax1.set_ylabel ('Native State Moving Average')
    ax1.set_title ('Native State Prob Over '+ str(len(native))  + ' Steps at '\
                   + environ +' ' + str(T) +'K')
    ax1.set_xticks(np.arange(0,len(native),iteration))
    ax1.grid(True)
    plt.legend (handles = [mpatches.Patch(color= 'blue', label= 'RA = ' \
                                          + str(r))])
    plt.show()
    
    
    from scipy.stats import norm
    import matplotlib.mlab as mlab
    sns.distplot (yMA, color = 'blue', label = str(T) +'K')#, kde = True)
    plt.legend (handles = legend)
    
    (mu,sigma) = norm.fit(yMA)
#    n, bins, patches = plt.hist(yMA, 60,normed=1, facecolor= 'blue', alpha=0.5)
#    y = mlab.normpdf( bins, mu, sigma)
#    plt.plot(bins, y, 'r--', linewidth=2)
    
    plt.xlabel ('Probability')
    plt.ylabel ('Frequency')
    plt.title ('Native State Distribution'+ str(len(native))  + ' Steps at ' + environ +' ' + str(T) +'K')
    plt.legend (handles = [mpatches.Patch(color= 'blue', label= 'RA = ' + str(r))])
    
    plt.show()
    
    print('Mean = ', mu , ' Standard Dev =' ,sigma)
    
    
def plotContactPair (protein_name,a, b, T, p_traj, period, amplitude, r = [100000]):
    """
    Functions plots the contact pair probabaility over a trajectory with legend
    displaying the running average windows. 
    
    protein_name is the name of the protein expressed as a string. See proteins
    a and b are integers describing the Python indices of the contact pair
    T is the temperature of the simulation expressed as an int
    p_traj is the length of the probability traj as an int
    period is an int describing the period of the temperature oscillation
    amplitude is an int describing the amplitude of the temperature oscillation
    r is the running average window size expressed as a list of integers
    
    """
    
    legend = []
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    x = np.array(range(len(p_traj)))
    y = np.array(p_traj)
    x2 = np.linspace(0, len(p_traj), period-1)
    y2 = [square_func(i, T, amplitude, p = period) for i in x2]
#                plt.scatter(x,y, c= colors[z], marker = '.',  s = .01)
    ax2.plot (x2, y2, c='red') #plots square function                
    
    plt.title(str(a) + '-' +str(b) + ' Contact ' \
              + 'Probability of ' + protein_name + ' in ' + str(T) + 'K')
    ax1.set_xlabel('MC Step')
    ax1.set_ylabel('Probability', color = 'b')
    ax2.set_ylabel('Temperature (K)', color='r')
    
    ax2.set_ylim(T- amplitude - 10, T + amplitude + 10)
    
    for i in r:
        c=np.random.rand(3,)
        yMA = running_avg(y, i)
        xMA = x[len(x)-len(yMA):]
        legend_bar = mpatches.Patch(color= c, label= str(i) + ' RA')
        legend.append(legend_bar)   
        ax1.scatter(xMA,yMA,c=c, s=.01, marker='.')
            
    plt.legend (handles = legend)
    plt.show()