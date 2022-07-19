# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:48:35 2021

@author: felixS
"""

from scipy.stats import linregress
import math
from statistics import mean
from time import time
import base_py.basic_func as basic_func


####### BHT CORRECTION METHODS ######
####
####
####
####    
#### Horner Plot
def Horner_improved(BHT_list, t_list, tc):
    
    if len(BHT_list) != len(t_list):
        print('ERROR: Input Listen nicht gleich lang!')
    
    tau_list = []
    for each_t in t_list:
        tau = (tc+each_t)/each_t
        tau_list.append(tau)
    b, a, r, p, std = linregress(tau_list,BHT_list)
    SFT = b*1+a
    return SFT
####
####
####
####
####
#### Brennand Method
def Brennand(BHT_list, t_list, tc):
    
    if len(BHT_list) != len(t_list):
        print('ERROR: Input Listen nicht gleich lang!')
    
    Bt_list = []
    for each_t in t_list:  
        Bt = 1/(each_t+0.785*tc)
        Bt_list.append(Bt)
    b, a, r, p, std = linregress(Bt_list,BHT_list)
    SFT = b*0+a
    return SFT	
####
####
####
####    
#### Lachenbruch&Brewer Plot
def Lachenbruch(BHT_list, t_list):

    if len(BHT_list) != len(t_list):
        print('ERROR: Input Listen nicht gleich lang!')
        
    Lt_list = []
    for each_t in t_list:
        Lt = 1/each_t
        Lt_list.append(Lt)
    b, a, r, p, std = linregress(Lt_list,BHT_list)
    SFT = b*0+a
    return SFT
####
####
####
####
#### Zekiri 1-BHT
def SingleBHT(BHT, t, a, Tm, k):
    ft = math.exp(-(a**2/(4*k*t)))
    if Tm < BHT:
        SFT = (BHT+Tm*(ft-1))/ft
    else:
        SFT = (BHT+(BHT-0.1)*(ft-1))/ft
    return SFT
####
####
####
####
####
#### forward modelling 2 or more BHT
def forward_modelling(BHT_list, t_list, a, k):

    if len(BHT_list) != len(t_list):
        print('ERROR: Input Listen nicht gleich lang!')

    ft_list = []
    for each_t in t_list:
        ft = math.exp(-(a**2/(4*k*each_t)))
        ft_list.append(ft)
    
    SFT_list = []  
    
    #generate list of indices one entry smaller than BHT_list     
    for i in list(range(len(BHT_list)-1)):
        SFT_single = BHT_list[i] + (BHT_list[i+1]-BHT_list[i])*(1-ft_list[i])/(ft_list[i+1]-ft_list[i])
        SFT_list.append(SFT_single)
    
    #find average of list
    SFT = mean(SFT_list)
    
    return SFT
