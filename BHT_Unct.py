#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 10:12:41 2019
@author: felixS
"""
import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
import os 
from time import time
import base_py.bht_functions as bht
import base_py.Iterativ_Goetzl as It
import base_py.basic_func as basic_func
import base_py.plotting_stuff as pls

# if possible, multiprocessing can be activated by uncomment lines 19, .... Multiprocessing is highly recommended when LM method is applied
# import multiprocessing as mp

# Functions are located in extra function python file: bht_functions.py and Iterative_Goetzl.py
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Insert name of well:
well = 'Name of Well'                   
method = 'BM'                               # specify correction method  '1BHTM' (automatically specified when n = 1)
                                            #                            'LM' (low confidence in shut-in time, n>2)
                                            #                            'FM','BM' (high confidence in shut-in time, n>1)
TVD = 1000                                # depth of BHT measurement [m TVD]
radius = 0.07620                            # borehole radius [m]
BHT_list = [70,72]                          # in-situ temperatures [degC]
t_list = [8,10.5]                           # shut-in times [h]
Tmud = 54                                   # estimate Tm as Outflow temperature, if such is available, otherwise: Tm = 0, initial range for Tm is [20,80]
quality_BHT = 'low'                         # confidence in quality of BHT and shut-in time 'low', 'med', 'high'
quality_ts = 'low'

sample_variation = [10,100,500,5000]           # specify the sample number
n = len(BHT_list)                           # number of BHT measurements in a row
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if n == 1:
    case = '1BHTM'
else:
    case = method
tc = 7.2 * TVD
kappa = 3E-7
if Tmud > BHT_list[0]:
    Tm = BHT_list[0]-0.5
else:
    Tm = Tmud
ts = [t_list[i]*3600 for i in range(len(BHT_list))]
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
###############################################################################
#Define Parameter Sets appropriate to quality of BHT data
###############################################################################
if quality_BHT == 'low':
    BHTZek = [BHT_list[0] * 0.97,BHT_list[0] * 1.03]
    BHT_lmin = [BHT_list[i]*0.97 for i in range(len(BHT_list))]
    BHT_lmax = [BHT_list[i]*1.03 for i in range(len(BHT_list))]
    for x in range(len(BHT_lmin)-1):
        if BHT_lmin[x+1]<=BHT_lmax[x]:
            f = (BHT_lmax[x]+BHT_lmin[x+1])/2
            BHT_lmin[x+1]=f+0.05
            BHT_lmax[x]=f-0.05
if quality_BHT == 'med':
    BHTZek = [BHT_list[0] * 0.98,BHT_list[0] * 1.02]
    BHT_lmin = [BHT_list[i]*0.98 for i in range(len(BHT_list))]
    BHT_lmax = [BHT_list[i]*1.02 for i in range(len(BHT_list))]
    for x in range(len(BHT_lmin)-1):
        if BHT_lmin[x+1]<=BHT_lmax[x]:
            f = (BHT_lmax[x]+BHT_lmin[x+1])/2
            BHT_lmin[x+1]=f+0.05
            BHT_lmax[x]=f-0.05
if quality_BHT == 'high':
    BHTZek = [BHT_list[0] * 0.99,BHT_list[0] * 1.01]
    BHT_lmin = [BHT_list[i]*0.99 for i in range(len(BHT_list))]
    BHT_lmax = [BHT_list[i]*1.01 for i in range(len(BHT_list))]
    for x in range(len(BHT_lmin)-1):
        if BHT_lmin[x+1]<=BHT_lmax[x]:
            f = (BHT_lmax[x]+BHT_lmin[x+1])/2
            BHT_lmin[x+1]=f+0.05
            BHT_lmax[x]=f-0.05   
###############################################################################
a_list = [radius, radius + 0.0254] 
k_list = [3E-7*0.5, 3E-7*1.5]
tc_min = tc-7200
if tc_min < 3600:
    tc_min = 3600
tc_list =[tc_min, tc+7200]
Tm_max = BHT_lmin[0]-0.5
if Tm_max > 80:
    Tm_max = 80
Tm_list = [20,Tm_max]
###############################################################################
if quality_ts == 'low':
    tZek = [t_list[0]*3600 - 3600,t_list[0]*3600 + 3600]
    t_lmin = [t_list[i]*3600-3600 for i in range(len(BHT_list))]
    t_lmax = [t_list[i]*3600+3600 for i in range(len(BHT_list))]
    for x in range(len(t_lmin)-1):
        if t_lmin[x+1]<=t_lmax[x]:
            f = (t_lmax[x]+t_lmin[x+1])/2
            t_lmin[x+1]=f+300
            t_lmax[x]=f-300 
if quality_ts == 'med':
    tZek = [t_list[0]*3600 - 1800,t_list[0]*3600 + 1800]
    t_lmin = [t_list[i]*3600-1800 for i in range(len(BHT_list))]
    t_lmax = [t_list[i]*3600+1800 for i in range(len(BHT_list))]
    for x in range(len(t_lmin)-1):
        if t_lmin[x+1]<=t_lmax[x]:
            f = (t_lmax[x]+t_lmin[x+1])/2
            t_lmin[x+1]=f+300
            t_lmax[x]=f-300       
if quality_ts == 'high':
    tZek = [t_list[0]*3600 - 900,t_list[0]*3600 + 900]
    t_lmin = [t_list[i]*3600-900 for i in range(len(BHT_list))]
    t_lmax = [t_list[i]*3600+900 for i in range(len(BHT_list))]
    for x in range(len(t_lmin)-1):
        if t_lmin[x+1]<=t_lmax[x]:
            f = (t_lmax[x]+t_lmin[x+1])/2
            t_lmin[x+1]=f+300
            t_lmax[x]=f-300        
###############################################################################
if (radius + 0.0254)*(radius + 0.0254)/(4*k_list[0]*t_lmin[0])>1:
    print(' ')
    print(' !!! stability criterion not fulfilled !!! ')
    print(' ')   
    
###############################################################################    
###############################################################################
###############################################################################
###############################################################################

# define the problem
if case == '1BHTM': 
    p_names = ['BHT','t_shut','radius','T_Mud','kappa']
    no_var = 5
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [a_list[1]-a_list[0],0.01],
                          [Tm_list[1]-Tm_list[0],0.5],
                          k_list],
               'dists':['unif','unif','triang','triang','unif']}
if case == 'BM' and len(BHT_list) == 2:                             #Brennand
    p_names = ['1st_BHT','t1','2nd_BHT','t2','tc']
    no_var = 5
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [tc_list[1]-tc_list[0],0.5]],
               'dists':['unif','unif','unif','unif','triang']}
if case == 'BM' and len(BHT_list) == 3:                             #Brennand
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','tc']
    no_var = 7
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],                         
                          [tc_list[1]-tc_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','triang']}                          
if case == 'BM' and len(BHT_list) == 4:                             #Brennand
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','tc']
    no_var = 9
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],    
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],                               
                          [tc_list[1]-tc_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','triang']}                 
if case == 'BM' and len(BHT_list) == 5:                             #Brennand
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','5th_BHT','t5',
               'tc']
    no_var = 11
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],    
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],     
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],                          
                          [tc_list[1]-tc_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang']}                            
if case == 'BM' and len(BHT_list) == 6:                             #Brennand
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','5th_BHT','t5',
               '6th_BHT','t6','tc']
    no_var = 13
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],  
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]], 
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],
                          [BHT_lmin[5],BHT_lmax[5]],
                          [t_lmin[5],t_lmax[5]],                            
                          [tc_list[1]-tc_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang']}     
if case == 'LM' and len(BHT_list) == 3:                            #Goetzl
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','radius','Tm']
    no_var = 8
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [a_list[1]-a_list[0],0.01],
                          [Tm_list[1]-Tm_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','triang','triang']}                               
if case == 'LM' and len(BHT_list) == 4:                            #Goetzl
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','radius','Tm']
    no_var = 10
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],                          
                          [a_list[1]-a_list[0],0.01],
                          [Tm_list[1]-Tm_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','triang','triang']}                            
if case == 'LM' and len(BHT_list) == 5:                            #Goetzl
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','5th_BHT','t5',
               'radius','Tm']
    no_var = 12
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],                          
                          [a_list[1]-a_list[0],0.01],
                          [Tm_list[1]-Tm_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang','triang']}                          
if case == 'LM' and len(BHT_list) == 6:                            #Goetzl
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','5th_BHT','t5',
               '6th_BHT','t6','radius','Tm']
    no_var = 14
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],
                          [BHT_lmin[5],BHT_lmax[5]],
                          [t_lmin[5],t_lmax[5]],                           
                          [a_list[1]-a_list[0],0.01],
                          [Tm_list[1]-Tm_list[0],0.5]],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang','triang']} 
if case == 'FM' and len(BHT_list) == 2:                             #Forward 
    p_names = ['1st_BHT','t1','2nd_BHT','t2','radius','kappa']
    no_var = 6
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],                    
                          [a_list[1]-a_list[0],0.01],
                          k_list],
               'dists':['unif','unif','unif','unif','triang','unif']}                          
if case == 'FM' and len(BHT_list) == 3:                             #Forward 
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','radius','kappa']
    no_var = 8
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],                          
                          [a_list[1]-a_list[0],0.01],
                          k_list],
               'dists':['unif','unif','unif','unif','unif','unif','triang','unif']}  
if case == 'FM' and len(BHT_list) == 4:                             #Forward 
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','radius','kappa']
    no_var = 10
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],                             
                          [a_list[1]-a_list[0],0.01],
                          k_list],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','triang','unif']}                                                     
if case == 'FM' and len(BHT_list) == 5:                             #Forward 
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
               't3','4th_BHT','t4','5th_BHT','t5',
               'radius','kappa']
    no_var = 12
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],                          
                          [a_list[1]-a_list[0],0.01],
                          k_list],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang','unif']}
if case == 'FM' and len(BHT_list) == 6:                             #Forward 
    p_names = ['1st_BHT','t1','2nd_BHT','t2','3rd_BHT',
              't3','4th_BHT','t4','5th_BHT','t5',
               '6th_BHT','t6','radius','kappa']
    no_var = 14
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[BHT_lmin[0],BHT_lmax[0]],
                          [t_lmin[0],t_lmax[0]],
                          [BHT_lmin[1],BHT_lmax[1]],
                          [t_lmin[1],t_lmax[1]],
                          [BHT_lmin[2],BHT_lmax[2]],
                          [t_lmin[2],t_lmax[2]],
                          [BHT_lmin[3],BHT_lmax[3]],
                          [t_lmin[3],t_lmax[3]],
                          [BHT_lmin[4],BHT_lmax[4]],
                          [t_lmin[4],t_lmax[4]],
                          [BHT_lmin[5],BHT_lmax[5]],
                          [t_lmin[5],t_lmax[5]],                          
                          [a_list[1]-a_list[0],0.01],
                          k_list],
               'dists':['unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','unif','triang','unif']}
# calculate the base value
if case == '1BHTM':
   base = bht.SingleBHT(BHT_list[0],ts[0],radius,Tm,kappa)
elif case == 'BM':
   base = bht.Brennand(BHT_list,ts,tc)
elif case == 'FM':
   base = bht.forward_modelling(BHT_list,ts,radius,kappa)
elif case == 'LM' and len(BHT_list) == 3:
   base = It.iterative_Goetzl_3(BHT_list[0],BHT_list[1],BHT_list[2],ts[0],
                                ts[1],ts[2],radius,Tm)
elif case == 'LM' and len(BHT_list) == 4:
   base = It.iterative_Goetzl_4(BHT_list[0],BHT_list[1],BHT_list[2],BHT_list[3],
                              ts[0],ts[1],ts[2],ts[3],radius,Tm)
elif case == 'LM' and len(BHT_list) == 5:
   base = It.iterative_Goetzl_5(BHT_list[0],BHT_list[1],BHT_list[2],BHT_list[3],
                              BHT_list[4],ts[0],ts[1],ts[2],
                              ts[3],ts[4],radius,Tm)
elif case == 'LM' and len(BHT_list) == 6:
   base = It.iterative_Goetzl_6(BHT_list[0],BHT_list[1],BHT_list[2],BHT_list[3],
                              BHT_list[4],BHT_list[5],ts[0],ts[1],
                              ts[2],ts[3],ts[4],ts[5],radius,Tm)     
###############################################################################
# do the sensitivity analysis for multiple samples
############################################################################### 

# specify the save path here:
savepath = os.path.abspath('')   
savepath_plots = os.path.abspath('plots')

# generate empty dataframes for Sobol indices sotring:
df_ST = pd.DataFrame()
df_S1 = pd.DataFrame()
df_S2 = pd.DataFrame()
df_ST['names']=df_S1['names']=p_names

df_ST_conf = pd.DataFrame()
df_S1_conf = pd.DataFrame()
df_S2_conf = pd.DataFrame()
df_ST_conf['names']=df_S1_conf['names']=p_names

# starting the actual Sensitivity Analysis here:
for eachsamplenumber in sample_variation:
    
    param_values = saltelli.sample(problem, eachsamplenumber)    #generates N*(2*D+2)samples D = number of inputs = 6; N = userinput
    no_samples = (no_var*2+2)*eachsamplenumber
    Y = np.zeros([param_values.shape[0]])
    #just some console output and convenience stuff:
    print
    print(str(no_samples)+' satelli samples generated!')
    print('Start calculating models for these samples:')
    print
    l = len(param_values)
    t1 = time()

# print_progress(0, l, prefix = 'Progress:', suffix = 'Complete', bar_length = 50)
    for i, X in enumerate(param_values):
        # Calculate your model with parameter combination X here: 2BHT       
        if case == '1BHTM':
            BHT = X[0]
            t = X[1]
            a = X[2]+a_list[0]
            Tm = X[3]+Tm_list[0]
            k = X[4]
            SFT = bht.SingleBHT(BHT, t, a, Tm, k)        
        elif len(BHT_list) == 2 and case == 'BM':
            BHT_l = [X[0],
                     X[2]]
            t_l = [X[1],
                   X[3]]
            tc = X[4]+tc_list[0]
            SFT = bht.Brennand(BHT_l, t_l, tc)
        elif len(BHT_list) == 3 and case == 'BM':
            BHT_l = [X[0],
                     X[2],
                     X[4]]
            t_l = [X[1],
                   X[3],
                   X[5]]
            tc = X[6]+tc_list[0]
            SFT = bht.Brennand(BHT_l, t_l, tc)
        elif len(BHT_list) == 4 and case == 'BM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7]]
            tc = X[8]+tc_list[0]
            SFT = bht.Brennand(BHT_l, t_l, tc)            
        elif len(BHT_list) == 5 and case == 'BM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6],
                     X[8]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7],
                   X[9]]
            tc = X[10]+tc_list[0]
            SFT = bht.Brennand(BHT_l, t_l, tc)
        elif len(BHT_list) == 6 and case == 'BM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6],
                     X[8],
                     X[10]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7],
                   X[9],
                   X[11]]
            tc = X[12]+tc_list[0]
            SFT = bht.Brennand(BHT_l, t_l, tc)                               
        elif len(BHT_list) == 3 and case == 'LM':
            BHT1 = X[0]
            BHT2 = X[2]
            BHT3 = X[4]
            t1 = X[1]
            t2 = X[3]
            t3 = X[5]
            a = X[6]+a_list[0]
            Tm = X[7]+Tm_list[0]
            SFT = It.iterative_Goetzl_3(BHT1,BHT2,BHT3,t1,t2,t3,a,Tm)
        elif len(BHT_list) == 4 and case == 'LM':
            BHT1 = X[0]
            BHT2 = X[2]
            BHT3 = X[4]
            BHT4 = X[6]            
            t1 = X[1]
            t2 = X[3]
            t3 = X[5]
            t4 = X[7]            
            a = X[8]+a_list[0]
            Tm = X[9]+Tm_list[0]
            SFT = It.iterative_Goetzl_4(BHT1,BHT2,BHT3,BHT4,t1,t2,t3,t4,a,Tm)
        elif len(BHT_list) == 5 and case == 'LM':
            BHT1 = X[0]
            BHT2 = X[2]
            BHT3 = X[4]
            BHT4 = X[6]
            BHT5 = X[8]            
            t1 = X[1]
            t2 = X[3]
            t3 = X[5]
            t4 = X[7]
            t5 = X[9]             
            a = X[10]+a_list[0]
            Tm = X[11]+Tm_list[0]
            SFT = It.iterative_Goetzl_5(BHT1,BHT2,BHT3,BHT4,BHT5,t1,t2,t3,t4,t5,a,Tm)
        elif len(BHT_list) == 6 and case == 'LM':
            BHT1 = X[0]
            BHT2 = X[2]
            BHT3 = X[4]
            BHT4 = X[6]
            BHT5 = X[8]
            BHT6 = X[10]             
            t1 = X[1]
            t2 = X[3]
            t3 = X[5]
            t4 = X[7]
            t5 = X[9]
            t6 = X[11]             
            a = X[12]+a_list[0]
            Tm = X[13]+Tm_list[0]
            SFT = It.iterative_Goetzl_6(BHT1,BHT2,BHT3,BHT4,BHT5,BHT6,t1,t2,t3,t4,t5,t6,a,Tm)            
        elif len(BHT_list) == 2 and case == 'FM':
            BHT_l = [X[0],
                     X[2]]
            t_l = [X[1],
                   X[3]]
            a = X[4]+a_list[0]
            k = X[5]
            SFT = bht.forward_modelling(BHT_l, t_l, a, k)
        elif len(BHT_list) == 3 and case == 'FM':
            BHT_l = [X[0],
                     X[2],
                     X[4]]
            t_l = [X[1],
                   X[3],
                   X[5]]
            a = X[6]+a_list[0]
            k = X[7]
            SFT = bht.forward_modelling(BHT_l, t_l, a, k)
        elif len(BHT_list) == 4 and case == 'FM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7]]
            a = X[8]+a_list[0]
            k = X[9]
            SFT = bht.forward_modelling(BHT_l, t_l, a, k)            
        elif len(BHT_list) == 5 and case == 'FM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6],
                     X[8]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7],
                   X[9]]
            a = X[10]+a_list[0]
            k = X[11]
            SFT = bht.forward_modelling(BHT_l, t_l, a, k)
        elif len(BHT_list) == 6 and case == 'FM':
            BHT_l = [X[0],
                     X[2],
                     X[4],
                     X[6],
                     X[8],
                     X[10]]
            t_l = [X[1],
                   X[3],
                   X[5],
                   X[7],
                   X[9],
                   X[11]]
            a = X[12]+a_list[0]
            k = X[13]
            SFT = bht.forward_modelling(BHT_l, t_l, a, k)     
#        elif len(BHT_list) == 2 and case == 'Horn':
#            BHT_l = [X[0],
#                     X[2]]
#            t_l = [X[1],
#                   X[3]]
#            tc = X[4]
#            SFT = bht.Horner_improved(BHT_l, t_l, tc) 
#        elif len(BHT_list) == 2 and case == 'Lach':
#            BHT_l = [X[0],
#                     X[2]]
#            t_l = [X[1],
#                   X[3]]
#            SFT = bht.Lachenbruch(BHT_l, t_l)
        #Save it to array:
        Y[i] = SFT
        
###############################################################################
    # Do the Sobol Analysis here:
###############################################################################             
    Si = sobol.analyze(problem, Y)
    
    #storing sobol indices
    df_S1[no_samples]=Si['S1']
    df_ST[no_samples]=Si['ST']
    df_secondorder = pd.DataFrame(Si['S2'])
    df_secondorder['sample_no']=no_samples
    df_S2 = pd.concat([df_S2,df_secondorder])
    
    #storing confidence intervals:
    df_S1_conf[no_samples]=Si['S1_conf']
    df_ST_conf[no_samples]=Si['ST_conf']
    df_secondorder_conf = pd.DataFrame(Si['S2_conf'])
    df_secondorder_conf['sample_no']=no_samples
    df_S2_conf = pd.concat([df_S2_conf,df_secondorder_conf])

    #save this information to csv
    #for every iteration of this for loop the saved files will be updated/overridden:
    if os.path.isdir(os.path.join(savepath,'SobolInformation')) == False:
            os.mkdir(os.path.join(savepath,'SobolInformation'))
    
    soboldir = os.path.join(savepath,'SobolInformation')

    #saving sobol indices to csv:
    df_S1.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_S1.csv'))
    df_S2.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_S2.csv'))
    df_ST.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_ST.csv'))
    
    #saving sobol confidence intervals to csv:
    df_S1_conf.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_S1_conf.csv'))
    df_S2_conf.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_S2_conf.csv'))
    df_ST_conf.to_csv(os.path.join(soboldir,well+' '+str(case)+'_df_conf.csv'))
    
    #just some console output and convenience stuff:
    t2 = time()
    t_s=basic_func.tidy(t2-t1,4)
    t_m=basic_func.tidy((t2-t1)/60,3)
    t_h=basic_func.tidy((t2-t1)/3600,2)
    print
    print(str(t_s)+'s/'+str(t_m)+'m/'+str(t_h)+'h elapsed to finish sample batch')
    print

############################################################
#Plot the Sensitivity Analysis results:
############################################################

#now plot the Sobol information and save it:
if case == '1BHTM':
    indv_plotname='1BHTM ('+str(len(BHT_list))+'no. BHT values)'
elif case == 'BM':
    indv_plotname='BM ('+str(len(BHT_list))+'no. BHT values)'
elif case == 'FM':
    indv_plotname='FM ('+str(len(BHT_list))+'no. BHT values)'
elif case == 'LM':
    indv_plotname='LM ('+str(len(BHT_list))+'no. BHT values)'

#this function will generate a new folder inside savepath
#and save a sobol plot in 3 different formats into it
#using indv_plotname as filename
#you can turn on/off total and first order Sobol indices    
pls.plot_sobol(df_S1,
               df_S1_conf,
               df_ST,
               df_ST_conf,
               indv_plotname,
               savepath_plots,
               plottitle=indv_plotname,
               plot_ST=True,
               plot_S1=False,
               width=1000,
               height=500)

#just some console output and convenience stuff:
print 
print('Everything done! Sensitivity Analysis finished :)')
print

p10=np.round(np.percentile(Y, 10),1)
p20=np.percentile(Y, 20)
p90=np.round(np.percentile(Y, 90),1)
p50=np.round(np.percentile(Y, 50),1)
p70=np.percentile(Y, 70)
f = max(Y)
mean = np.mean(Y)
median = np.median(Y)
dp = p90-p10
base=np.round(base,1)

print('p10 case: ', p10)
print('p90 case: ', p90)
print('p50 case: ', p50)
print('Base value: ', base)

############################################################
#Plot the Statistics:
############################################################
pls.SFT_plot(well,
             quality_BHT,
             quality_ts,
             Y,
             base,
             p10,
             p50,
             p90,
             case,
             TVD)

std = np.std(Y)
std


