# -*- coding: utf-8 -*-
"""
Spyder Editor

@author: felixS
"""
import math
from numba import jit
import numpy as np
#import time



@jit(nopython=True)              # Set "nopython" mode for best performance, equivalent to @njit



def iterative_Goetzl_3(BHT1, BHT2, BHT3, t1, t2, t3, a, Tm):
    ###  erstelle leere Listen
    BHT_list = [BHT1,BHT2,BHT3]
    t_list = [t1,t2,t3] 
#    BHT_list = [90,95,102]
#    t_list = [18000,20000,24000]
#    a = 0.008
#    Tm = 70
    Tf = np.zeros(shape=(1000))
    Tf_opt = np.zeros(shape=(150))
    sum_square = np.zeros(shape=(1000))
    sum_square_opt = np.zeros(shape=(150))
    l = len(BHT_list)
    k_start = 1e-7      #initiale Temp Leitfähigkeit
    kappa = k_start 
    step = 0.2
    residuum = np.zeros(shape=(3,1000))
    residuen_opt = np.zeros(shape=(3,150))
    signum = np.zeros(shape=(1000))
        
    ###   große k-Schleife, Erhöhung kappa
    for k in range(0,150):
        temp_delta = Tm*0.5    #initiale Temperaturstörung
        j = 0
        Tf[j]=Tm+temp_delta
        sum_square[j] = 0
        
        # 1. Inversionsdurchgang
        for i in list(range(l)):
            exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
            BHT_calc = Tf[j] + temp_delta*exponent_fl
            residuum_single = BHT_calc - BHT_list[i]
            residuum[i,j] = residuum_single
            sum_square[j] = sum_square[j]+(residuum[i,j]**2)
#            time.sleep(0.2)
            
        sum_square_min = sum_square[j]
        j_opt = j
        signum[j] = 1
        temp_delta = temp_delta + step *math.sqrt(sum_square[j])
#        time.sleep(0.2)
    
        
        #   j-Schleife - Inverse Optimierung der Formationstemperatur
        for j in range(1,1000):
            Tf[j] = (Tm + temp_delta)
            sum_square[j] = 0
            
            for i in list(range(l)):
                exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
                BHT_calc = Tf[j] + temp_delta*exponent_fl
                residuum[i,j] = BHT_calc - BHT_list[i]
                sum_square[j] = sum_square[j] + residuum[i,j]**2
#                time.sleep(0.2)
            
            if sum_square[j] < sum_square_min:
                j_opt = j
                sum_square_min = sum_square[j]
#                time.sleep(0.2)
            
            signum[j] = (sum_square[j-1] - sum_square[j])/(abs(sum_square[j-1] - sum_square[j])) * signum[j-1]
            
            temp_delta = temp_delta + step * signum[j]
#            time.sleep(0.2)
            
        sum_square_opt[k] = sum_square[j_opt]
        Tf_opt[k] = Tf[j_opt]
        
        for i in list(range(l)):
            residuen_opt[i,k] = residuum[i,j_opt]
#            time.sleep(0.2)
            ### Erhöhung kappa
        kappa = kappa + 5E-9
    
    sum_square_result = sum_square_opt[0]
    k_opt = 0
    
    for k in range(1,150):
        if sum_square_opt[k] < sum_square_result:
            k_opt = k
            sum_square_result = sum_square_opt[k]
    
    sum_square_result = sum_square_opt[k_opt]
    SFT = Tf_opt[k_opt]
    return SFT

def iterative_Goetzl_4(BHT1, BHT2, BHT3, BHT4, t1, t2, t3, t4, a, Tm):
    BHT_list = [BHT1,BHT2,BHT3,BHT4]
    t_list = [t1,t2,t3,t4] 
    
    ###  erstelle leere Listen    
    Tf = np.zeros(shape=(1000))
    Tf_opt = np.zeros(shape=(150))
    sum_square = np.zeros(shape=(1000))
    sum_square_opt = np.zeros(shape=(150))
    l = len(BHT_list)
    k_start = 1e-7      #initiale Temp Leitfähigkeit
    kappa = k_start 
    step = 0.2
    residuum = np.zeros(shape=(4,1000))
    residuen_opt = np.zeros(shape=(4,150))
    signum = np.zeros(shape=(1000))

    ###   große k-Schleife, Erhöhung kappa
    for k in range(0,150):
        temp_delta = Tm*0.5    #initiale Temperaturstörung
        j = 0
        Tf[j]=Tm+temp_delta
        sum_square[j] = 0
        
        # 1. Inversionsdurchgang
        for i in list(range(l)):
            exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
            BHT_calc = Tf[j] + temp_delta*exponent_fl
            residuum_single = BHT_calc - BHT_list[i]
            residuum[i,j] = residuum_single
            sum_square[j] = sum_square[j]+(residuum[i,j]**2)
            
        sum_square_min = sum_square[j]
        j_opt = j
        signum[j] = 1
        temp_delta = temp_delta + step *math.sqrt(sum_square[j])
        
        #   j-Schleife - Inverse Optimierung der Formationstemperatur
        for j in range(1,1000):
            Tf[j] = (Tm + temp_delta)
            sum_square[j] = 0
            
            for i in list(range(l)):
                exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
                BHT_calc = Tf[j] + temp_delta*exponent_fl
                residuum[i,j] = BHT_calc - BHT_list[i]
                sum_square[j] = sum_square[j] + residuum[i,j]**2
            
            if sum_square[j] < sum_square_min:
                j_opt = j
                sum_square_min = sum_square[j]
            
            signum[j] = (sum_square[j-1] - sum_square[j])/(abs(sum_square[j-1] - sum_square[j])) * signum[j-1]
            temp_delta = temp_delta + step * signum[j]
            
        sum_square_opt[k] = sum_square[j_opt]
        Tf_opt[k] = Tf[j_opt]
        
        for i in list(range(l)):
            residuen_opt[i,k] = residuum[i,j_opt]
            
            ### Erhöhung kappa
        kappa = kappa + 5E-9
    
    sum_square_result = sum_square_opt[0]
    k_opt = 0
    
    for k in range(1,150):
        if sum_square_opt[k] < sum_square_result:
            k_opt = k
            sum_square_result = sum_square_opt[k]
    
    sum_square_result = sum_square_opt[k_opt]
    SFT = Tf_opt[k_opt]
    return SFT

def iterative_Goetzl_5(BHT1, BHT2, BHT3, BHT4, BHT5, t1, t2, t3, t4, t5, a, Tm):
    ###  erstelle leere Listen
    BHT_list = [BHT1,BHT2,BHT3,BHT4,BHT5]
    t_list = [t1,t2,t3,t4,t5] 
#    BHT_list = [90,95,102]
#    t_list = [18000,20000,24000]
#    a = 0.008
#    Tm = 70
    Tf = np.zeros(shape=(1000))
    Tf_opt = np.zeros(shape=(150))
    sum_square = np.zeros(shape=(1000))
    sum_square_opt = np.zeros(shape=(150))
    l = len(BHT_list)
    k_start = 1e-7      #initiale Temp Leitfähigkeit
    kappa = k_start 
    step = 0.2
    residuum = np.zeros(shape=(5,1000))
    residuen_opt = np.zeros(shape=(5,150))
    signum = np.zeros(shape=(1000))

    ###   große k-Schleife, Erhöhung kappa
    for k in range(0,150):
        temp_delta = Tm*0.5    #initiale Temperaturstörung
        j = 0
        Tf[j]=Tm+temp_delta
        sum_square[j] = 0
        
        # 1. Inversionsdurchgang
        for i in list(range(l)):
            exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
            BHT_calc = Tf[j] + temp_delta*exponent_fl
            residuum_single = BHT_calc - BHT_list[i]
            residuum[i,j] = residuum_single
            sum_square[j] = sum_square[j]+(residuum[i,j]**2)
            
        sum_square_min = sum_square[j]
        j_opt = j
        signum[j] = 1
        temp_delta = temp_delta + step *math.sqrt(sum_square[j])
    
        
        #   j-Schleife - Inverse Optimierung der Formationstemperatur
        for j in range(1,1000):
            Tf[j] = (Tm + temp_delta)
            sum_square[j] = 0
            
            for i in list(range(l)):
                exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
                BHT_calc = Tf[j] + temp_delta*exponent_fl
                residuum[i,j] = BHT_calc - BHT_list[i]
                sum_square[j] = sum_square[j] + residuum[i,j]**2
            
            if sum_square[j] < sum_square_min:
                j_opt = j
                sum_square_min = sum_square[j]
            
            signum[j] = (sum_square[j-1] - sum_square[j])/(abs(sum_square[j-1] - sum_square[j])) * signum[j-1]
            
            temp_delta = temp_delta + step * signum[j]
            
        sum_square_opt[k] = sum_square[j_opt]
        Tf_opt[k] = Tf[j_opt]
        
        for i in list(range(l)):
            residuen_opt[i,k] = residuum[i,j_opt]
            
            ### Erhöhung kappa
        kappa = kappa + 5E-9
    
    sum_square_result = sum_square_opt[0]
    k_opt = 0
    
    for k in range(1,150):
        if sum_square_opt[k] < sum_square_result:
            k_opt = k
            sum_square_result = sum_square_opt[k]
    
    sum_square_result = sum_square_opt[k_opt]
    SFT = Tf_opt[k_opt]
    return SFT

def iterative_Goetzl_6(BHT1, BHT2, BHT3, BHT4, BHT5, BHT6, t1, t2, t3, t4, t5, t6, a, Tm):
    ###  erstelle leere Listen
    BHT_list = [BHT1,BHT2,BHT3,BHT4,BHT5,BHT6]
    t_list = [t1,t2,t3,t4,t5,t6] 
#    BHT_list = [90,95,102]
#    t_list = [18000,20000,24000]
#    a = 0.008
#    Tm = 70
    Tf = np.zeros(shape=(1000))
    Tf_opt = np.zeros(shape=(150))
    sum_square = np.zeros(shape=(1000))
    sum_square_opt = np.zeros(shape=(150))
    l = len(BHT_list)
    k_start = 1e-7      #initiale Temp Leitfähigkeit
    kappa = k_start 
    step = 0.2
    residuum = np.zeros(shape=(6,1000))
    residuen_opt = np.zeros(shape=(6,150))
    signum = np.zeros(shape=(1000))
    
    ###   große k-Schleife, Erhöhung kappa
    for k in range(0,150):
        temp_delta = Tm*0.5    #initiale Temperaturstörung
        j = 0
        Tf[j]=Tm+temp_delta
        sum_square[j] = 0
        
        # 1. Inversionsdurchgang
        for i in list(range(l)):
            exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
            BHT_calc = Tf[j] + temp_delta*exponent_fl
            residuum_single = BHT_calc - BHT_list[i]
            residuum[i,j] = residuum_single
            sum_square[j] = sum_square[j]+(residuum[i,j]**2)
            
        sum_square_min = sum_square[j]
        j_opt = j
        signum[j] = 1
        temp_delta = temp_delta + step *math.sqrt(sum_square[j])
    
        
        #   j-Schleife - Inverse Optimierung der Formationstemperatur
        for j in range(1,1000):
            Tf[j] = (Tm + temp_delta)
            sum_square[j] = 0
            
            for i in list(range(l)):
                exponent_fl = (math.exp(-(a**2/(4*kappa*t_list[i])))-1)
                BHT_calc = Tf[j] + temp_delta*exponent_fl
                residuum[i,j] = BHT_calc - BHT_list[i]
                sum_square[j] = sum_square[j] + residuum[i,j]**2
            
            if sum_square[j] < sum_square_min:
                j_opt = j
                sum_square_min = sum_square[j]
            
            signum[j] = (sum_square[j-1] - sum_square[j])/(abs(sum_square[j-1] - sum_square[j])) * signum[j-1]
            
            temp_delta = temp_delta + step * signum[j]
            
        sum_square_opt[k] = sum_square[j_opt]
        Tf_opt[k] = Tf[j_opt]
        
        for i in list(range(l)):
            residuen_opt[i,k] = residuum[i,j_opt]
            
            ### Erhöhung kappa
        kappa = kappa + 5E-9
    
    sum_square_result = sum_square_opt[0]
    k_opt = 0
    
    for k in range(1,150):
        if sum_square_opt[k] < sum_square_result:
            k_opt = k
            sum_square_result = sum_square_opt[k]
    
    sum_square_result = sum_square_opt[k_opt]
    SFT = Tf_opt[k_opt]
    return SFT

        