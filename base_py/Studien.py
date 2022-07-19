# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 10:15:34 2021

@author: ga56qan
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 18 11:27:45 2021

@author: ga56qan
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 10:12:41 2019

@author: hombre
"""
import sys, os
sys.path.append('E:/Daten/Temperaturen/Studies')

import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
from time import time
import base_py.Iterativ_Goetzl as It
import base_py.bht_functions as bht
import base_py.basic_func as basic_func
import base_py.plotting_stuff as pls
from base_py.progressbar import print_progress
from SALib.plotting.bar import plot as barplot

############################################################
#Functions are located in extra function python file: bht_functions.py
############################################################

BHT_Soll = 84.57
############################################################
#Set up Sensitivity problem:
############################################################

# Define the model inputs
#1 BHT
#section = '_6.25inch'
analyse_typ = 'Brennand' #Brennand, Horner, Lachenbruch, Goetzl, Forward
n=2
#for n in [2,3]:
#           'bounds': [[58,62],
#                      [20,62],
#                      [0.079375,0.10795],
#                      [1.5E-7,1.64E-6],
#                      [19800,125000]]
           
if n == 2 and analyse_typ == 'Horner' or n == 2 and analyse_typ == 'Brennand': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'tc']
    no_var = 5
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[48.5,51.5],
                          [43800,45600],
                          [1.85,8.15],
                          [9720, 13320],
                          [7572.74, 21972.74]]
               }

if n == 2 and analyse_typ == 'Lachenbruch': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1']
    no_var = 4
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[48.5,51.5],
                          [43800,45600],
                          [1.85,8.15],
                          [9720, 13320]]
               }

if n == 3 and analyse_typ == 'Horner' or n == 3 and analyse_typ == 'Brennand' or n == 3 and analyse_typ == 'Lachenbruch':
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'delta_2ndBHT',
               'delta_t2',
               'tc']
    no_var = 7
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[59.17,62.83],
                          [18500,170000],
                          [0.1,36],
                          [1800, 125000],
                          [0.1,36],
                          [1800, 125000],
                          [3600, 144000]]
               }
    
if n == 2 and analyse_typ == 'Goetzl': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'radius',
               'Tm']
    no_var = 6
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[48.5,51.5],
                          [43800,45600],
                          [1.85,8.15],
                          [9720, 13320],
                          [0.071628,0.082042],
                          [40,60]]
               }
    
if n == 3 and analyse_typ == 'Goetzl': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'delta_2ndBHT',
               'delta_t2',
               'radius',
               'Tm']
    no_var = 8
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[59.17,62.83],
                          [18500,170000],
                          [0.1,36],
                          [1800, 125000],
                          [0.1,36],
                          [1800, 125000],
                          [0.079375,0.104775],
                          [24,60]]
               }
    
if n == 2 and analyse_typ == 'Forward': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'radius',
               'kappa']
    no_var = 6
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[48.5,51.5],
                          [43800,45600],
                          [1.85,8.15],
                          [9720, 13320],
                          [0.071628,0.082042],
                          [1.5E-7,6.8E-7]]
               }
    
if n == 3 and analyse_typ == 'Forward': 
    p_names = ['1st_BHT',
               't1',
               'delta_1stBHT',
               'delta_t1',
               'delta_2ndBHT',
               'delta_t2',
               'radius',
               'kappa']
    no_var = 8
    problem = {'num_vars': no_var,
               'names': p_names,
               'bounds': [[59.17,62.83],
                          [18500,170000],
                          [0.1,36],
                          [1800, 125000],
                          [0.1,36],
                          [1800, 125000],
                          [0.079375,0.104775],
                          [1.5E-7,6.8E-7]]
               }
#do the sensetivity analysis for multiple sample sizes to test convergence:
# graphische
#sample_variation = [100,1000,5000,10000,50000,100000]
#sample_variation = [100,1000,5000,10000]    
#Götzl 
#sample_variation = [100,5000,12000,20000]

sample_variation = [100,5000]
#analytische
#sample_variation = [100,5000,8000,12000,20000,75000]

#specify the save path here, a new folder will be generated inside:
savepath = os.getcwd()



#generate empty dataframes for Sobol indices sotring:
df_ST = pd.DataFrame()
df_S1 = pd.DataFrame()
df_S2 = pd.DataFrame()
df_ST['names']=df_S1['names']=p_names

df_ST_conf = pd.DataFrame()
df_S1_conf = pd.DataFrame()
df_S2_conf = pd.DataFrame()
df_ST_conf['names']=df_S1_conf['names']=p_names



#starting the actual Sensitivity Analysis here:
for eachsamplenumber in sample_variation:
    
    param_values = saltelli.sample(problem, eachsamplenumber)#generates N*(2*D+2)samples D = number of inputs = 6; N = userinput
    no_samples = (no_var*2+2)*eachsamplenumber
    Y = np.zeros([param_values.shape[0]])
    
    
    #just some console output and convenience stuff:
    print
    print(str(no_samples)+' satelli samples generated!')
    print('Start calculating models for these samples:')
    print
    l = len(param_values)
    t1 = time()
    print_progress(0, l, prefix = 'Progress:', suffix = 'Complete', bar_length = 50)
    
    #zum testen parameter input gekürzt:
    #param_values=param_values[0:1]
    
    for i, X in enumerate(param_values):
        
        #Calculate your model with parameter combination X here: 2BHT
        if n == 2 and analyse_typ == 'Horner':
            BHT_list = [X[0],X[0]+X[2]]
            t_list = [X[1],X[1]+X[3]]
            tc = X[4]
            SFT = bht.Horner_improved(BHT_list, t_list, tc)
        elif n == 3 and analyse_typ == 'Horner':
            BHT_list = [X[0],X[0]+X[2],X[0]+X[2]+X[4]]
            t_list = [X[1],X[1]+X[3],X[1]+X[3]+X[5]]
            tc = X[6]
            SFT = bht.Horner_improved(BHT_list, t_list, tc)
        elif n == 2 and analyse_typ == 'Brennand':
            BHT_list = [X[0],X[0]+X[2]]
            t_list = [X[1],X[1]+X[3]]
            tc = X[4]
            SFT = bht.Brennand(BHT_list, t_list, tc)
        elif n == 3 and analyse_typ == 'Brennand':
            BHT_list = [X[0],X[0]+X[2],X[0]+X[2]+X[4]]
            t_list = [X[1],X[1]+X[3],X[1]+X[3]+X[5]]
            tc = X[6]
            SFT = bht.Brennand(BHT_list, t_list, tc)
        elif n == 2 and analyse_typ == 'Lachenbruch':
            BHT_list = [X[0],X[0]+X[2]]
            t_list = [X[1],X[1]+X[3]]
            SFT = bht.Lachenbruch(BHT_list, t_list)
        elif n == 3 and analyse_typ == 'Lachenbruch':
            BHT_list = [X[0],X[0]+X[2],X[0]+X[2]+X[4]]
            t_list = [X[1],X[1]+X[3],X[1]+X[3]+X[5]]
            SFT = bht.Lachenbruch(BHT_list, t_list) 
            
        elif n == 2 and analyse_typ == 'Goetzl':
            BHT1 = X[0]
            BHT2 = X[0]+X[2]
            t1 = X[1]
            t2 = X[1]+X[3]
            a = X[4]
            Tm = X[5]
            SFT = It.iterative_Goetzl_2BHT(BHT1,BHT2,t1,t2,a,Tm)
        elif n == 3 and analyse_typ == 'Goetzl':
            BHT1 = X[0]
            BHT2 = X[0]+X[2]
            BHT3 = X[0]+X[2]+X[4]
            t1 = X[1]
            t2 = X[1]+X[3]
            t3 = X[1]+X[3]+X[5]
            a = X[6]
            Tm = X[7]
            SFT = It.iterative_Goetzl(BHT1,BHT2,BHT3,t1,t2,t3,a,Tm)             
        elif n == 2 and analyse_typ == 'Forward':
            BHT_list = [X[0],X[0]+X[2]]
            t_list = [X[1],X[1]+X[3]]
            a = X[4]
            k = X[5]
            SFT = bht.forward_modelling(BHT_list, t_list, a, k)
        elif n == 3 and analyse_typ == 'Forward':
            BHT_list = [X[0],X[0]+X[2],X[0]+X[2]+X[4]]
            t_list = [X[1],X[1]+X[3],X[1]+X[3]+X[5]]
            a = X[6]
            k = X[7]
            SFT = bht.forward_modelling(BHT_list, t_list, a, k)
            
        
        #Save it to array:
        Y[i] = SFT
        

        
        #just some console output and convenience stuff:
        #update progressbar:
        print_progress(i + 1, l, prefix = 'Progress:', suffix = 'Complete', bar_length = 50)
    
    
    #Do the Sobol Analysis here:
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
    df_S1.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_S1.csv'))
    df_S2.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_S2.csv'))
    df_ST.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_ST.csv'))
    
    
    #saving sobol confidence intervals to csv:
    df_S1_conf.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_S1_conf.csv'))
    df_S2_conf.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_S2_conf.csv'))
    df_ST_conf.to_csv(os.path.join(soboldir,analyse_typ+'_AschheimTh2_df_conf.csv'))
    
    
    
    
    
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
if analyse_typ == 'Horner':
    indv_plotname='Horner SA ('+str(n)+'no. BHT values)'
elif analyse_typ == 'Brennand':
    indv_plotname='Brennand SA ('+str(n)+'no. BHT values)'
elif analyse_typ == 'Lachenbruch':
    indv_plotname='Lachenbruch SA ('+str(n)+'no. BHT values)'
elif analyse_typ == 'Goetzl':
    indv_plotname='Iterative Method SA_TestParallel ('+str(n)+'no. BHT values)'
elif analyse_typ == 'Forward':
    indv_plotname='Forward Modelling SA_2BHT ('+str(n)+'no. BHT values)'


#this function will generate a new folder inside savepath
#and save a sobol plot in 3 different formats into it
#using indv_plotname as filename
#you can turn on/off total and first order Sobol indices    
pls.plot_sobol(df_S1,
               df_S1_conf,
               df_ST,
               df_ST_conf,
               indv_plotname,
               savepath,
               plottitle=indv_plotname,
               plot_ST=True,
               plot_S1=False,
               width=1000,
               height=500)



#just some console output and convenience stuff:
print 
print('Everything done! Sensitivity Analysis finished :)')
print
    
p10=np.percentile(Y, 10)
p90=np.percentile(Y, 90)
p50=np.percentile(Y, 50)
f = max(Y)
min_diff = BHT_Soll - min(Y)
max_diff = max(Y) - BHT_Soll
mean_diff = BHT_Soll - np.mean(Y)

print('negative Abweichung zu T-Log: ', min_diff)
print('positive Abweichung zu T-Log: ', max_diff)
print('mittlere Abweichung zu T-Log: ', mean_diff)

print('p10 case: ', p10)
print('p90 case: ', p90)
print('p50 case: ', p50)
print('Mean value: ', np.mean(Y))
print('Median: ', np.median(Y))
import seaborn as sns
from matplotlib import pyplot as plt
sns.set(rc={"figure.figsize": (5, 5)}); np.random.seed(0)
sns.set(font_scale = 1.2)
ax = sns.distplot(Y)
plt.axvline(np.mean(Y), color="k", linestyle="--")
plt.axvline(p10, color="k", linestyle="--")
plt.axvline(p90, color="k", linestyle="--")
plt.axvline(BHT_Soll, color="r", linestyle="--")
plt.xlabel('Temperature [°C]', size =15)
ax.annotate(analyse_typ,
            xy=(1.05, 0.8), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            size = 18,
            horizontalalignment='right',
            verticalalignment='bottom')
std = np.std(Y)
std
#-----------------------------------------------------------------------------
