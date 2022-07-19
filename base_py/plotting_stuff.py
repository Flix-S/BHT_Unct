#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 16:46:44 2019

@author: hombre

define commonly used plot stuff from plotly
with same layout

"""


import plotly.graph_objs as go
import plotly
import plotly.io as pio
import os
import numpy as np
import pandas as pd
import ntpath
import math
import matplotlib.cm as cm
import base_py.basic_func as basic_func
import matplotlib.pyplot as plt
import seaborn as sns

#calculate timesteps and times once:
t_list = [0.0,1.0]
dt=1
t=1
for i in range(125):
    dt=dt*1.1
    t+=dt
    t_list.append(t)
t_list.append(1.8e6)

def gen_xaxis_dict(title,typ='log',rang=[-14,-9], dtick=1,showticklabels=True,showgrid=True):
    xaxis =dict(
                title=title,
                titlefont=dict(
                        size=15
                        ),
                type=typ,
                showgrid=showgrid,
                mirror=True,
                showline=True,
                zeroline=False,
                gridcolor='#000000',
                ticks='outside',
                showticklabels=showticklabels,
                tickcolor='#000000',
                autorange=False,
                range = rang, #log
                dtick = dtick,
                tickformat="g",
                tickfont=dict(
                        size=13
                        ),
                )
    return xaxis

def gen_yaxis_dict(title,typ='log',rang=[-17,-12], dtick=1,tickformat="1.2f",ticks='outside',showticklabels=True,showgrid=True):
    yaxis =dict(
                title=title,
                titlefont=dict(
                        size=15
                        ),
                type=typ,
                showgrid=showgrid,
                mirror=True,
                showline=True,
                zeroline=False,
                gridcolor='#000000',
                
                ticks=ticks,
                showticklabels=showticklabels,
                tickcolor='#000000',
                autorange=False,
                range = rang, #log
                dtick = dtick,
                tickformat=tickformat,
                tickfont=dict(
                        size=13
                        ),
                )
    return yaxis

def saveplot_routine(fig,indv_plotname,plotdir,newfoldername,auto_open=True,pdf=True,png=True,html=True):
    #generating plot dir for semilog plots:
    os.chdir(plotdir)
    if os.path.isdir('./'+newfoldername) == False:
        os.mkdir('./'+newfoldername)
        print('directory '+newfoldername+' created')
    plotpath = os.path.join(plotdir,newfoldername)
    
    filename_html = indv_plotname + '.html'
    filename_png = indv_plotname + '.png'
    filename_pdf = indv_plotname + '.pdf'
    
    filenamepath = os.path.join(plotpath,filename_html)
    
    if html == True:
        plotly.offline.plot(fig, filename=filenamepath, auto_open=auto_open)
    if png == True:
        pio.write_image(fig, os.path.join(plotpath,filename_png))
    if pdf == True:
        pio.write_image(fig, os.path.join(plotpath,filename_pdf))


def layout_f_plotly(xaxis,yaxis,width=850,height=850,annotations=[],fontsize=13):
    x = 0.5
    y = 0.8
    layout=dict(font=dict(family='Tahoma', 
                          size=fontsize, 
                          color='#000000'),
                width=width,
                height=height,
                showlegend=True,
                plot_bgcolor='rgb(255,255,255)',
                #paper_bgcolor='rgb(155,155,155)',
                legend=dict(x=x,
                            y=y,
                            tracegroupgap=20,
                            bgcolor='rgba(255,255,255,0.8)',
                            bordercolor='rgba(0,0,0,1.0)',
                            borderwidth=1),
                dragmode='pan',
                autosize=False)
    #append axes and annotations to dict:
    layout['annotations']=annotations
    layout['xaxis']=xaxis
    layout['yaxis']=yaxis
    
    return layout



def make_scatter_traces_withconfidence(traces,x,y,y_conf,name,rgb,style='solid'):
    p_trace=go.Scatter(showlegend=True,
                        x=x,
                        y=y,
                        error_y=dict(type='data',
                                     array=np.array(y_conf)/2,
                                     color=rgb,
                                     visible=True),
                        name=name,
                        hoverinfo = 'y+name',
                        mode = 'lines+markers',
                        line=dict(color=rgb,
                                  width=2,
                                  dash = style),
                        marker=dict(size=1,
                                    opacity=0,
                                    color=rgb,
                                    ),
                        )
    traces.append(p_trace)
    return traces

def SFT_plot(well,
             quality_BHT,
             quality_ts,
             Y,
             base,
             p10,
             p50,
             p90,
             case,
             TVD):
    
    f, (ax_box, ax_hist) = plt.subplots(2, 
       sharex=True,
       gridspec_kw={"height_ratios": (0.2, 1)},
       figsize=(5,4),
       dpi=150)
    
    flierprops = dict(markerfacecolor='1',
                      markersize=0.5,
                      linestyle='none',
                      alpha = 0.3)
    
    for style in ['ticks']:
        sns.set_style(style)  
        sns.boxplot(Y, 
                    ax=ax_box,
                    width=0.8,
                    linewidth=1,
                    flierprops=flierprops,
                    palette="Set3")
    #    ax_box.axvline(BHT_Soll,
    #                   color='r',
    #                   linewidth = 1.5,
    #                   linestyle="--")
        ax = sns.distplot(Y,
                          ax=ax_hist,
                          color="gray")
        x = ax.lines[0].get_xdata()
        y = ax.lines[0].get_ydata()
        mainSFT = np.round(x[np.argmax(y)],1)
        plt.axvline(x[np.argmax(y)],
                      color='grey')
        plt.axvline(base, color="b",
                    linestyle=":",
                    alpha = 0.8,
                    linewidth = 1.5)
        plt.axvline(p10,
                    color="k",
                    linestyle="--",
                    alpha = 0.8,
                    linewidth = 1.5)
        plt.text(p10,1.03,'p10',transform=ax.get_xaxis_transform())
        plt.axvline(p50,
                    color="k",
                    linestyle="--",
                    alpha = 0.8,
                    linewidth = 1.5)
        plt.text(p50,1.03,'p50',transform=ax.get_xaxis_transform())
        plt.axvline(p90, color="k",
                    linestyle="--",
                    alpha = 0.8,
                    linewidth = 1.5)
        plt.text(p90,1.03,'p90',transform=ax.get_xaxis_transform())
    #    plt.axvline(BHT_Soll,
    #                color="r",
    #                linestyle="--")
        plt.xlabel('Temperature [°C]',
                   size =15)
        plt.ylabel('Density',
                   size =15)
        ax.annotate('no. samples:',
                    xy=(0.01, 0.975), xycoords='axes fraction',
                    size = 5)
        ax.annotate(str(len(Y)),
                    xy=(0.03, 0.94), xycoords='axes fraction',
                    size = 5)        
        ax.annotate('p50 value = '+str(p50)+' °C',
                    xy=(0.98, 0.6), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate('base value: '+str(base)+' °C',
                    xy=(0.98, 0.5), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    color='blue',
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
    #    ax.annotate('SFT: '+str(BHT_Soll)+' °C',
    #                xy=(0.98, 0.5), xycoords='axes fraction',
    #                xytext=(0, 20), textcoords='offset pixels',
    #                size = 11,
    #                color='red',
    #                horizontalalignment='right',
    #                verticalalignment='bottom',
    #                alpha = 0.8,
    #                bbox=dict(facecolor='white',
    #                          alpha=0.7))
        ax.annotate('modal value: '+str(mainSFT)+' °C',
                    xy=(0.98, 0.4), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    color='grey',
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate('p10 value = '+str(p10)+' °C',
                    xy=(0.98, 0.8), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate('p90 value = '+str(p90)+' °C',
                    xy=(0.98, 0.7), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate(well + ' - ' + str(case),
                    xy=(0.98, 0.1), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate('depth = ' + str(TVD) + ' [mTVD]',
                    xy=(0.98, 0), xycoords='axes fraction',
                    xytext=(0, 20), textcoords='offset pixels',
                    size = 11,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    alpha = 0.8,
                    bbox=dict(facecolor='white',
                              alpha=0.7))
        ax.annotate("",
                    xy=(p10,0.01),
                    xycoords='data',
                    xytext=(p50,0.01),
                    textcoords='data',
                    arrowprops=dict(arrowstyle='<->',color="0"))
        dp1 = np.round(p50-p10,1)
        dp2 = np.round(p90-p50,1)
        ax.annotate('  ' + str(dp1) + '°C',
                    xy=(p10,0.01),
                    xycoords='data',
                    xytext=(p10,0.006),
                    textcoords='data',
                    size=8)
        ax.annotate("",
                    xy=(p50,0.01),
                    xycoords='data',
                    xytext=(p90,0.01),
                    textcoords='data',
                    arrowprops=dict(arrowstyle='<->',color="0"))
        ax.annotate('  ' + str(dp2) + '°C',
                    xy=(p50,0.01),
                    xycoords='data',
                    xytext=(p50,0.006),
                    textcoords='data',
                    size=8)                      
        ax.patch.set_facecolor('white')
        ax.patch.set_alpha(0.5)       
    plt.savefig('SFT '+case+' '+well+' '+quality_ts+' quality data set.png',bbox_inches='tight')
    



def plot_sobol(df_S1,
               df_S1_conf,
               df_ST,
               df_ST_conf,
               indv_plotname,
               workingdir,
               plottitle='Global Sensitivity Analysis',
               plot_ST=True,
               plot_S1=False,
               width=850,
               height=850):
    
    #read x values = sample size which is in headers:
    x = df_S1.columns.tolist()
    del x[0:2]
    x = [float(i) for i in x]
    
    x_sorted = sorted(x)
    
    sortingindex = []
    for each in x_sorted:
        sortingindex.append(x.index(each))
        
    #get variable names:
    l_names = df_S1['names'].tolist()
    
    traces = []
    cma = cm.rainbow
    
    for index,eachname in enumerate(l_names):
        y_S1 = df_S1.loc[df_S1['names']==eachname,:].drop(df_S1.columns[[0,1]], axis=1).values.flatten().tolist()
        y_S1 = [y_S1[i] for i in sortingindex]
        
        y_S1_conf = df_S1_conf.loc[df_S1_conf['names']==eachname,:].drop(df_S1_conf.columns[[0,1]], axis=1).values.flatten().tolist()
        y_S1_conf = [y_S1_conf[i] for i in sortingindex]
        
        y_ST = df_ST.loc[df_ST['names']==eachname,:].drop(df_ST.columns[[0,1]], axis=1).values.flatten().tolist()
        y_ST = [y_ST[i] for i in sortingindex]
        
        y_ST_conf = df_ST_conf.loc[df_ST_conf['names']==eachname,:].drop(df_ST_conf.columns[[0,1]], axis=1).values.flatten().tolist()
        y_ST_conf = [y_ST_conf[i] for i in sortingindex]
        
        cn=basic_func.normalize(index,0,len(l_names)-1 )
        c = cma(cn,1)
        rgb = c[:3]
        rgb = "rgb(%s, %s, %s)" % (int(rgb[0]*255),int(rgb[1]*255),int(rgb[2]*255))
        
        if plot_S1 == True:
            traces = make_scatter_traces_withconfidence(traces,x_sorted,y_S1,y_S1_conf,eachname+'_S1',rgb,style='dash')
        if plot_ST == True: 
            traces = make_scatter_traces_withconfidence(traces,x_sorted,y_ST,y_ST_conf,eachname+'_ST',rgb,style='solid')
    
    
    xaxis =dict(
                title='no. of samples',
                titlefont=dict(family='Tahoma',
                               size=17),
                type='linear',
                showgrid=True,
                mirror=True,
                showline=True,
                zeroline=False,
                gridcolor='#000000',
                ticks='outside',
                showticklabels=True,
                tickcolor='#000000',
                autorange=True,
                tickformat="g",
                tickfont=dict(family='Tahoma',
                              size=15),
                )
    
    
    yaxis =dict(
                title='sensitivity index',
                titlefont=dict(family='Tahoma',
                               size=17),
                type='linear',
                showgrid=True,
                mirror=True,
                showline=True,
                zeroline=False,
                gridcolor='#000000',
                ticks='outside',
                showticklabels=True,
                tickcolor='#000000',
                autorange=True,
                tickformat="1.2f",
                tickfont=dict(family='Tahoma',
                              size=15),
                )
    
    annotations=[]
    an = dict(x=0.40,
              y=0.95,
              showarrow=False,
              text=plottitle,
              bgcolor='rgb(255,255,255)',
              opacity=0.7,
              bordercolor='#000000',
              borderwidth=2,
              align='center',
              xref='paper',
              yref='paper',
              font =dict(family='Tahoma', 
                         size=17, 
                         color='#000000'))
    annotations.append(an)
    
    
    layout = layout_f_plotly(xaxis,yaxis,width=width,height=height,annotations=annotations)
    
    fig1 = dict(data=traces, layout=layout)
    
    saveplot_routine(fig1,indv_plotname,workingdir,'sobol_plots')          



def SalibBarplot(Si_df):
    
    fig, ax = plt.subplots(1)
    
    indices = Si_df[['S1','ST']]
    err = Si_df[['S1_conf','ST_conf']]
    
    indices.plot.bar(yerr=err.values.T,ax=ax)
    fig.set_size_inches(8,4)
    indv_plotname = "BarPlot"
    plt.show()
    saveplot_routine(fig,indv_plotname,workingdir,'sobol_plots')          


