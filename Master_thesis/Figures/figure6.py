# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 13:28:58 2022

@author: Yujia Ao
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import h5py
import pandas as pd
import math
import scipy.stats as sci


def cm2inch(value): 
    return value/2.54
def pearson(x,y):
    r,p=sci.pearsonr(x,y)

def set_plotparam(ylabels,i,ifre):
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    if ifre == 2:
        plt.xticks(np.arange(0.14,0.28,0.04),fontname='Times New Roman',fontsize=tick_size)   
    if ifre == 3:
        plt.xticks(np.arange(0.2,0.38,0.04),fontname='Times New Roman',fontsize=tick_size)   
    if ifre == 4:
        plt.xticks(np.arange(0.26,0.46,0.04),fontname='Times New Roman',fontsize=tick_size)  
    if i == 2 or i == 3:
        plt.yticks(np.arange(0.8,1.01,0.05),fontname='Times New Roman',fontsize=tick_size)  
    plt.xlabel('全局效率',fontname='SimSun',fontsize=label_size)
    plt.ylabel(ylabels[i],fontname='SimSun',fontsize=label_size)

main_dir = ['E:\HCP\Results_task\slow-1',
        'E:\HCP\Results_task\slow-2',
        'E:\HCP\Results_task\slow-3',
        'E:\HCP\Results_task\slow-4',
        'E:\HCP\Results_task\slow-5'];

ylabels = ['反应时','反应时','正确率','正确率']


for ifre in [2,3,4]:  
  
    d = main_dir[ifre]
    data_dir = os.path.join(d,'RT_ACC_Eg.mat')
    data_file = scio.loadmat(data_dir)
    RT0 = np.squeeze(data_file['RT0'])
    RT2 = np.squeeze(data_file['RT2'])
    ACC0 = np.squeeze(data_file['ACC0'])
    ACC2 = np.squeeze(data_file['ACC2'])
    Eg =np.squeeze( data_file['Eg'])
    
    Data = pd.DataFrame(data={'RT0': RT0, 'RT2': RT2,'ACC0': ACC0,'ACC2': ACC2,
                              'Eg': Eg})
    
    # 画图参数
    title_size = 8
    label_size = 7
    tick_size = 6
    scatter_size = 2
    line_size = 1
    
    
    # fig = plt.figure(figsize=(cm2inch(4.5), cm2inch(3.5)))
    grid = sns.jointplot(x='Eg',y='RT0',data=Data,ratio=4,
                scatter_kws={'s':scatter_size},line_kws={"color": "tab:red",'linewidth' : line_size},kind='reg')

    grid.annotate(sci.pearsonr,template='r={val:.2f}, p={p:.3f}',loc='best',
                  frameon=False,prop={'family' : 'Times New Roman','size':6})
    # r,p = pearson(Data.Eg,Data.RT0)

    # grid.annotate(pearson,template='r:{val:.2f},p:{p:.3f}')
    grid.fig.set_figwidth(cm2inch(4.5))
    grid.fig.set_figheight(cm2inch(3.5))
    plt.xlabel('全局效率',fontname='SimSun',fontsize=label_size,labelpad=1)
    plt.ylabel('')
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    plt.xticks(fontname='Times New Roman',fontsize=tick_size)
    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Fig6')

    plt.savefig('RT0_Slow-'+str(ifre+1)+'.tif',dpi=300)
    plt.show()
    
    
    grid2 = sns.jointplot(x='Eg',y='RT2',data=Data,ratio=4,
                scatter_kws={'s':scatter_size},line_kws={"color": "tab:red",'linewidth' : line_size},kind='reg')
    grid2.annotate(sci.pearsonr,template='r={val:.2f}, p={p:.3f}',loc='best',
                  frameon=False,prop={'family' : 'Times New Roman','size':6})
    grid2.fig.set_figwidth(cm2inch(4.5))
    grid2.fig.set_figheight(cm2inch(3.5))

    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    plt.xlabel('全局效率',fontname='SimSun',fontsize=label_size,labelpad=1)
    plt.ylabel('')
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    plt.xticks(fontname='Times New Roman',fontsize=tick_size)
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Fig6')

    plt.savefig('RT2_Slow-'+str(ifre+1)+'.tif',dpi=300)
    plt.show()
    
 
    grid3 = sns.jointplot(x='Eg',y='ACC0',data=Data,ratio=4,
                scatter_kws={'s':scatter_size},line_kws={"color": "tab:red",'linewidth' : line_size},kind='reg')
    grid3.annotate(sci.pearsonr,template='r={val:.2f}, p={p:.3f}',loc='best',
                  frameon=False,prop={'family' : 'Times New Roman','size':6})
    grid3.fig.set_figwidth(cm2inch(4.5))
    grid3.fig.set_figheight(cm2inch(3.5))

    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    plt.xlabel('全局效率',fontname='SimSun',fontsize=label_size,labelpad=1)
    plt.ylabel('')
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    plt.xticks(fontname='Times New Roman',fontsize=tick_size)
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Fig6')

    plt.savefig('ACC0_Slow-'+str(ifre+1)+'.tif',dpi=300)
    plt.show()
    
 
    grid4 = sns.jointplot(x='Eg',y='ACC2',data=Data,ratio=4,
                scatter_kws={'s':scatter_size},line_kws={"color": "tab:red",'linewidth' : line_size},kind='reg')
    grid4.annotate(sci.pearsonr,template='r={val:.2f}, p={p:.3f}',loc='best',
                  frameon=False,prop={'family' : 'Times New Roman','size':6})
    grid4.fig.set_figwidth(cm2inch(4.5))
    grid4.fig.set_figheight(cm2inch(3.5))

    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    plt.xlabel('全局效率',fontname='SimSun',fontsize=label_size,labelpad=1)
    plt.ylabel('')
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    plt.xticks(fontname='Times New Roman',fontsize=tick_size)
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Fig6')
    plt.savefig('ACC2_Slow-'+str(ifre+1)+'.tif',dpi=300)
    plt.show()
    
    
    
    
    
