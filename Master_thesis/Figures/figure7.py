# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 16:59:12 2022

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

def cm2inch(value): 
    return value/2.54

main_dir = ['E:\HCP\Results_PL\slow-1',
            'E:\HCP\Results_PL\slow-3',
            'E:\HCP\Results_PL\slow-4',
            'E:\HCP\Results_PL\slow-5'];
minv = np.array([-20,-31,-31,-24])
maxv = -minv

phases = ['FWE_T_peak', 'FWE_T_trough', 'FWE_T_ascend', 'FWE_T_descend']

for i in [0,1,2,3]:
    for j in phases:
        
        d = main_dir[i]
        
        T_dir = os.path.join(d,'FWE_T.mat')
        T_file = scio.loadmat(T_dir)
        T = T_file[j]
        
        """热图"""
        # 画图参数
        title_size = 8
        label_size = 7
        tick_size = 6
    
        x_ticks = []

        y_ticks = []
        
        # 画图
        fig = plt.figure(figsize=(cm2inch(3.2), cm2inch(2.6)))
        h = plt.axes([0.05,0.05,0.85,0.87])
        
        sns.heatmap(T,cmap='RdBu_r',
                    cbar_kws={'label':'T'},vmin=minv[i],vmax=maxv[i])
        h.spines['right'].set_visible(True)
        h.spines['left'].set_visible(True)
        h.spines['bottom'].set_visible(True)
        h.spines['top'].set_visible(True)
        
        
        
        # 坐标轴
        plt.xticks(x_ticks,rotation=360,fontname='Times New Roman',fontsize=tick_size)
        plt.yticks(y_ticks,fontname='Times New Roman',fontsize=tick_size)
        
        
        width = 0.5
        plt.axhline(69,lw=width,color='black')
        plt.axhline(125,lw=width,color='black')
        plt.axhline(163,lw=width,color='black')
        plt.axhline(175,lw=width,color='black')
        plt.axhline(189,lw=width,color='black')
        plt.axhline(211,lw=width,color='black')
        
        plt.axvline(69,lw=width,color='black')
        plt.axvline(125,lw=width,color='black')
        plt.axvline(163,lw=width,color='black')
        plt.axvline(175,lw=width,color='black')
        plt.axvline(189,lw=width,color='black')
        plt.axvline(211,lw=width,color='black')
        
        # 颜色条
        cb = h.collections[0].colorbar
        
        for l in cb.ax.yaxis.get_ticklabels():
            l.set_family('Times New Roman')
        # cb.set_ticks([-30,-20,-10,0,10,20,30])
        cb.ax.tick_params(labelsize=tick_size,length=2,pad=1)
        cb.set_label('T',fontdict={'family':'Times New Roman','size':label_size},labelpad=-2)
        
       
        
        os.chdir('E:\HCP\Results_PL\FiguresV0.2\Supplement_Fig1')
        # plt.savefig(img_name[i])
        plt.savefig(d[22:]+j+'.tif',dpi=300)
        plt.show()
    