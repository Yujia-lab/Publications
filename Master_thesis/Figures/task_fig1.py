# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 16:52:05 2021

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
"""图二，不同相位的空间相似性"""
main_dir = ['E:\HCP\Results_task\slow-1',
        'E:\HCP\Results_task\slow-2',
        'E:\HCP\Results_task\slow-3',
        'E:\HCP\Results_task\slow-4',
        'E:\HCP\Results_task\slow-5'];
title_name = ['Slow-1','Slow-2','Slow-3','Slow-4','Slow-5'];
img_name = ['CVC_Slow1.pdf','CVC_Slow2.pdf',
            'CVC_Slow3.pdf','CVC_Slow4.pdf',
            'CVC_Slow5.pdf']
img_name2 = ['CVC_Slow1.tif','CVC_Slow2.tif',
            'CVC_Slow3.tif','CVC_Slow4.tif',
            'CVC_Slow5.tif']

for i in range(0,5):
    d = main_dir[i]
    # 读取数据
    intensity_phase_dir = os.path.join(d,'intensity_phase.mat')
    intensity_phase_file = scio.loadmat(intensity_phase_dir)
    GScorr_CVC = intensity_phase_file['CVC_phase']
    mean_CVC_r_dir = os.path.join(d,'mean_cl_r.mat')
    mean_cl_r_file = scio.loadmat(mean_CVC_r_dir)
    mean_cl_r = mean_cl_r_file['mean_cl_r']
    """热图"""
    # 画图参数
    title_size = 8
    label_size = 7
    tick_size = 6

    
    x = np.arange(0,361,90)
    x_ticks = ['-π','-π/2','0','π/2','π']
    y = np.array([0,90,180,270,360])
    y_ticks = ['-π','-π/2','0','π/2','π']
    
    # 画图
    fig = plt.figure(figsize=(cm2inch(4.5), cm2inch(4)))
    h = plt.axes([0.14,0.09,0.8,0.67])
    
    sns.heatmap(GScorr_CVC,cmap='RdBu_r',vmax=1,vmin=-1,
                cbar_kws={'label':'R'})
    h.spines['right'].set_visible(True)
    h.spines['left'].set_visible(True)
    h.spines['bottom'].set_visible(True)
    
    # 坐标轴
    plt.xticks(x,x_ticks,rotation=360,fontname='Times New Roman',fontsize=tick_size)
    plt.yticks(y,y_ticks,fontname='Times New Roman',fontsize=tick_size)
    plt.xlabel('相位',fontdict={'fontsize':label_size,
    'fontfamily':'SimSun'})
    plt.ylabel('相位',labelpad=-1.5,fontdict={'fontsize':label_size,
    'fontfamily':'SimSun'})
    
    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    
    # 颜色条
    cb = h.collections[0].colorbar
    
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_family('Times New Roman')
    cb.set_ticks([-1,-.5,0,.5,1])
    cb.ax.tick_params(labelsize=tick_size,length=2,pad=1)
    cb.set_label('R',fontdict={'family':'Times New Roman','size':label_size},labelpad=0)
    
    """折线图"""  
    mean_cl_r2 = np.transpose(mean_cl_r)
    l = plt.axes([0.14,0.76,0.64,0.16])
    plt.plot(mean_cl_r2,linewidth=1)
    
    # 坐标轴刻度
    plt.xlim(0,360)
    plt.xticks([])
    plt.yticks(fontname='Times New Roman',fontsize=tick_size)
    plt.ylabel('R',fontdict={'family':'Times New Roman','size':label_size},labelpad=0.5)
    plt.tick_params(axis='x',length=2,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    if i == 1:
        plt.yticks([0.2,0.3],fontname='Times New Roman',fontsize=tick_size)

    if i == 4:
        plt.yticks(np.arange(0.3,0.9,0.2),fontname='Times New Roman',fontsize=tick_size)
             
        
    # 显著性线
    min_r = 0.26*np.ones((360,1))
    plt.plot(min_r,linewidth=0.8,color='grey',linestyle='dashed')


    if i == 0:
        line1 = 65
        line2 = 114
        line3 = 249
        line4 = 300
        plt.yticks([0.3,0.5],fontname='Times New Roman',fontsize=tick_size)
        
    elif i ==1:
        line1 = 72
        line2 = 119
        line3 = 250
        line4 = 291
        
    elif i == 2:
        line1 = 64
        line2 = 114
        line3 = 250
        line4 = 297
        plt.yticks([0.3,0.6,0.9],fontname='Times New Roman',fontsize=tick_size)
    elif i == 3:
        line1 = 71
        line2 = 116
        line3 = 246
        line4 = 296
        plt.yticks([0.3,0.6,0.9],fontname='Times New Roman',fontsize=tick_size)
    elif i == 4:
        line1 = 65
        line2 = 116
        line3 = 245
        line4 = 299
        plt.yticks([0.3,0.6,0.9],fontname='Times New Roman',fontsize=tick_size)
        
        
    if i != 1 and i != 0:    
        plt.axvline(line1,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line2,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line3,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line4,color='black',lw=0.5,linestyle='-.')
    
    
    plt.title(title_name[i],
              fontdict={'fontsize':title_size,'fontfamily':'Times New Roman','weight':'bold'},
              pad=0)
    
    
    plt.sca(h)
    if i != 1 and i != 0:  
        plt.axvline(line1,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line2,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line3,color='black',lw=0.5,linestyle='-.')
        plt.axvline(line4,color='black',lw=0.5,linestyle='-.')
    
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Task_Fig1')
    # plt.savefig(img_name[i])
    plt.savefig(img_name2[i],dpi=300)
    plt.show()

############
x = np.linspace(-np.pi, np.pi, 256,endpoint=True)
C = np.cos(x)
pi = np.pi

tick_size = 6

fig = plt.figure(figsize=(cm2inch(4.5), cm2inch(4)))
ax = plt.axes([0.02,0.1,0.95,0.8])

plt.plot(x,C,linewidth=1,color='black')




ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_position(('data',0))  
ax.spines['bottom'].set_position(('data',0))  
# 坐标轴刻度

x_ticks = ['-π','-π/2','0','π/2','π']
plt.xticks([-pi,-pi/2,0,pi/2,pi],x_ticks,fontname='Times New Roman',fontsize=tick_size)
plt.yticks([-1,-0.5,0.5,1],fontname='Times New Roman',fontsize=tick_size)


os.chdir('E:\HCP\Results_PL\FiguresV0.2\Task_Fig1')
# plt.savefig(img_name[i])
plt.savefig('cosine.tif',dpi=300)
plt.show()    