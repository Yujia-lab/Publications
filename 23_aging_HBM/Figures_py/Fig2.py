# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 16:57:47 2022

@author: Yujia Ao
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 17:18:32 2021

@author: Yujia Ao
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py

def cm2inch(value): 
    return value/2.54
def set_para1():
    plt.sca(g.ax_joint)
    plt.xlabel('C',fontsize=label_size,fontname='Arial')
    plt.ylabel('Z',fontsize=label_size,fontname='Arial')
    plt.xticks(fontsize=tick_size,fontname='Arial')
    plt.yticks(fontsize=tick_size,fontname='Arial')
    g.fig.set_figwidth(cm2inch(5.5))
    g.fig.set_figheight(cm2inch(4.5))
    plt.tick_params(axis='x',length=3,pad=1)
    plt.tick_params(axis='y',length=3,pad=1)
    plt.sca(g.ax_marg_y)
    plt.tick_params(length=3,pad=1)
    plt.sca(g.ax_marg_x)
    plt.tick_params(length=3,pad=1)
    
def set_para2():
    plt.xticks(my_xticks,fontname='Arial',fontsize=tick_size)
    plt.yticks(fontname='Arial',fontsize=tick_size)
    plt.xlabel('C',fontname='Arial',fontsize=label_size,labelpad=1)
    plt.ylabel('Z',fontname='Arial',fontsize=label_size,labelpad=1)
    plt.tick_params(axis='both',length=3,pad=1)
   
label_size = 8
tick_size = 7

"""时间去分化的Hexbin图"""
# 整理数据
Cohere = scio.loadmat('E:\GS Coherence\Deconv\ROI_Cohere.mat')
CORR = scio.loadmat('E:\GS Coherence\Deconv\Age_GScorr.mat')
spatial_dediff = scio.loadmat('E:\GS Coherence\Deconv\spatial_dediff.mat')
C = np.transpose(Cohere['Mean'])
R = np.transpose(CORR['Z'])
R_spatial = np.squeeze(spatial_dediff['R_tempfdr'])
x = C.reshape((257*246,))
y = R.reshape((257*246,))

# 区分正负
C_pos = C[R_spatial>0,:]
C_neg = C[R_spatial<0,:]
x2 = C_pos.reshape(257*np.size(C_pos,0),)
x3 = C_neg.reshape(257*np.size(C_neg,0),)
R_pos = R[R_spatial>0,:]
R_neg = R[R_spatial<0,:]
y2 = R_pos.reshape(257*np.size(R_pos,0),)
y3 = R_neg.reshape(257*np.size(R_neg,0),)

# 画图
sns.set(style="ticks")
g = sns.jointplot(x, y, kind="hex", color="b",xlim=[0.1,0.65],ylim=[-0.45,0.35],
            joint_kws={'gridsize':15})
sns.regplot(x, y, ax=g.ax_joint, scatter=False,color='r' )
set_para1()
g.fig.set_figwidth(cm2inch(5.5))
g.fig.set_figheight(cm2inch(4.5))
os.chdir('E:\GS Coherence\FIG_revision\Fig2')
plt.savefig('temp_dediff1.tif',dpi=300)
plt.show()

sns.set(style="ticks")
g = sns.jointplot(x2, y2, kind="hex", color="orange",xlim=[0.1,0.65],ylim=[-0.45,0.35],
            joint_kws={'gridsize':15})
sns.regplot(x2, y2, ax=g.ax_joint, scatter=False,color='r' )
set_para1()
plt.savefig('temp_dediff2.tif',dpi=300)
plt.show()

sns.set(style="ticks")
g = sns.jointplot(x3, y3, kind="hex", color="b",xlim=[0.1,0.65],ylim=[-0.45,0.35],
            joint_kws={'gridsize':15})
sns.regplot(x3, y3, ax=g.ax_joint, scatter=False,color='r' )
set_para1()
plt.savefig('temp_dediff3.tif',dpi=300)
plt.show()


"""空间去分化的折线图"""
#整理数据
spatial_dediff = scio.loadmat('E:\GS Coherence\Deconv\spatial_dediff.mat')
fp = h5py.File('E://GS Coherence//result//f.mat')
R_betw = np.squeeze(spatial_dediff['R_spatial'][:])
R_FDR = np.squeeze(spatial_dediff['R_spatialfdr'][:])
f = np.squeeze(fp['f'][:])

C1 = np.mean(C[:,:15],1)
Z1 = np.mean(R[:,:15],1)
C2 = np.mean(C[:,146:],1)
Z2 = np.mean(R[:,146:],1)

# 画图

xtick = ['0','0.05','0.1','0.15','0.2','0.25']
fig = plt.figure(figsize=(cm2inch(5.5), cm2inch(4.5)))
h = plt.axes([0.18,0.16,0.77,0.8])
sns.lineplot(f,R_betw)
plt.xlim((0, 0.25))
plt.plot(f[:],np.ones(257,)*0.145,color='r',linewidth=0.5,linestyle='--')
plt.plot(f[:],np.ones(257,)*(-0.145),color='r',linewidth=0.5,linestyle='--')
my_xticks1=np.arange(0,0.26,0.05)
plt.xticks(my_xticks1,xtick,fontname='Arial',fontsize=tick_size)
my_yticks=np.arange(-0.6,0.3,0.2)
plt.yticks(my_yticks,fontname='Arial',fontsize=tick_size)
plt.xlabel('Frequency (Hz)',fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel('R',fontname='Arial',fontsize=label_size,labelpad=1)
plt.tick_params(axis='both',length=3,pad=1)
os.chdir('E:\GS Coherence\FIG_revision\Fig2')
plt.savefig('spatial_dediff1.tif',dpi=300)
plt.show()

## 图B2
fig = plt.figure(figsize=(cm2inch(5.5), cm2inch(4.5)))
h = plt.axes([0.18,0.16,0.77,0.8])
sns.regplot(C1, Z1, color='b', 
            scatter_kws={'s':10},
            line_kws={"color": "r"})
plt.xlim(0.22,0.64)
plt.ylim(-0.26,0.17)
my_xticks = np.arange(0.2, 0.61, 0.1)
set_para2()

plt.savefig('spatial_dediff2.tif',dpi=300)
plt.show()


## 图B3 
fig = plt.figure(figsize=(cm2inch(5.5), cm2inch(4.5)))
h = plt.axes([0.18,0.16,0.77,0.8])
sns.regplot(C2, Z2, color='b', 
            scatter_kws={'s':10},
            line_kws={"color": "r"})
#设置坐标轴范围
plt.xlim(0.12,0.58)
plt.ylim(-0.25,0.15)
#设置坐标轴刻度
my_xticks = np.arange(0.2, 0.51, 0.1)

set_para2()
plt.savefig('spatial_dediff3.tif',dpi=300)
plt.show()