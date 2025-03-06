#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:26:42 2023

@author: yujia
"""



import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py

in_path = 'For seaborn2'
out_path = 'edge_frequency'
subj_num = 178

# 导入数据集
img_name = ['phase_mean_band1','phase_mean_band2','phase_mean_band3','phase_mean_band4']


# 绘图
x_ticks = ['-π','-π/2','0','π/2','π']
x_num = [-180,-90,0,90,180]
tick_size = 16
label_size = 20

path = '/Users/yujia/Data/Data/HCP/Figures/'+in_path+'/phase_mean_dy.mat'
PV_file = scio.loadmat(path)
x = PV_file['x']

y = np.squeeze(PV_file['y'])

for i in [0,1,2,3]:
    x1 = x[:,i]
    
    
    Data = pd.DataFrame(data={'Mean frequency': x1, 'Phase angle': y})    
    ax = sns.relplot(x="Phase angle", y="Mean frequency", kind="line", ci="sd", data=Data)
    plt.xticks(x_num,x_ticks,fontname='Arial',fontsize=tick_size)
    plt.xlim(-180,180)
    plt.yticks(fontname='Arial',fontsize=tick_size)
    plt.xlabel("Phase angle",fontname='Arial',fontsize=label_size,labelpad=1)
    plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)
    sns.despine(top=False,right=False)
    ax.figure.set_size_inches(6.1,5.1)
    
    os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
    plt.savefig(img_name[i],dpi=300)
    plt.show()
    
    

# bar figure    
path = '/Users/yujia/Data/Data/HCP/Figures/'+in_path+'/figure3_bar.mat'
PV_file = scio.loadmat(path)

x2 = np.squeeze(PV_file['x2'])

band1 = subj_num*['0.005-0.01 Hz']
band2 = subj_num*['0.01-0.03 Hz']
band3 = subj_num*['0.03-0.05 Hz']
band4 = subj_num*['0.05-0.07 Hz']
band5 = subj_num*['0.07-0.09 Hz']

Frequency = band1+band2+band3+band4+band5

Data_mean = pd.DataFrame(data={'Frequency Mean': x2, 'Frequency': Frequency})

ax = plt.axes([0.14,0.09,0.8,0.67])
sns.barplot(x="Frequency", y="Frequency Mean", data=Data_mean)
ax.set_yscale('log')
ax.figure.set_size_inches(7,5) 
plt.ylim(0.000000001,0.001)
plt.xticks(fontname='Arial',fontsize=tick_size,rotation=45)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Band",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("KL divergence",fontname='Arial',fontsize=label_size,labelpad=1)

os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('mean_bar',dpi=300)
plt.show()

  
    
