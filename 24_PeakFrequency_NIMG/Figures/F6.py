#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:16:11 2023

@author: yujia
"""


import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py
import matplotlib.ticker as ticker
# 绘图

tick_size = 12
label_size = 14
sub_num = 179

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure6_kl.mat'
PV_file = scio.loadmat(path)
kl_data = np.squeeze(PV_file['x'])

band1 = sub_num*['0.01-0.03 Hz']
band2 = sub_num*['0.03-0.05 Hz']
band3 = sub_num*['0.05-0.07 Hz']
band4 = sub_num*['0.07-0.09 Hz']    
Frequency = 2*(band1+band2+band3+band4)

rest = 4*sub_num*['Rest']
task = 4*sub_num*['Task']
state = rest+task

Data = pd.DataFrame(data={'KL divergence': kl_data, 'Frequency': Frequency,'state':state})

ax = sns.barplot(x='Frequency',y='KL divergence',data=Data,hue='state')
ax.set_yscale('log')
plt.ylim(0.000000001,0.0001)
ax.legend().set_title('')
font1 = {'family' : 'Arial',
'weight' : 'normal',
'size'   : tick_size}
plt.legend(prop=font1)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Band",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel('KL divergence',fontname='Arial',fontsize=label_size,labelpad=1)

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 6')
plt.savefig('KL_bar',dpi=300)
plt.show()

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure6_fs.mat'
PV_file = scio.loadmat(path)
fs_data = np.squeeze(PV_file['x'])

phase1 = 2*sub_num*['Peak']
phase2 = 2*sub_num*['Trough']
phase3 = 2*sub_num*['Rise']
phase4 = 2*sub_num*['Fall']

Phase = phase1+phase2+phase3+phase4
rest = sub_num*['Rest']
task = sub_num*['Task']
state = 4*(rest+task)

Data = pd.DataFrame(data={'FS': fs_data, 'Phase': Phase,'state':state})

ax = sns.barplot(x='Phase',y='FS',data=Data,hue='state')
plt.ylim(0.08,0.095)
ax.legend().set_title('')
font1 = {'family' : 'Arial',
'weight' : 'normal',
'size'   : tick_size}
plt.legend(prop=font1)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Phase",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel('PF',fontname='Arial',fontsize=label_size,labelpad=1)

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 6')
plt.savefig('FS_bar',dpi=300)
plt.show()

# KL的小图
in_path = 'For seaborn 7T task'
out_path = 'Figure 3 7T task'
subj_num = 179

# 导入数据集
img_name = ['phase_mean_band1','phase_mean_band2','phase_mean_band3','phase_mean_band4']


# 绘图

x_ticks = ['-π','-π/2','0','π/2','π']
x_num = [-180,-90,0,90,180]

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn 7T/phase_mean_dy.mat'
PV_file = scio.loadmat(path)
x_rest = PV_file['x']

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn 7T task/phase_mean_dy.mat'
PV_file = scio.loadmat(path)
x_task = PV_file['x']

y = np.squeeze(PV_file['y'])


x1 = x_rest[:,0]
x2 = x_task[:,0]

x = [x1+x2]
x = np.concatenate((x1,x2))

y2 = np.concatenate((y,y))

rest = 64440*['Rest']
task = 64440*['Task']
state = rest+task



Data = pd.DataFrame(data={'Mean frequency': x, 'Phase angle': y2,'state':state})      
ax=sns.lineplot(x="Phase angle", y="Mean frequency", hue='state',ci=0,legend=False,data=Data)

plt.xticks(x_num,x_ticks,fontname='Arial',fontsize=tick_size)

plt.xlim(-180,180)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Phase angle",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("FS",fontname='Arial',fontsize=label_size,labelpad=1)
sns.despine(top=False,right=False)
ax.figure.set_size_inches(6,5)




os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 6')
plt.savefig('phase_mean',dpi=300)
plt.show()








