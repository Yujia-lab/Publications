#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:44:28 2023

@author: yujia
"""
import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py
from statannotations.Annotator import Annotator


in_path = 'For seaborn 7T'
out_path = 'Figure 2 7T'
subj_num = 178

tick_size = 16
label_size = 20
imag_name = ['F_band1','F_band2','F_band3','F_band4']
# bar figure    
path = '/Users/yujia/Data/Data/HCP/Figures/'+in_path+'/figure2_bar.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x_mean_phase'])
x2 = np.squeeze(PV_file['x_mean_frequency'])

band1 = subj_num*['0.01-0.03 Hz']
band2 = subj_num*['0.03-0.05 Hz']
band3 = subj_num*['0.05-0.07 Hz']
band4 = subj_num*['0.07-0.09 Hz']

Frequency = band1+band2+band3+band4

phase1 = subj_num*['Peak']
phase2 = subj_num*['Trough']
phase3 = subj_num*['Rise']
phase4 = subj_num*['Fall']

Phase = phase1+phase2+phase3+phase4

Data_frequency = pd.DataFrame(data={'Frequency Mean': x2, 'Frequency': Frequency})
Data_phase = pd.DataFrame(data={'Frequency Mean': x1, 'Phase': Phase})


ax = plt.axes([0.15,0.15,0.8,0.8])
sns.violinplot(x="Phase", y="Frequency Mean", data=Data_phase)


ax.figure.set_size_inches(7,5)

plt.ylim(0.218,0.226)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Phase",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('mean_phase',dpi=300)
plt.show()



ax = plt.axes([0.15,0.15,0.8,0.8])
sns.violinplot(x="Frequency", y="Frequency Mean", data=Data_frequency)
ax.figure.set_size_inches(7,5) 
plt.ylim(0.07,0.39)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Band",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('mean_fre',dpi=300)
plt.show()







os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('var_fre',dpi=300)
plt.show()

x1 = np.squeeze(PV_file['x_band1'])
x2 = np.squeeze(PV_file['x_band2'])
x3 = np.squeeze(PV_file['x_band3'])
x4 = np.squeeze(PV_file['x_band4'])
x = [x1,x2,x3,x4]

for i in [0,1,2,3]:
    Data = pd.DataFrame(data={'Frequency Mean': x[i], 'Phase': Phase})
    ax = plt.axes([0.15,0.15,0.8,0.8])
    sns.violinplot(x="Phase", y="Frequency Mean", data=Data)
    ax.figure.set_size_inches(7,5)
    if i==0:
        plt.ylim(0.078,0.102)
    if i==1:
        plt.ylim(0.170,0.185)
    if i==2:
        plt.ylim(0.259,0.273)
    if i==3:
        plt.ylim(0.349,0.362)
    plt.xticks(fontname='Arial',fontsize=tick_size)
    plt.yticks(fontname='Arial',fontsize=tick_size)
    plt.xlabel("Phase",fontname='Arial',fontsize=label_size,labelpad=1)
    plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)
    os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
    plt.savefig(imag_name[i],dpi=300)
    
    plt.show()
    
## 热力图
label = ['peak','trough','rise','fall']
label1 = 4*label

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/correlation_mat.mat'
PV_file = scio.loadmat(path)
x = np.squeeze(PV_file['correlation_map'])

ax = sns.heatmap(x,xticklabels=label1,yticklabels=label1,annot_kws={'fontstyle':'Arial'})

plt.rcParams['font.family'] = 'Arial'

os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('heatmap',dpi=300)

plt.show()










