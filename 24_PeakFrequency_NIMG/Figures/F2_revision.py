# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:32:20 2023

@author: Yujia Ao
"""

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



in_path = 'For seaborn'
out_path = 'line chart'
subj_num = 181

tick_size = 16
label_size = 20
imag_name = ['F_band1','F_band2','F_band3','F_band4']
# bar figure    
path = 'E:\\Phase dynamics and INT\\Figure\\'+in_path+'\\figure2_bar_revision.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x_mean_phase'])



phase1 = subj_num*['Peak']
phase2 = subj_num*['Trough']
phase3 = subj_num*['Rise']
phase4 = subj_num*['Fall']

Phase = phase1+phase2+phase3+phase4


Data_phase = pd.DataFrame(data={'Frequency Mean': x1, 'Phase': Phase})


ax = plt.axes([0.15,0.15,0.8,0.8])
sns.violinplot(x="Phase", y="Frequency Mean", data=Data_phase)


ax.figure.set_size_inches(7,5)

plt.ylim(0.218,0.226)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Phase",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)


os.chdir('E:\\Phase dynamics and INT\\Figure\\'+out_path)
plt.savefig('mean_phase',dpi=300)
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
    os.chdir('E:\\Phase dynamics and INT\\Figure\\'+out_path)
    plt.savefig(imag_name[i],dpi=300)
    
    plt.show()
    
    



import matplotlib.ticker as ticker


in_path = 'For seaborn'
out_path = 'line chart'
subj_num = 181
# 导入数据集
img_name = ['mean_scatter_band1','mean_scatter_band2','mean_scatter_band3','mean_scatter_band4',
            'mean_scatter_band5','mean_scatter_band6','mean_scatter_band7','mean_scatter_band8']


# 绘图

tick_size = 16
label_size = 20

path = 'E:\\Phase dynamics and INT\\Figure\\' + in_path + '\\figure5_scatter.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x2'])
y1 = np.squeeze(PV_file['y1'])




for i in range(0,4):
        
    x = np.squeeze(x1[:,i])
    
    ax = sns.regplot(x=x,y=y1)
   
    plt.yticks(fontname='Arial',fontsize=tick_size)
    plt.xticks(fontname='Arial',fontsize=tick_size)
    plt.xlabel("Frequency power",fontname='Arial',fontsize=label_size,labelpad=1)
    plt.ylabel("ACW",fontname='Arial',fontsize=label_size,labelpad=1)
    sns.despine(top=False,right=False)
    ax.figure.set_size_inches(5.5,5)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    os.chdir('E:\\Phase dynamics and INT\\Figure\\'+out_path)
    plt.savefig(img_name[i],dpi=300)
    plt.show()


    
path = 'E:\\Phase dynamics and INT\\Figure\\'+in_path+'\\figure5_bar.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x1'])

band1 = ['0.01-0.03 Hz']
band2 = ['0.03-0.05 Hz']
band3 = ['0.05-0.07 Hz']
band4 = ['0.07-0.09 Hz']    
Frequency = band1+band2+band3+band4
Data_mean = pd.DataFrame(data={'CRC': x1, 'Frequency': Frequency})


ax = sns.barplot(x="Frequency", y="CRC", data=Data_mean)
ax.figure.set_size_inches(7,5) 
plt.ylim(-1,0.6)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Band",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("R",fontname='Arial',fontsize=label_size,labelpad=1)
plt.axhline(y=0)
os.chdir('E:\\Phase dynamics and INT\\Figure\\'+out_path)
plt.savefig('mean_bar',dpi=300)
plt.show()    
    



