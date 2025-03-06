# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:20:34 2023

@author: Yujia Ao
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
subj_num = 178

# 导入数据集
img_name = ['phase_mean_band1','phase_mean_band2','phase_mean_band3','phase_mean_band4']

# 绘图
x_ticks = ['-π','-π/2','0','π/2','π']
x_num = [-180,-90,0,90,180]
tick_size = 16
label_size = 20

path = 'E:\\Phase dynamics and INT\\Figure\\'+in_path+'\\phase_mean_dy.mat'
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
    
    os.chdir('E:\\Phase dynamics and INT\\Figure\\'+out_path)
    plt.savefig(img_name[i],dpi=300)
    plt.show()
    
    

  
    
