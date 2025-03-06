#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 26 16:53:47 2023

@author: yujia
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py

tick_size = 14
label_size = 16

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure7_bar.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x1'])

band1 = ['0.01-0.03 Hz']
band2 = ['0.03-0.05 Hz']
band3 = ['0.05-0.07 Hz']
band4 = ['0.07-0.09 Hz']    
Frequency = band1+band2+band3+band4
Data_mean = pd.DataFrame(data={'R': x1, 'Frequency': Frequency})


ax = sns.barplot(x="Frequency", y="R", data=Data_mean)
ax.figure.set_size_inches(7,5) 
plt.ylim(-1,0)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Frequnecy",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("CRC",fontname='Arial',fontsize=label_size,labelpad=1)
os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 7')
plt.savefig('mean_bar',dpi=300)
plt.show()