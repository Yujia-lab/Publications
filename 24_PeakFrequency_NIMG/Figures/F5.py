#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:02:10 2023

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

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure5_bar.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x1'])
x2 = np.squeeze(PV_file['x2'])
x3 = np.squeeze(PV_file['x3'])


rest = 35*['Rest']
task = 35*['Task']

condition = rest+task

Data_mean = pd.DataFrame(data={'condition': condition, 'value': x1})
Data_var = pd.DataFrame(data={'condition': condition, 'value': x2})
Data_acw = pd.DataFrame(data={'condition': condition, 'value': x3})

ax = sns.violinplot(x="condition", y="value", data=Data_mean)

ax.figure.set_size_inches(7,5) 
plt.ylim(-0.00002,0.00014)
plt.rcParams["font.family"] = 'Arial'
ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)

plt.xlabel(" ",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("MI",fontname='Arial',fontsize=label_size,labelpad=1)

plt.rcParams['ytick.labelsize'] = tick_size

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 5')
plt.savefig('mean_bar',dpi=300)
plt.show()



ax = sns.violinplot(x="condition", y="value", data=Data_var)

ax.figure.set_size_inches(7,5) 
plt.ylim(-0.0001,0.0005)
plt.rcParams["font.family"] = 'Arial'
ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)

plt.xlabel(" ",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("MI",fontname='Arial',fontsize=label_size,labelpad=1)

plt.rcParams['ytick.labelsize'] = tick_size

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 5')
plt.savefig('var_bar',dpi=300)
plt.show()


ax = sns.violinplot(x="condition", y="value", data=Data_acw)

ax.figure.set_size_inches(7,5) 
plt.ylim(6.4,11)
plt.rcParams["font.family"] = 'Arial'
ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)

plt.xlabel(" ",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("ACW-0",fontname='Arial',fontsize=label_size,labelpad=1)

plt.rcParams['ytick.labelsize'] = tick_size

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 5')
plt.savefig('acw_bar',dpi=300)
plt.show()

