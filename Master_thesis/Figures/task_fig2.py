# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:29:43 2022

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
from matplotlib.patches import Patch
from statannotations.Annotator import Annotator

def cm2inch(value): 
    return value/2.54
main_dir = ['E:\HCP\Results_task\slow-1',
        'E:\HCP\Results_task\slow-2',
        'E:\HCP\Results_task\slow-3',
        'E:\HCP\Results_task\slow-4',
        'E:\HCP\Results_task\slow-5'];


T0_all = []
A0_all = []
P0_all = []
D0_all = []

# 画图参数
x = np.arange(0,25,5)
width = 0.7
colors = ['#459EF9','#86A1A8','#484F5F','#F2D8BF']
x_ticks = ['Slow-3','Slow-4','Slow-5']

title_size = 8
label_size = 7
tick_size = 6

legend_elements = [Patch(facecolor=colors[0], edgecolor='black',label='波谷'),
                   Patch(facecolor=colors[1], edgecolor='black',label='波峰'),
                   Patch(facecolor=colors[2], edgecolor='black',label='波峰'),
                   Patch(facecolor=colors[3], edgecolor='black',label='下降')]
    
for i in [2,3,4]:
    d = main_dir[i]
    # 读取数据
    Eg_dir = os.path.join(d,'Efficiency.mat')
    Eg_file = scio.loadmat(Eg_dir)
    Eg = Eg_file['Eg_trough']
    T0_all = np.append(T0_all,Eg[0,:])
    Eg = Eg_file['Eg_ascend']
    A0_all = np.append(A0_all,Eg[0,:])
    Eg = Eg_file['Eg_peak']
    P0_all = np.append(P0_all,Eg[0,:])
    Eg = Eg_file['Eg_descend']
    D0_all = np.append(D0_all,Eg[0,:])

Eg = np.hstack((T0_all,A0_all,P0_all,D0_all))

Trough = 82*3*['波谷']
Rise = 82*3*['上升']
Peak = 82*3*['波峰']
Fall = 82*3*['下降']
All_phase = Trough+Rise+Peak+Fall


Slow3 = 82*['Slow-3']
Slow4 = 82*['Slow-4']
Slow5 = 82*['Slow-5']
All_frequency = 4*(Slow3+Slow4+Slow5)

Data = pd.DataFrame(data={'Phase': All_phase, 'Frequency': All_frequency,'Eg':Eg})

# 画图
fig = plt.figure(figsize=(cm2inch(9), cm2inch(6)))
ax = sns.boxplot(x='Frequency',y='Eg',hue='Phase',data=Data,palette=colors,showfliers=False,
                 width=0.8,linewidth=0.5)


plt.xlabel('频率',fontname='SimSun',fontsize=label_size,labelpad=1)
plt.ylabel('全局效率',fontname='SimSun',fontsize=label_size,labelpad=1)

plt.legend(loc=1,ncol=4,
           prop={'family' : 'SimSun','size':5},columnspacing=1,
          handletextpad=0.2,frameon=False)

box_pairs = [(('Slow-3','波峰'),('Slow-3','上升')),
             (('Slow-3','波峰'),('Slow-3','下降')),
             (('Slow-3','波谷'),('Slow-3','上升')),
             (('Slow-3','波谷'),('Slow-3','下降')),                   
             (('Slow-4','波峰'),('Slow-4','上升')),
             (('Slow-4','波峰'),('Slow-4','下降')),
             (('Slow-4','波谷'),('Slow-4','上升')),
             (('Slow-4','波谷'),('Slow-4','下降')),            
             (('Slow-5','波峰'),('Slow-5','上升')),
             (('Slow-5','波峰'),('Slow-5','下降')),
             (('Slow-5','波谷'),('Slow-5','上升')),
             (('Slow-5','波谷'),('Slow-5','下降')),
             (('Slow-5','波谷'),('Slow-5','波峰'))]
annotator =  Annotator(ax, data=Data, x="Frequency",y="Eg",hue="Phase",
                      pairs=box_pairs)
annotator.configure(test='t-test_paired', text_format='star',line_height=0.01,
                    line_width=0.5,comparisons_correction='Bonferroni',text_offset=-0.5)
annotator.apply_and_annotate()

plt.xlim(-0.7,2.7)
plt.ylim(0.1,0.55)
plt.yticks(np.arange(0.1,0.6,0.1),fontname='Times New Roman',fontsize=tick_size)
plt.xticks(range(3),x_ticks,fontname='Times New Roman',fontsize=tick_size)
plt.tick_params(axis='x',length=0,pad=1)
plt.tick_params(axis='y',length=2,pad=1)



os.chdir('E:\HCP\Results_PL\FiguresV0.2\Task_Fig2')

plt.savefig('Eg.tif',dpi=300)
plt.show()