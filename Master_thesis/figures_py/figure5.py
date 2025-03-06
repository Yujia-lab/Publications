# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 18:36:14 2022

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




# 画图参数
x = np.arange(0,10,5)
width = 0.7
colors = ['#459EF9','#86A1A8','#484F5F','#F2D8BF']
y_ticks = np.arange(0,1400,200)
x_ticks = ['0-back','2-back']
title_size = 8
label_size = 7
tick_size = 6

legend_elements = [Patch(facecolor=colors[0], edgecolor='black',label='波谷'),
                   Patch(facecolor=colors[1], edgecolor='black',label='上升'),
                   Patch(facecolor=colors[2], edgecolor='black',label='波峰'),
                   Patch(facecolor=colors[3], edgecolor='black',label='下降')]

  
for ifre in [0,2,3,4]:
    
    T0_all = []
    A0_all = []
    P0_all = []
    D0_all = []
    T2_all = []
    A2_all = []
    P2_all = []
    D2_all = []

    d = main_dir[ifre]
    # 读取数据
    RT_dir = os.path.join(d,'RT_phase_ACC0.mat')
    RT_file = scio.loadmat(RT_dir)
    RT = RT_file['T0']
    T0_all = np.append(T0_all,RT[ifre,:])
    RT = RT_file['A0']
    A0_all = np.append(A0_all,RT[ifre,:])
    RT = RT_file['P0']
    P0_all = np.append(P0_all,RT[ifre,:])
    RT = RT_file['D0']
    D0_all = np.append(D0_all,RT[ifre,:])
    
    
    d = main_dir[ifre]
    # 读取数据
    RT_dir = os.path.join(d,'RT_phase_ACC2.mat')
    RT_file = scio.loadmat(RT_dir)
    RT = RT_file['T2']
    T2_all = np.append(T2_all,RT[ifre,:])
    RT = RT_file['A2']
    A2_all = np.append(A2_all,RT[ifre,:])
    RT = RT_file['P2']
    P2_all = np.append(P2_all,RT[ifre,:])
    RT = RT_file['D2']
    D2_all = np.append(D2_all,RT[ifre,:])
    
    RT = np.hstack((T0_all,A0_all,P0_all,D0_all,T2_all,A2_all,P2_all,D2_all))
    
    Trough = 82*['波谷']
    Rise = 82*['上升']
    Peak = 82*['波峰']
    Fall = 82*['下降']
    All_phase = 2*(Trough+Rise+Peak+Fall)
    
    task0 = 328*['0-back']
    task2 = 328*['2-back']
    
    All_task = task0+task2
    
    Data = pd.DataFrame(data={'Phase': All_phase, 'Task': All_task,'RT':RT})
    
    # 画图
    fig = plt.figure(figsize=(cm2inch(4.5), cm2inch(3.5)))
    plt.axes([0.18,0.16,0.8,0.76])
    ax = sns.boxplot(x='Task',y='RT',hue='Phase',data=Data,palette=colors,showfliers=False,
                     width=0.8,linewidth=0.5)
    
    
    plt.yticks(np.arange(0.75,1.01,0.05),fontname='Times New Roman',fontsize=tick_size)
    plt.xticks(range(2),x_ticks,fontname='Times New Roman',fontsize=tick_size)
    
    plt.xlabel('任务',fontname='SimSun',fontsize=label_size,labelpad=1)
    plt.ylabel('正确率',fontname='SimSun',fontsize=label_size,labelpad=1)
    
    plt.legend(loc=1,ncol=4,
               prop={'family' : 'SimSun','size':5},columnspacing=0.3,
              handletextpad=0.4,frameon=False)
    
    # if ifre == 2:
    #     box_pairs = []
    # elif ifre == 3:
    #     box_pairs = [(('0-back','波谷'),('0-back','上升')),
    #                  (('0-back','波峰'),('0-back','上升')),
    #                  (('0-back','下降'),('0-back','上升'))]
    #     annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
    #                       pairs=box_pairs)
    #     annotator.configure(test='t-test_paired', text_format='star',line_height=0.02,line_width=1.5)
    #     annotator.apply_and_annotate()
    # elif ifre == 4:
    #     box_pairs = [(('0-back','下降'),('0-back','上升'))]
    #     annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
    #                       pairs=box_pairs)
    #     annotator.configure(test='t-test_paired', text_format='star',line_height=0.02,line_width=1.5)
    #     annotator.apply_and_annotate()
    if ifre == 0:
        box_pairs = [(('0-back','波谷'),('0-back','波峰')),
                 (('0-back','上升'),('0-back','波峰')),
                 (('0-back','波峰'),('0-back','下降')),
                 (('2-back','波谷'),('2-back','波峰')),
                 (('2-back','上升'),('2-back','波峰')),
                 (('2-back','波峰'),('2-back','下降'))]
        annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
                               pairs=box_pairs)
        annotator.configure(test='t-test_paired', text_format='star',line_height=0.01,
                    line_width=0.5,comparisons_correction='Bonferroni',text_offset=-0.5)
        annotator.apply_and_annotate()
    if ifre == 2:
        box_pairs = [(('0-back','波谷'),('0-back','上升')),
                 (('0-back','波谷'),('0-back','波峰')),
                 (('0-back','波谷'),('0-back','下降')),
                 (('0-back','上升'),('0-back','波峰')),
                 (('0-back','波峰'),('0-back','下降')),
                 (('2-back','波谷'),('2-back','上升')),
                 (('2-back','波谷'),('2-back','波峰')),
                 (('2-back','波谷'),('2-back','下降')),
                 (('2-back','上升'),('2-back','波峰')),
                 (('2-back','波峰'),('2-back','下降'))]
        annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
                               pairs=box_pairs)
        annotator.configure(test='t-test_paired', text_format='star',line_height=0.01,
                    line_width=0.5,comparisons_correction='Bonferroni',text_offset=-0.5)
        annotator.apply_and_annotate()
    elif ifre == 3:
        box_pairs = [(('0-back','波谷'),('0-back','上升')),
             (('0-back','波谷'),('0-back','下降')),
             (('0-back','上升'),('0-back','波峰')),
             (('0-back','波峰'),('0-back','下降')),
             (('2-back','波谷'),('2-back','上升')),
             (('2-back','波谷'),('2-back','下降')),
             (('2-back','上升'),('2-back','波峰')),
             (('2-back','波峰'),('2-back','下降'))]
        annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
                               pairs=box_pairs)
        annotator.configure(test='t-test_paired', text_format='star',line_height=0.01,
                    line_width=0.5,comparisons_correction='Bonferroni',text_offset=-0.5)
        annotator.apply_and_annotate()
    elif ifre == 4:
        box_pairs = [(('0-back','波谷'),('0-back','上升')),
             (('0-back','波谷'),('0-back','下降')),
             (('0-back','上升'),('0-back','波峰')),
             (('0-back','波峰'),('0-back','下降')),
             (('2-back','波谷'),('2-back','上升')),
             (('2-back','波谷'),('2-back','波峰')),
             (('2-back','波谷'),('2-back','下降')),
             (('2-back','上升'),('2-back','波峰')),
             (('2-back','波峰'),('2-back','下降'))]
        annotator =  Annotator(ax, data=Data, x='Task',y="RT",hue="Phase",
                               pairs=box_pairs)
        annotator.configure(test='t-test_paired', text_format='star',line_height=0.01,
                    line_width=0.5,comparisons_correction='Bonferroni',text_offset=-0.5)
        annotator.apply_and_annotate()
        
    
    plt.ylim(0.8,1.15)
    if ifre == 0:
        plt.ylim(0.75,1.15)
    plt.xlim(-0.5,1.5)
    plt.tick_params(axis='x',length=0,pad=1)
    plt.tick_params(axis='y',length=2,pad=1)
    plt.title('Slow-'+str(ifre+1),fontdict={'fontsize':title_size,'fontfamily':'Times New Roman','weight':'bold'},
              pad=1)
    
    os.chdir('E:\HCP\Results_PL\FiguresV0.2\Fig5')

    plt.savefig('ACC_Slow-'+str(ifre+1)+'.tif',dpi=300)
    plt.show()




