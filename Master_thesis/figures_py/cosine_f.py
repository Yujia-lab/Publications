# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 16:12:12 2022

@author: Yujia Ao
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

import pandas as pd
import math

def cm2inch(value): 
    return value/2.54
 
x = np.linspace(-np.pi, np.pi, 256,endpoint=True)
C = np.cos(x)
pi = np.pi

tick_size = 6

fig = plt.figure(figsize=(cm2inch(6), cm2inch(4))) 
ax = plt.axes([0.02,0.02,0.95,0.92])


plt.plot(x[0:15],C[0:15],linewidth=3,color='#D48640')
plt.plot(x[15:50],C[15:50],linewidth=3,color='black')
plt.plot(x[50:78],C[50:78],linewidth=3,color='#539045')
plt.plot(x[78:114],C[78:114],linewidth=3,color='black')
plt.plot(x[114:143],C[114:143],linewidth=3,color='#44729D')
plt.plot(x[143:180],C[143:180],linewidth=3,color='black')
plt.plot(x[178:204],C[178:204],linewidth=3,color='#B14743')
plt.plot(x[205:242],C[205:242],linewidth=3,color='black')
plt.plot(x[242:],C[242:],linewidth=3,color='#D48640')




ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_position(('data',0))  
ax.spines['bottom'].set_position(('data',0))  
# 坐标轴刻度

x_ticks = ['-π','-π/2','0','π/2','π']
plt.xticks([-pi,-pi/2,0,pi/2,pi],x_ticks,fontname='Arial',fontsize=tick_size)
plt.yticks([-1,-0.5,0.5,1],fontname='Arial',fontsize=tick_size)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
# plt.savefig(img_name[i])
plt.savefig('cosine.tif',dpi=300)
plt.show()

 
tick_size = 6

fig = plt.figure(figsize=(cm2inch(6), cm2inch(4))) 
ax = plt.axes([0.02,0.02,0.95,0.92])


plt.plot(x,C,linewidth=3,color='black')




ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_position(('data',0))  
ax.spines['bottom'].set_position(('data',0))  
# 坐标轴刻度

x_ticks = ['-π','-π/2','0','π/2','π']
plt.xticks([-pi,-pi/2,0,pi/2,pi],x_ticks,fontname='Arial',fontsize=tick_size)
plt.yticks([-1,-0.5,0.5,1],fontname='Arial',fontsize=tick_size)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
# plt.savefig(img_name[i])
plt.savefig('cosine_2 .tif',dpi=300)
plt.show()

