#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:32:36 2023

@author: yujia
"""


import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py

y_ticks = ['-π','-π/2','0','π/2','π']
tick_size = 16
label_size = 20
marker_size = 4
line_size = 1
y = np.arange(-180,181,90)
y_ticks = ['-π','-π/2','0','π/2','π']

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/original_data.mat'
PV_file = scio.loadmat(path)
original_data = np.squeeze(PV_file['x'])


path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/filtered_data.mat'
PV_file = scio.loadmat(path)
filtered_data = np.squeeze(PV_file['data'])
time = np.squeeze(PV_file['time'])

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/Phase_data.mat'
PV_file = scio.loadmat(path)
phase_data = np.squeeze(PV_file['data_phase2'])

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/diff_data.mat'
PV_file = scio.loadmat(path)
diff_data = np.squeeze(PV_file['data_diff'])
filt_diff_data = np.squeeze(PV_file['filt_data_diff'])

# original data
ax = sns.lineplot(x=time, y=original_data)
ax.figure.set_size_inches(8,4) 

plt.xticks([])
plt.yticks([])

plt.ylabel("Full frequency",fontname='Arial',fontsize=label_size,labelpad=1)
os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('original',dpi=300)
plt.show()

# filtered data
ax = sns.lineplot(x=time, y=filtered_data)
ax.figure.set_size_inches(8,4) 

plt.xticks([])
plt.yticks([])
plt.ylabel("0.01-0.03 Hz",fontname='Arial',fontsize=label_size,labelpad=1)
os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('filtered',dpi=300)
plt.show()

# phase data
ax = sns.scatterplot(x=time[0:1180], y=phase_data)
ax.figure.set_size_inches(8,4) 

plt.xticks([])
plt.yticks([-180,-90,0,90,180],y_ticks,fontname='Arial',fontsize=tick_size)
plt.ylabel("Phase angle",fontname='Arial',fontsize=label_size,labelpad=1)
os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('phase',dpi=300)
plt.show()


# sliding frequency
ax = sns.lineplot(x=time[0:1179], y=diff_data,color='grey',alpha=0.6,lw=2)
ax = sns.lineplot(x=time[0:1179], y=filt_diff_data,lw=3)
ax.figure.set_size_inches(8,4) 

plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.ylabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)
os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('diff',dpi=300)
plt.show()


# step 5
Peak = [];
x_Peak = [];
all_phase=[]
x_all=[]
x1 = np.array([63])
x2 = []

for i in np.arange(0,1180,1):
    if phase_data[i] > -15 and phase_data[i] < 15:
        Peak=np.append(Peak,phase_data[i])
        x_Peak=np.append(x_Peak,i)
    else:
        all_phase=np.append(all_phase,phase_data[i])
        x_all=np.append(x_all,i)
        
for i in range(1,100):
    if x_Peak[i]-x_Peak[i-1]>2:
        x1 = np.append(x1,x_Peak[i])
        x2 = np.append(x2,x_Peak[i-1])
        
x2 = np.append(x2,1180)        
        
# 画图
fig = plt.figure(figsize=(8, 4))
ax = plt.axes()
ax.plot(filt_diff_data,color='black')

ax2 = ax.twinx()
ax2.plot(x_Peak,Peak,'o',label = 'Peak',markersize=marker_size)
ax2.plot(x_all,all_phase,'o',markersize=marker_size-1,alpha=0.5)

plt.sca(ax)
ax.set_ylim(-0.05,0.55)
plt.xlim(0,1180)
plt.yticks(np.arange(0,0.21,0.05),fontname='Arial',fontsize=tick_size)

plt.xticks([],fontname='Arial',fontsize=tick_size)

plt.sca(ax2)
ax2.yaxis.set_ticks_position('left')
plt.yticks(y,y_ticks,fontname='Arial',fontsize=tick_size)
ax2.set_ylim(-600,200)
plt.tick_params(axis='y',length=2)

for i in range(0,17):
    plt.axvspan(xmin=x1[i], xmax=x2[i], facecolor="grey", alpha=0.5)

plt.plot(np.arange(1,1181,1),-220*np.ones(1180,),'--',linewidth=1,color='black')
plt.sca(ax2)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('frequency mean',dpi=300)
plt.show()
      


# for figure 2
Peak = [];
x_Peak = [];
x_Trough = [];
Trough = [];
Rise = [];
x_Rise = [];
Fall = [];
x_Fall = [];
all_phase=[]
x_all=[]


for i in np.arange(0,1180,1):
    if phase_data[i] > -15 and phase_data[i] < 15:
        Peak=np.append(Peak,phase_data[i])
        x_Peak=np.append(x_Peak,i)
    if phase_data[i] > 165 or phase_data[i] < -165:
        Trough=np.append(Trough,phase_data[i])
        x_Trough=np.append(x_Trough,i)
    if phase_data[i] > -105 and phase_data[i] < -75:
        Rise=np.append(Rise,phase_data[i])
        x_Rise=np.append(x_Rise,i)
    if phase_data[i] > 75 and phase_data[i] < 105:
        Fall=np.append(Fall,phase_data[i])
        x_Fall=np.append(x_Fall,i)
    else:
        all_phase=np.append(all_phase,phase_data[i])
        x_all=np.append(x_all,i)

x1_Peak = []
x2_Peak = []
x1_Trough = []
x2_Trough = []
x1_Rise = []
x2_Rise = []
x1_Fall = []
x2_Fall = []
        
for i in range(0,len(x_Peak)-1):
    if x_Peak[i]-x_Peak[i-1]>2:
        x1_Peak = np.append(x1_Peak,x_Peak[i])
        x2_Peak = np.append(x2_Peak,x_Peak[i-1])
for i in range(1,len(x_Trough)-1):
    if x_Trough[i]-x_Trough[i-1]>2:
        x1_Trough = np.append(x1_Trough,x_Trough[i])
        x2_Trough = np.append(x2_Trough,x_Trough[i-1])
for i in range(1,len(x_Rise)-1):
    if x_Rise[i]-x_Rise[i-1]>2:
        x1_Rise = np.append(x1_Rise,x_Rise[i])
        x2_Rise = np.append(x2_Rise,x_Rise[i-1])
for i in range(1,len(x_Fall)-1):
    if x_Fall[i]-x_Fall[i-1]>2:
        x1_Fall = np.append(x1_Fall,x_Fall[i])
        x2_Fall = np.append(x2_Fall,x_Fall[i-1])
        
x2 = np.append(x2,1180)        
        
# 画图
fig = plt.figure(figsize=(8, 2))

ax2 = plt.axes([0.06,0.13,0.8,0.85])
ax2.plot(x_Peak,Peak,'o',label = 'Peak',markersize=marker_size)
ax2.plot(x_Trough,Trough,'o',label = 'Trough',markersize=marker_size)
ax2.plot(x_Rise,Rise,'o',label = 'Rise',markersize=marker_size)
ax2.plot(x_Fall,Fall,'o',label = 'Fall',markersize=marker_size)
ax2.plot(x_all,all_phase,'o',color='grey',markersize=marker_size-1,alpha=0.5)



plt.xticks([],fontname='Arial',fontsize=tick_size)


plt.yticks(y,y_ticks,fontname='Arial',fontsize=tick_size)



# 图例
plt.sca(ax)
plt.legend(bbox_to_anchor=(1, 0.3),prop={'family' : 'Arial','size':12},
           frameon=False)
plt.sca(ax2)
plt.legend(loc=2,bbox_to_anchor=(0.97, 0.94),ncol=1,prop={'family' : 'Arial','size':tick_size},columnspacing=0.5,
          handletextpad=0,frameon=False)




os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('For F21',dpi=300)
plt.show()

fig = plt.figure(figsize=(8, 2))

ax = plt.axes([0.06,0.13,0.8,0.85])
x = x_Peak.astype(int)
ax.plot(x_Peak,filt_diff_data[x_Peak.astype(int)],'o',markersize=marker_size)
ax.plot(x_Trough,filt_diff_data[x_Trough.astype(int)],'o',markersize=marker_size)
ax.plot(x_Rise,filt_diff_data[x_Rise.astype(int)],'o',markersize=marker_size)
ax.plot(x_Fall,filt_diff_data[x_Fall.astype(int)-1],'o',markersize=marker_size)
ax.plot(x_all,filt_diff_data[x_all.astype(int)],'o',markersize=marker_size,color='grey',alpha=0.3)

plt.xticks(fontname='Arial',fontsize=tick_size)


plt.yticks(fontname='Arial',fontsize=tick_size)

os.chdir('/Users/yujia/Data/Data/HCP/Figures/Schematic figure')
plt.savefig('For F22',dpi=300)
plt.show()

