# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 15:05:46 2021

@author: Yujia Ao
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import math

def cm2inch(value): 
    return value/2.54
"""滑窗相位法示意图"""
# 读取GS数据
GS_dir = 'E:\HCP\Results_PL\slow-5\GS.mat'
phase_dir = 'E:\HCP\Results_PL\slow-5\GS_phase.mat'

GS_file = h5py.File(GS_dir)
phase_file = h5py.File(phase_dir)
GS = GS_file['GS'][:]
phase = phase_file['GS_phase'][:]

# 读取ROI数据
ROI_dir = 'E:\HCP\Results_PL\slow-5\ROISignals.mat';
ROI_file = h5py.File(ROI_dir)
ROI = ROI_file['ROISignals'][:]
ROI1 = np.squeeze(ROI[0,:,0])

# 调整数据
GS1 = GS[0,:]
pi = math.pi
phase1 = phase[0,:]/pi*180
x_phase = np.arange(0,1200,1)

y = np.arange(-180,181,90)
y_ticks = ['-π','-π/2','0','π/2','π']

trough = [];
x_trough = [];
Rise = [];
x_Rise = [];
Peak = [];
x_Peak = [];
Fall = [];
x_Fall = [];
Peak_GS = []
Peak_ROI = []

for i in np.arange(0,1200,1):
    if phase1[i] < -135 or phase1[i] > 135:
        trough=np.append(trough,phase1[i])
        x_trough=np.append(x_trough,i)
    if phase1[i] > -135 and phase1[i] < -45:
        Rise=np.append(Rise,phase1[i])
        x_Rise=np.append(x_Rise,i)
    if phase1[i] > -45 and phase1[i] < 45:
        Peak=np.append(Peak,phase1[i])
        Peak_GS = np.append(Peak_GS,GS1[i])
        Peak_ROI = np.append(Peak_ROI,ROI1[i])
        x_Peak=np.append(x_Peak,i)
    if phase1[i] > 45 and phase1[i] < 135:
        Fall=np.append(Fall,phase1[i])
        x_Fall=np.append(x_Fall,i)
 
x1 = np.array([63])
x2 = []        
for i in range(1,317):
    if x_Peak[i]-x_Peak[i-1]>2:
        x1 = np.append(x1,x_Peak[i])
        x2 = np.append(x2,x_Peak[i-1])
        
x2 = np.append(x2,1180)
        

# 画图参数
colors = ['#459EF9','#86A1A8','#484F5F','#F2D8BF']
label_size = 7
tick_size = 6
marker_size = 2
line_size = 1.5
# 画图
fig = plt.figure(figsize=(cm2inch(14), cm2inch(5)))
ax = plt.axes([0.06,0.13,0.87,0.85])
ax.plot(GS1,color='black',linewidth=line_size)

ax2 = ax.twinx()
ax2.plot(x_trough,trough,'o',label = '波谷',color=colors[0],markersize=marker_size)
ax2.plot(x_Rise,Rise,'o',label = '上升',color=colors[1],markersize=marker_size)
ax2.plot(x_Peak,Peak,'o',label = '波峰',color=colors[2],markersize=marker_size)
ax2.plot(x_Fall,Fall,'o',label = '下降',color=colors[3],markersize=marker_size)
        
# 设置刻度和字体
plt.sca(ax)
ax.set_ylim(-30,100)
plt.yticks(range(-20,21,10),fontname='Times New Roman',fontsize=tick_size)
plt.xlim(0,1200)
plt.xticks(fontname='Times New Roman',fontsize=tick_size)

plt.sca(ax2)
ax2.yaxis.set_ticks_position('left')
plt.yticks(y,y_ticks,fontname='Times New Roman',fontsize=tick_size)
ax2.set_ylim(-600,200)
plt.tick_params(axis='y',length=2)

# 坐标轴label
plt.sca(ax)
plt.xlabel('时间点',fontname='SimSun',fontsize=label_size,labelpad=1)
plt.text(-80, -5, "信号", fontname='SimSun',size = label_size,rotation=90)
plt.text(-80, 60, "相位", fontname='SimSun',size = label_size,rotation=90)
plt.tick_params(axis='x',length=2)
plt.tick_params(axis='y',length=2)

# 阴影和虚线
for i in range(0,16):
    plt.axvspan(xmin=x1[i], xmax=x2[i], facecolor="grey", alpha=0.5)

plt.plot(np.arange(1,1201,1),27*np.ones(1200,),'--',linewidth=1,color='black')
plt.sca(ax2)



# 图例
plt.sca(ax)
plt.legend(bbox_to_anchor=(1, 0.3),prop={'family' : 'SimSun','size':tick_size},
           frameon=False)
plt.sca(ax2)
plt.legend(loc=2,bbox_to_anchor=(0.98, 0.95),ncol=1,prop={'family' : 'SimSun','size':tick_size},columnspacing=0.5,
          handletextpad=0.1,frameon=False)

os.chdir('E:\HCP\Results_PL\Figures_new\PSWFC')
plt.savefig('PSWFC.tif',dpi=300)
plt.show()

"""新的时间序列"""
# GS
fig = plt.figure(figsize=(cm2inch(5), cm2inch(2)))
ax = fig.add_subplot(111)
plt.plot(Peak_GS,linewidth=1,color='black')

plt.xticks([])
plt.yticks([])

ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('GS波峰信号',fontname='SimSun',fontsize=label_size,pad=-0.5)

os.chdir('E:\HCP\Results_PL\Figures_new\PSWFC')
plt.savefig('GS_peak.tif',dpi=300)
plt.show()

# ROI
fig = plt.figure(figsize=(cm2inch(5), cm2inch(2)))
ax = fig.add_subplot(111)
plt.plot(Peak_ROI,linewidth=1,color='black')

plt.xticks([])
plt.yticks([])

ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('ROI信号',fontname='SimSun',fontsize=label_size,pad=0.8)

os.chdir('E:\HCP\Results_PL\Figures_new\PSWFC')
plt.savefig('ROI_peak.tif',dpi=300)
plt.show()



