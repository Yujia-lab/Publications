# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 17:00:29 2022

@author: Yujia Ao
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 16:27:14 2022

@author: Yujia Ao
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import h5py
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def cm2inch(value): 
    return value/2.54
# 读取数据
Cohere_path = 'E:\GS Coherence\Deconv\ROI_Cohere.mat'
CORR_path = 'E:\GS Coherence\Deconv\Age_GScorr.mat'
Cohere_file = scio.loadmat(Cohere_path)
CORR_file = scio.loadmat(CORR_path)

Cohere = Cohere_file['Mean']
CORR = CORR_file['Z']
CORR_FDR = CORR_file['Z_FDR']

# 整理数据
Cohere = np.transpose(Cohere)
CORR = np.transpose(CORR)
CORR_FDR = np.transpose(CORR_FDR)
CORR_FDR_abs = np.absolute(CORR_FDR)
CORR_FDR_abs_non0 = CORR_FDR_abs.nonzero()
minZ = CORR_FDR_abs[CORR_FDR_abs_non0].min()
maxZ = CORR_FDR_abs.max()
"""图A热力"""
# 画图参数
x_ticks = ['0','0.05','0.1','0.15','0.2','0.25']
x = [0,5/25*257,10/25*257,15/25*257,20/25*257,257]
tick_size = 7
label_size = 8
bwidth = 0.5
lwidth = 0.6
brain_regions = [69,125,163,175,189,211]
cbar = 'coolwarm'
hm = [0.15,0.14,0.85,0.64]
lc = [0.15,0.78,0.68,0.18]
fig_size = (cm2inch(6), cm2inch(5))
# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(Cohere,cmap=cbar,vmax=0.62,vmin=0.12,
                cbar_kws={'label':'C'},center=0.12)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([0.2,0.4,0.6])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'



"""图A折线"""
# 整理数据
Mean_Cohere = np.mean(Cohere,axis=0)

# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('C',fontname='Arial',fontsize=label_size,labelpad=0)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_Deconv.tif',dpi=300)
plt.show()

"""图B热力"""

# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(CORR,cmap=cbar,vmax=maxZ,vmin=-maxZ,
                cbar_kws={'label':'C'},center=0)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([-0.2,0,0.2])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'

"""图B折线"""
# 整理数据
Mean_Cohere = np.mean(CORR,axis=0)

# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('',fontname='Arial',fontsize=label_size,labelpad=0)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_Deconv2.tif',dpi=300)
plt.show()

"""图C热力"""
fig = plt.figure(figsize=(fig_size))
h3 = plt.axes(hm)
# 画图
sns.heatmap(CORR_FDR,cmap=cbar,cbar_kws={'label':'C'},vmin=-maxZ,vmax=maxZ)
# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)

# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 加框
h3.spines['top'].set_visible(True)
h3.spines['right'].set_visible(True)
h3.spines['left'].set_visible(True)
h3.spines['bottom'].set_visible(True)
h3.spines['top'].set_linewidth(bwidth)
h3.spines['right'].set_linewidth(bwidth)
h3.spines['left'].set_linewidth(bwidth)
h3.spines['bottom'].set_linewidth(bwidth)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)

# 颜色条字体
cb = h3.collections[0].colorbar
cb.set_ticks([-0.2,0,0.2])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=1)
plt.rcParams['font.family'] = 'Arial'

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1_Deconv3.tif',dpi=300)
plt.show()


# ####################################2
# # 读取数据
Cohere_path = 'E:\GS Coherence\Regression\ROI_Cohere.mat'
CORR_path = 'E:\GS Coherence\Regression\Age_GScorr.mat'
Cohere_file = scio.loadmat(Cohere_path)
CORR_file = scio.loadmat(CORR_path)

Cohere = Cohere_file['Mean']
CORR = CORR_file['Z']
CORR_FDR = CORR_file['Z_FDR']

# 整理数据
Cohere = np.transpose(Cohere)
CORR = np.transpose(CORR)
CORR_FDR = np.transpose(CORR_FDR)
CORR_FDR_abs = np.absolute(CORR_FDR)
CORR_FDR_abs_non0 = CORR_FDR_abs.nonzero()
minZ = CORR_FDR_abs[CORR_FDR_abs_non0].min()
maxZ = CORR_FDR_abs.max()
"""图A热力"""
# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(Cohere,cmap=cbar,vmax=0.62,vmin=0.12,
                cbar_kws={'label':'C'},center=0.12)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([0.2,0.4,0.6])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'



"""图A折线"""
# 整理数据
Mean_Cohere = np.mean(Cohere,axis=0)

# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('C',fontname='Arial',fontsize=label_size,labelpad=1)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_regression.tif',dpi=300)
plt.show()

"""图B热力"""

# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(CORR,cmap=cbar,vmax=maxZ,vmin=-maxZ,
                cbar_kws={'label':'C'},center=0)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([-0.2,0,0.2])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'

"""图B折线"""
# 整理数据
Mean_Cohere = np.mean(CORR,axis=0)
# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('',fontname='Arial',fontsize=label_size,labelpad=0)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_regression2.tif',dpi=300)
plt.show()


"""图C热力"""
fig = plt.figure(figsize=(fig_size))
h3 = plt.axes(hm)
# 画图
sns.heatmap(CORR_FDR,cmap=cbar,cbar_kws={'label':'C'},vmin=-maxZ,vmax=maxZ)
# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)

# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 加框
h3.spines['top'].set_visible(True)
h3.spines['right'].set_visible(True)
h3.spines['left'].set_visible(True)
h3.spines['bottom'].set_visible(True)
h3.spines['top'].set_linewidth(bwidth)
h3.spines['right'].set_linewidth(bwidth)
h3.spines['left'].set_linewidth(bwidth)
h3.spines['bottom'].set_linewidth(bwidth)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)

# 颜色条字体
cb = h3.collections[0].colorbar
cb.set_ticks([-0.2,0,0.2])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=1)
plt.rcParams['font.family'] = 'Arial'

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1_regression3.tif',dpi=300)
plt.show()

# ####################################3
# 读取数据
Cohere_path = 'E:\GS Coherence\Deconv_noreg\ROI_Cohere.mat'
CORR_path = 'E:\GS Coherence\Deconv_noreg\Age_GScorr.mat'
Cohere_file = scio.loadmat(Cohere_path)
CORR_file = scio.loadmat(CORR_path)

Cohere = Cohere_file['Mean']
CORR = CORR_file['Z']
CORR_FDR = CORR_file['Z_FDR']

# 整理数据
Cohere = np.transpose(Cohere)
CORR = np.transpose(CORR)
CORR_FDR = np.transpose(CORR_FDR)
CORR_FDR_abs = np.absolute(CORR_FDR)
CORR_FDR_abs_non0 = CORR_FDR_abs.nonzero()
minZ = CORR_FDR_abs[CORR_FDR_abs_non0].min()
maxZ = CORR_FDR_abs.max()

"""图A热力"""
# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(Cohere,cmap=cbar,
                cbar_kws={'label':'C'},center=0.12)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([0.2,0.4,0.6,0.8])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'



"""图A折线"""
# 整理数据
Mean_Cohere = np.mean(Cohere,axis=0)

# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('C',fontname='Arial',fontsize=label_size,labelpad=1)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_noregression.tif',dpi=300)
plt.show()

"""图B热力"""

# 画图
fig = plt.figure(figsize=(fig_size))
h1 = plt.axes(hm)
sns.heatmap(CORR,cmap=cbar,vmax=maxZ,vmin=-maxZ,
                cbar_kws={'label':'C'},center=0)

# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)
# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
          linewidth=lwidth)


# 加框
h1.spines['right'].set_visible(True)
h1.spines['left'].set_visible(True)
h1.spines['bottom'].set_visible(True)
h1.spines['right'].set_linewidth(bwidth)
h1.spines['left'].set_linewidth(bwidth)
h1.spines['bottom'].set_linewidth(bwidth)


# 颜色条字体
cb = h1.collections[0].colorbar
cb.set_ticks([-0.4,-0.2,0,0.2,0.4])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=2)
plt.rcParams['font.family'] = 'Arial'

"""图B折线"""
# 整理数据
Mean_Cohere = np.mean(CORR,axis=0)

# 画图
l1 = plt.axes(lc)
plt.plot(Mean_Cohere,linewidth=1,color='grey')
# 坐标轴标题
plt.ylabel('',fontname='Arial',fontsize=label_size,labelpad=0)
# 坐标轴刻度
plt.xticks([])
plt.yticks(fontname='Arial',fontsize=tick_size)
# 刻度线
plt.tick_params(axis='y',length=2,pad=1)
# 边框线
l1.spines['right'].set_linewidth(bwidth)
l1.spines['left'].set_linewidth(bwidth)
l1.spines['bottom'].set_linewidth(bwidth)
l1.spines['top'].set_linewidth(bwidth)

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1.pdf')
plt.savefig('Fig1_noregression2.tif',dpi=300)
plt.show()

"""图C热力"""
fig = plt.figure(figsize=(fig_size))
h3 = plt.axes(hm)
# 画图
sns.heatmap(CORR_FDR,cmap=cbar,cbar_kws={'label':'C'},vmin=-maxZ,vmax=maxZ)
# 坐标轴刻度
plt.xticks(x,x_ticks,rotation=360,fontname='Arial',fontsize=tick_size)
plt.yticks([])
plt.xlabel('Frequency (Hz)',fontdict={'family':'Arial','size':label_size},labelpad=1)

# 刻度线
plt.tick_params(axis='x',length=2,pad=1)

# 加框
h3.spines['top'].set_visible(True)
h3.spines['right'].set_visible(True)
h3.spines['left'].set_visible(True)
h3.spines['bottom'].set_visible(True)
h3.spines['top'].set_linewidth(bwidth)
h3.spines['right'].set_linewidth(bwidth)
h3.spines['left'].set_linewidth(bwidth)
h3.spines['bottom'].set_linewidth(bwidth)

# 分脑区的线
plt.plot(np.arange(0,257),brain_regions[0]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[1]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[2]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[3]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[4]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)
plt.plot(np.arange(0,257),brain_regions[5]*np.ones((257,1)),'--',color='black',
         linewidth=lwidth)

# 颜色条字体
cb = h3.collections[0].colorbar
cb.set_ticks([-0.4,-0.2,0,0.2,0.4])
cb.ax.tick_params(length=2,labelsize=tick_size,pad=1)
cb.set_label('',fontdict={'family':'Arial','size':label_size},
             labelpad=1)
plt.rcParams['font.family'] = 'Arial'

os.chdir('E:\GS Coherence\FIG_revision\Fig1')
plt.savefig('Fig1_noregression3.tif',dpi=300)
plt.show()