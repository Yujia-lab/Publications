# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 10:54:17 2021

@author: Yujia Ao
"""

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import h5py
import pandas as pd

"""人口学数据bar图"""
# 数据
x = np.arange(6)
y = np.arange(0, 130, 20)
male = [50,15,12,25,19,7]
female = [77,14,34,35,21,13]
#年龄范围label
age_range = ['19-29','30-39','40-49','50-59','60-69','70-80']
bar_width = 0.6

# 创建图像
plt.figure(figsize=(4,2.5))
plt.bar(x,male,bar_width,align='center',color='#6495ED',label='male')
plt.bar(x,female,bar_width,bottom=male,align='center',color='#FFA500',label='female')

# 显示label
plt.xlabel('Age range (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : 9})
plt.ylabel('Count',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : 9})
plt.xticks(x,age_range,fontname='Arial',fontsize=7)
plt.yticks(y,fontname='Arial',fontsize=7)

plt.legend(prop={"family" : "Arial",'size':9},
                    fontsize=9)
os.chdir('E:/lifespan/time_diff/figure_revision')
plt.savefig('democratic.pdf', bbox_inches = 'tight')

plt.show()



"""Fig1 ABC（原始数据，十岁一段功率谱，功率谱重心的散点图"""
# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec.mat'
specdeconv_path = 'E:/lifespan/time_diff/original/figure/mean_spec_deconv.mat'
specnonC_path = 'E:/lifespan/time_diff/original/figure/mean_spec_nonC.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'
f_path2 = 'E:/lifespan/time_diff/original/figure/f_original.mat'
orig_path = 'E:/lifespan/time_diff/original/spec.mat'
origdeconv_path = 'E:/lifespan/time_diff/original/spec_Deconv.mat'
orignonC_path = 'E:/lifespan/time_diff/original/spec_nonC.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])

spec_file = h5py.File(spec_path)
mean_osci = np.squeeze(spec_file['mean_osci'][:])

specdeconv_file = h5py.File(specdeconv_path)
mean_osci_deconv = np.squeeze(specdeconv_file['mean_osci'][:])

specnonC_file = h5py.File(specnonC_path)
mean_osci_nonC = np.squeeze(specnonC_file['mean_osci'][:])

mean_spec = np.squeeze(spec_file['mean_spec'][:])
mean_spec_deconv = np.squeeze(specdeconv_file['mean_spec'][:])
mean_spec_nonC = np.squeeze(specnonC_file['mean_spec'][:])

f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])
f_file2 = h5py.File(f_path2)
f2 = np.squeeze(f_file2['f'][:])



cmap_path = 'E:/GS Coherence/Colormap/jet.csv'
cmap = pd.read_csv(cmap_path)
cmap2 = np.array(cmap)


#字体大小
axis_numsize = 7
label_size = 9
rp_size = 7

#画图
plt.figure(figsize=(10,2.5))

A = plt.axes([0.0,0.1,0.26,0.9]) 
B = plt.axes([0.33,0.1,0.26,0.9])
C = plt.axes([0.66,0.1,0.26,0.9])

# 图A
plt.sca(A)
for i in range(0,6):
    plt.plot(f[0:248],mean_spec[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(0, 26, 5)

plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yscale('log')
A.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

A.set_ylabel('Power',labelpad=3)
A.set_xlabel('Frequency (Hz)',labelpad=2)

legend = plt.legend(age_range, 
                    loc=1, ncol=3 , 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':7},
                    fontsize=6,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 8)

# 图B
plt.sca(B)
for i in range(0,6):
    plt.plot(f[0:248],mean_spec_deconv[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
plt.yscale('log')
my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(0.1, 0.21, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
B.tick_params(length=2)


plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
B.set_ylabel('Power',labelpad=3)
B.set_xlabel('Frequency (Hz)',labelpad=2)

legend = plt.legend(age_range, 
                    loc=1, ncol=3 , 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':7},
                    fontsize=6,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 8)

# 图C
plt.sca(C)
for i in range(0,6):
    plt.plot(f[0:248],mean_spec_nonC[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
my_x_ticks = np.arange(0, 0.26, 0.05)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(fontname='Arial',fontsize=axis_numsize)
plt.yscale('log')
B.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
B.set_ylabel('Power',labelpad=3)
B.set_xlabel('Frequency (Hz)',labelpad=2)

legend = plt.legend(age_range, 
                    loc=1, ncol=3 , 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':7},
                    fontsize=6,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 8)



os.chdir('E:/lifespan/time_diff/figure_revision')
plt.savefig('Fig1_RGB.pdf', bbox_inches = 'tight')
plt.show()

"""Fig2"""
# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec.mat'
SC_path = 'E:/lifespan/time_diff/original/figure/SC.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
spec_file = h5py.File(spec_path)
mean_osci = np.squeeze(spec_file['mean_osci'][:])
osci_r = np.squeeze(spec_file['osci_r'][:])
SC_file = h5py.File(SC_path)
SC1 = np.squeeze(SC_file['SC1'][:])
SC2 = np.squeeze(SC_file['SC2'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

cmap_path = 'E:/GS Coherence/Colormap/jet.csv'
cmap = pd.read_csv(cmap_path)
cmap2 = np.array(cmap)


#字体大小
axis_numsize = 6
label_size = 7
rp_size = 6


#画图
plt.figure(figsize=(8,6.5))



A = plt.axes([0.02,0.7,0.195,0.22]) 
B = plt.axes([0.28,0.7,0.195,0.22])
C = plt.axes([0.54,0.7,0.195,0.22]) 
D = plt.axes([0.8,0.7,0.195,0.22]) 
E = plt.axes([0.02,0.4,0.195,0.22])
F = plt.axes([0.28,0.4,0.195,0.22]) 
G = plt.axes([0.54,0.4,0.195,0.22]) 
H = plt.axes([0.8,0.4,0.195,0.22]) 
I = plt.axes([0.02,0.1,0.195,0.22])
J = plt.axes([0.28,0.1,0.195,0.22]) 
K = plt.axes([0.54,0.1,0.195,0.22]) 
L = plt.axes([0.8,0.1,0.195,0.22]) 

# 图A
plt.sca(A)
for i in range(0,6):
    plt.plot(f,mean_osci[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
plt.ylim((-14, 8)) 

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-10, 8, 5)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
A.tick_params(length=2)


plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
A.set_ylabel('Power',labelpad=-0.5)
A.set_xlabel('Frequency (Hz)',labelpad=2)



# A=7*np.ones((387,1))
# plt.plot(f[0:2],A[0:2],color='red',linewidth=0.7)
# plt.plot(f[6:12],A[6:12],color='y',linewidth=0.7)
# plt.plot(f[26:36],A[26:36],color='red',linewidth=0.7)
# plt.plot(f[47:92],A[47:92],color='y',linewidth=0.7)
# plt.plot(f[118:245],A[118:245],color='red',linewidth=0.7)

# plt.plot(f[0]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)
# plt.plot(f[2]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)

# plt.plot(f[6]*np.ones((2,1)),[6.9,7.1],color='y',linewidth=0.7)
# plt.plot(f[12]*np.ones((2,1)),[6.9,7.1],color='y',linewidth=0.7)

# plt.plot(f[26]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)
# plt.plot(f[36]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)

# plt.plot(f[47]*np.ones((2,1)),[6.9,7.1],color='y',linewidth=0.7)
# plt.plot(f[92]*np.ones((2,1)),[6.9,7.1],color='y',linewidth=0.7)

# plt.plot(f[118]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)
# plt.plot(f[245]*np.ones((2,1)),[6.9,7.1],color='red',linewidth=0.7)

# plt.text(0.008, 6.8, "*", size = 5,color = 'black')
# plt.text(f[6]+0.001, 6.8, "*", size = 5,color = 'black')
# plt.text(f[28], 6.8, "*", size = 5,color = 'black')
# plt.text(f[67], 6.8, "*", size = 5,color = 'black')
# plt.text(f[176], 6.8, "*", size = 5,color = 'black')

legend = plt.legend(age_range, 
                    loc=4, ncol=3 , 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)

#图B
plt.sca(B)
plt.plot(f,osci_r,linewidth=1,color='orange')

plt.xlim((0, 0.25))
plt.ylim((-0.3, 0.4))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-0.3, 0.41, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
B.tick_params(length=2)


plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('R',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
B.set_xlabel('Age (year)',labelpad=2)
FDR_line = 0.115*np.ones((257,1))
plt.plot(f2[0:257],FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')
plt.plot(f2[0:257],-FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')

#图C
plt.sca(C)
sns.regplot(x=age, y=SC1, 
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.02, 0.038)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0.02, 0.035, 0.005)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
C.tick_params(length=2)


plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 1 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
C.set_xlabel('Age (year)',labelpad=2)
plt.text(51, 0.038-(2/19)*(0.038-0.02), "r = 0.15  p = 0.007 ", size = rp_size,color = 'black',
         fontname='Arial')



#图D
plt.sca(D)
sns.regplot(x=age, y=SC2,  
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.075, 0.12)) 

my_x_ticks = np.arange(20, 81, 20)
#my_y_ticks = np.arange(0.05, 0.12, 0.02)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(fontname='Arial',fontsize=axis_numsize)
D.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 2 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
D.set_xlabel('Age (year)',labelpad=2)
plt.text(51, 0.12-(2/19)*(0.12-0.075), "r = 0.15  p = 0.008 ", size = rp_size,color = 'black',
         fontname='Arial')


# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec_deconv.mat'
SC_path = 'E:/lifespan/time_diff/original/figure/SC_deconv.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
spec_file = h5py.File(spec_path)
mean_osci = np.squeeze(spec_file['mean_osci'][:])
osci_r = np.squeeze(spec_file['osci_r'][:])
SC_file = h5py.File(SC_path)
SC1 = np.squeeze(SC_file['SC1'][:])
SC2 = np.squeeze(SC_file['SC2'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

#图E
plt.sca(E)

for i in range(0,6):
    plt.plot(f,mean_osci_deconv[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
plt.ylim((-0.14, 0.16)) 

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-0.1, 0.15, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
E.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
E.set_ylabel('Power',labelpad=-0.5)
E.set_xlabel('Frequency (Hz)',labelpad=2)


# B=(0.16-0.3*1/22)*np.ones((387,1))
# plt.plot(f[29:43],B[29:43],color='red',linewidth=1)
# plt.plot(f[60:106],B[60:106],color='y',linewidth=1)



# plt.plot(f[28]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='red',linewidth=0.7)
# plt.plot(f[43]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='red',linewidth=0.7)

# plt.plot(f[60]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='y',linewidth=0.7)
# plt.plot(f[106]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='y',linewidth=0.7)


# plt.text(f[33], (0.16-0.3*1/22)-1/220*0.3, "*", size = 5,color = 'black')
# plt.text(f[82], (0.16-0.3*1/22)-1/220*0.3, "*", size = 5,color = 'black')

legend = plt.legend(age_range, 
                    loc=4, ncol=3, 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)

#图F
plt.sca(F)
plt.plot(f,osci_r,linewidth=1,color='orange')

plt.xlim((0, 0.25))
plt.ylim((-0.3, 0.3))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-0.3, 0.31, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
F.tick_params(length=2)


plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('R',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
F.set_xlabel('Age (year)',labelpad=2)

FDR_line = 0.136*np.ones((257,1))
plt.plot(f2[0:257],FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')
plt.plot(f2[0:257],-FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')

#图G
plt.sca(G)

sns.regplot(x=age, y=SC1,  
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.01, 0.04)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0.01, 0.05, 0.01)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
G.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 1 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
G.set_xlabel('Age (year)',labelpad=2)
plt.text(49, 0.04-(2/19)*(0.04-0.01), "r = 0.22  p < 0.0001 ", size = rp_size,color = 'black',
         fontname='Arial')

#图F
#plt.subplot(1,3,3)
plt.sca(H)

sns.regplot(x=age, y=SC2, 
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.09, 0.155)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0.09, 0.161, 0.02)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
H.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 2 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
H.set_xlabel('Age (year)',labelpad=2)
plt.text(51, 0.155-(2/19)*(0.155-0.09), "r = 0.10  p = 0.063 ", size = rp_size,color = 'black',
         fontname='Arial')

# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec_nonC.mat'
SC_path = 'E:/lifespan/time_diff/original/figure/SC_nonC.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
spec_file = h5py.File(spec_path)
mean_osci = np.squeeze(spec_file['mean_osci'][:])
osci_r = np.squeeze(spec_file['osci_r'][:])
SC_file = h5py.File(SC_path)
SC1 = np.squeeze(SC_file['SC1'][:])
SC2 = np.squeeze(SC_file['SC2'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

#图I
plt.sca(I)

for i in range(0,6):
    plt.plot(f,mean_osci[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=1)

plt.xlim((0, 0.25))
plt.ylim((-400, 230)) 

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-400, 230, 200)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
I.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
I.set_ylabel('Power',labelpad=-0.5)
I.set_xlabel('Frequency (Hz)',labelpad=2)




# B=(0.16-0.3*1/22)*np.ones((387,1))
# plt.plot(f[29:43],B[29:43],color='red',linewidth=1)
# plt.plot(f[60:106],B[60:106],color='y',linewidth=1)



# plt.plot(f[28]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='red',linewidth=0.7)
# plt.plot(f[43]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='red',linewidth=0.7)

# plt.plot(f[60]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='y',linewidth=0.7)
# plt.plot(f[106]*np.ones((2,1)),[(0.16-0.3*1/22)+1/220*0.3,(0.16-0.3*1/22)-1/220*0.3],color='y',linewidth=0.7)


# plt.text(f[33], (0.16-0.3*1/22)-1/220*0.3, "*", size = 5,color = 'black')
# plt.text(f[82], (0.16-0.3*1/22)-1/220*0.3, "*", size = 5,color = 'black')

legend = plt.legend(age_range, 
                    loc=4, ncol=3, 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)

#图J
plt.sca(J)
plt.plot(f,osci_r,linewidth=1,color='orange')

plt.xlim((0, 0.25))
plt.ylim((-0.3, 0.3))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(-0.3, 0.31, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
J.tick_params(length=2)


plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('R',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
J.set_xlabel('Age (year)',labelpad=2)

FDR_line = 0.124*np.ones((257,1))
plt.plot(f2[0:257],FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')
plt.plot(f2[0:257],-FDR_line[0:257],color='red',linewidth=0.7,linestyle='dashed')

#图K
plt.sca(K)

sns.regplot(x=age, y=SC1,  
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.01, 0.02)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0.01, 0.021, 0.005)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
K.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 1 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
K.set_xlabel('Age (year)',labelpad=2)
plt.text(52, 0.02-(2/19)*(0.02-0.01), "r = -0.06  p = 0.26 ", size = rp_size,color = 'black',
         fontname='Arial')

#图L
#plt.subplot(1,3,3)
plt.sca(L)

sns.regplot(x=age, y=SC2, 
            scatter_kws={'s':1.5,'color':'orange'},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0.065, 0.085)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0.07, 0.085, 0.01)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
L.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Spectral centroid 2 (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
L.set_xlabel('Age (year)',labelpad=2)
plt.text(52, 0.085-(2/19)*(0.085-0.065), "r = 0.00  p = 0.98 ", size = rp_size,color = 'black',
         fontname='Arial')


os.chdir('E:/lifespan/time_diff/figure_revision')
plt.savefig('Fig2_RGB.pdf', bbox_inches = 'tight')
plt.show()

"""趋势的系数"""
# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec.mat'
coef_path = 'E:/lifespan/time_diff/original/figure/GS_coef.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
trend_file = h5py.File(spec_path)
trend = np.squeeze(trend_file['mean_trend'][:])
coef_file = h5py.File(coef_path)
GS_a = np.squeeze(coef_file['GS_a'][:])
GS_b = np.squeeze(coef_file['GS_b'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

cmap_path = 'E:/GS Coherence/Colormap/jet.csv'
cmap = pd.read_csv(cmap_path)
cmap2 = np.array(cmap)

#字体大小
axis_numsize = 6
label_size = 7
rp_size = 6


#画图
plt.figure(figsize=(6,6))

A = plt.axes([0.02,0.7,0.26,0.22]) 
B = plt.axes([0.38,0.7,0.26,0.22])
C = plt.axes([0.74,0.7,0.26,0.22]) 

D = plt.axes([0.02,0.4,0.26,0.22])
E = plt.axes([0.38,0.4,0.26,0.22]) 
F = plt.axes([0.74,0.4,0.26,0.22]) 

G = plt.axes([0.02,0.1,0.26,0.22])
H = plt.axes([0.38,0.1,0.26,0.22]) 
I = plt.axes([0.74,0.1,0.26,0.22])

#图A
plt.sca(A)
for i in range(0,6):
    plt.plot(f,trend[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=0.7)
plt.xlim((0, 0.25))
plt.ylim((0, 40))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(0, 41, 10)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
A.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
A.set_ylabel('Power',labelpad=7.5)
A.set_xlabel('Frequency (Hz)',labelpad=2)
legend = plt.legend(age_range, 
                    loc=1, ncol=3, 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)



# 图B
#plt.subplot(1,3,2)
plt.sca(B)
sns.regplot(x=age, y=GS_a,  
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((0, 6.2)) 

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0, 6.2, 1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
B.tick_params(length=2)


plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. a (power law)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
B.set_xlabel('Age (year)',labelpad=2)
B.set_ylabel('Coef. a (power law)',labelpad=7)
plt.text(49, 6.2-(2/19)*(6.2), "r = -0.33  p < 0.0001 ", size = rp_size,color = 'black',
         fontname='Arial')



#图C
plt.sca(C)
sns.regplot(x=age, y=GS_b,  
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((-1, 0)) 

my_x_ticks = np.arange(20, 81, 20)
#my_y_ticks = np.arange(0.05, 0.12, 0.02)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(fontname='Arial',fontsize=axis_numsize)
C.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. b (power law)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
C.set_xlabel('Age (year)',labelpad=2)
plt.text(52, 0-(2/19)*(1), "r = 0.04  p = 0.44 ", size = rp_size,color = 'black',
         fontname='Arial')


# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec_deconv.mat'
coef_path = 'E:/lifespan/time_diff/original/figure/GS_coef_deconv.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
spec_file = h5py.File(spec_path)
trend = np.squeeze(spec_file['mean_trend'][:])
coef_file = h5py.File(coef_path)
GS_a = np.squeeze(coef_file['GS_a'][:])
GS_b = np.squeeze(coef_file['GS_b'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])


#图D
plt.sca(D)
for i in range(0,6):
    plt.plot(f,trend[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=0.7)
plt.xlim((0, 0.25))
plt.ylim((0.08, 0.25))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(0.1, 0.26, 0.05)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
D.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
D.set_xlabel('Frequency (Hz)',labelpad=2)
legend = plt.legend(age_range, 
                    loc=1, ncol=3, 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)




#图E
#plt.subplot(1,3,2)
plt.sca(E)

sns.regplot(x=age, y=GS_a, 
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((-1.2, 2.4))

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(-1, 2.4, 1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
E.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. a (linear)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
E.set_xlabel('Age (year)',labelpad=2)

plt.text(50, 2.4-(2/19)*(1.2+2.4), "r = 0.21  p < 0.0001 ", size = rp_size,color = 'black',
         fontname='Arial')

#图F
#plt.subplot(1,3,3)
plt.sca(F)

sns.regplot(x=age, y=GS_b,
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((-0.05, 0.35))

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(-0.05, 0.36, 0.1)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
F.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. b (linear)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
F.set_xlabel('Age (year)',labelpad=2)
plt.text(49, 0.35-(2/19)*(0.36+0.05), "r = -0.23  p < 0.0001 ", size = rp_size,color = 'black',
         fontname='Arial')


# 读取数据
age_path = 'E:/lifespan/time_diff/original/figure/age.mat'
spec_path = 'E:/lifespan/time_diff/original/figure/mean_spec_nonC.mat'
coef_path = 'E:/lifespan/time_diff/original/figure/GS_coef_nonC.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'

age_file = h5py.File(age_path)
age = np.squeeze(age_file['age'][:])
spec_file = h5py.File(spec_path)
trend = np.squeeze(spec_file['mean_trend'][:])
coef_file = h5py.File(coef_path)
GS_a = np.squeeze(coef_file['GS_a'][:])
GS_b = np.squeeze(coef_file['GS_b'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

#图G
plt.sca(G)
for i in range(0,6):
    plt.plot(f,trend[i,:],color=cmap2[15+i*40,:],label=age_range[i],linewidth=0.7)
plt.xlim((0, 0.25))
plt.ylim((0, 2000))

my_x_ticks = np.arange(0, 0.26, 0.05)
my_y_ticks = np.arange(0, 2001, 1000)
y_ticks_label = ['0','1,000','2,000']
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,y_ticks_label,fontname='Arial',fontsize=axis_numsize)
G.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
G.set_xlabel('Frequency (Hz)',labelpad=2)
legend = plt.legend(age_range, 
                    loc=1, ncol=3, 
                    title='Age range (year)',
                    prop={"family" : "Arial",'size':5},
                    fontsize=5,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 6)




#图H
plt.sca(H)

sns.regplot(x=age, y=GS_a, 
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((-0.2, 8))

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(0, 8.1, 2)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
H.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. a (power law)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
H.set_xlabel('Age (year)',labelpad=2)

plt.text(50, 8-(2/19)*(8), "r = -0.12  p = 0.03 ", size = rp_size,color = 'black',
         fontname='Arial')

#图I
plt.sca(I)

sns.regplot(x=age, y=GS_b,
            scatter_kws={'s':1.5,"color": "orange"},line_kws={"color": "black",'linewidth' : 1.2},ci=0)

plt.xlim((15, 82))
plt.ylim((-2, -0.5))

my_x_ticks = np.arange(20, 81, 20)
my_y_ticks = np.arange(-2, -0.49, 0.5)
plt.xticks(my_x_ticks,fontname='Arial',fontsize=axis_numsize)
plt.yticks(my_y_ticks,fontname='Arial',fontsize=axis_numsize)
I.tick_params(length=2)

plt.xlabel('Age (year)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
plt.ylabel('Coef. b (power law)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})
I.set_xlabel('Age (year)',labelpad=2)
plt.text(50, -0.5-(2/19)*(-0.5+2), "r = -0.11  p = 0.04 ", size = rp_size,color = 'black',
         fontname='Arial')


os.chdir('E:/lifespan/time_diff/figure_revision')
plt.savefig('Fig3_RGB.pdf', bbox_inches = 'tight')
plt.show()

"""图3 头动"""
# 读取数据
FD_path = 'E:/lifespan/time_diff/original/figure/meanFD_spec.mat'
f_path = 'E:/lifespan/time_diff/original/figure/f.mat'
FD_file = h5py.File(FD_path)
R_FD = np.squeeze(FD_file['mean_spec'][:])
f_file = h5py.File(f_path)
f = np.squeeze(f_file['f2'][:])

f2_path = 'E:/lifespan/time_diff/original/figure/f_original.mat'
f2_file = h5py.File(f2_path)
f2 = np.squeeze(f2_file['f'][:])


#年龄范围label
age_range = ['-25','25-30','30-35','35-40','40-45','45-50',
             '50-55','55-60','60-65','65-70','70-']

axis_numsize = 8
label_size = 9
rp_size = 6

plt.figure(figsize=(3.54330708661,3.54330708661*2/3))

ax=plt.subplot(1,1,1)
for i in range(0,11):
    plt.plot(f,R_FD[i,:],color=[i*0.09,1-i*0.1,i*0.01,1],label=age_range[i],linewidth=0.7)
plt.xlim((0, 0.25))
plt.ylim((0, 1.9))


plt.xticks(fontname='Arial',fontsize=axis_numsize)
plt.yticks(fontname='Arial',fontsize=axis_numsize)
ax.tick_params(length=2)

plt.xlabel('Frequency (Hz)',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

plt.ylabel('Power',
{'family' : 'Arial',
'weight' : 'normal',
'size'   : label_size})

B=(1.9-1.9*1/22)*np.ones((387,1))
plt.plot(f[0:19],B[0:19],color='r',linewidth=1)



plt.plot(f[0]*np.ones((2,1)),[(1.9-1.9*1/22)+1/220*1.9,(1.9-1.9*1/22)-1/220*1.9],color='r',linewidth=0.7)
plt.plot(f[19]*np.ones((2,1)),[(1.9-1.9*1/22)+1/220*1.9,(1.9-1.9*1/22)-1/220*1.9],color='r',linewidth=0.7)




plt.text(f[8], (1.9-1.9*1/22)-1/220*1.9, "*", size = 5,color = 'black')

legend = plt.legend(age_range, 
                    loc=1, ncol=3, 
                    title='Age Range',
                    prop={"family" : "Arial",'size':6},
                    fontsize=6,borderpad=0.3,handlelength=1,handletextpad=0.5,
                    columnspacing=0.7)
legend.get_title().set_fontsize(fontsize = 7)

os.chdir('E:/lifespan/time_diff/figure_new')
plt.savefig('Fig4_RGB.pdf', bbox_inches = 'tight')
plt.show()


