

import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py
import matplotlib.ticker as ticker


in_path = '\\For seaborn2'
out_path = '\\Figure 4 shuffled'
subj_num = 181
# 导入数据集
img_name = ['mean_scatter_band1','mean_scatter_band2','mean_scatter_band3','mean_scatter_band4',
            'mean_scatter_band5','mean_scatter_band6','mean_scatter_band7','mean_scatter_band8']


# 绘图

tick_size = 16
label_size = 20

path = 'E:\\Data\\HCP\\Figures' + in_path + '\\figure4_scatter.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x2'])
y1 = np.squeeze(PV_file['y1'])




for i in range(0,4):
        
    x = np.squeeze(x1[:,i])
    
    ax = sns.regplot(x=x,y=y1)
   
    plt.yticks(fontname='Arial',fontsize=tick_size)
    plt.xticks(fontname='Arial',fontsize=tick_size)
    plt.xlabel("PF",fontname='Arial',fontsize=label_size,labelpad=1)
    plt.ylabel("ACW",fontname='Arial',fontsize=label_size,labelpad=1)
    sns.despine(top=False,right=False)
    ax.figure.set_size_inches(5.5,5)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    os.chdir('E:\\Data\\HCP\\Figures'+out_path)
    plt.savefig(img_name[i],dpi=300)
    plt.show()


    
path = 'E:\\Data\\HCP\\Figures'+in_path+'\\figure4_bar.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x1'])

band1 = ['0.01-0.03 Hz']
band2 = ['0.03-0.05 Hz']
band3 = ['0.05-0.07 Hz']
band4 = ['0.07-0.09 Hz']    
Frequency = band1+band2+band3+band4
Data_mean = pd.DataFrame(data={'CRC': x1, 'Frequency': Frequency})


ax = sns.barplot(x="Frequency", y="CRC", data=Data_mean)
ax.figure.set_size_inches(7,5) 
plt.ylim(-1,0.6)
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Band",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("R",fontname='Arial',fontsize=label_size,labelpad=1)
plt.axhline(y=0)
os.chdir('E:\\Data\\HCP\\Figures'+out_path)
plt.savefig('mean_bar',dpi=300)
plt.show()

'''acf function'''
path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure4_acf.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['ACF'])
x2 = np.squeeze(PV_file['tp'])

ax = sns.lineplot(x=x2, y=x1,color='#FF9600')
ax.figure.set_size_inches(4,5) 

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_position(('data',0))  
ax.spines['bottom'].set_position(('data',0))  
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel(' ',fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel(' ',fontname='Arial',fontsize=label_size,labelpad=1)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/Figure 4')
plt.savefig('acf_30',dpi=300)
plt.show()

path = '/Users/yujia/Data/Data/HCP/Figures/For seaborn/figure4_acf2.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['acf'])
x2 = np.squeeze(PV_file['tp'])

ax = sns.lineplot(x=x2, y=x1,color='#FF5800')
ax.figure.set_size_inches(10,5) 

ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_position(('data',0))  
ax.spines['bottom'].set_position(('data',0))  
plt.xticks(fontname='Arial',fontsize=tick_size)
plt.yticks(fontname='Arial',fontsize=tick_size)
plt.xlabel("Lags (second)",fontname='Arial',fontsize=label_size,labelpad=1)
plt.ylabel("Autocorrelation",fontname='Arial',fontsize=label_size,labelpad=1)


os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
plt.savefig('acf_1200',dpi=300)
plt.show()
