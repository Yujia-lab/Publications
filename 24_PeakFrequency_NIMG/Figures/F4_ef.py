



import scipy.io as scio
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pandas as pd
import h5py
import matplotlib.ticker as ticker


in_path = 'For seaborn2'
out_path = 'edge_frequency'
subj_num = 179
# 导入数据集
img_name = ['mean_scatter_band1','mean_scatter_band2','mean_scatter_band3','mean_scatter_band4',
            'mean_scatter_band5','mean_scatter_band6','mean_scatter_band7','mean_scatter_band8']


# 绘图

tick_size = 16
label_size = 20

path = '/Users/yujia/Data/Data/HCP/Figures/'+in_path+'/figure4_scatter.mat'
PV_file = scio.loadmat(path)
x1 = np.squeeze(PV_file['x2'])
y1 = np.squeeze(PV_file['y1'])




for i in range(0,5):
        
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
    os.chdir('/Users/yujia/Data/Data/HCP/Figures/'+out_path)
    plt.savefig(img_name[i],dpi=300)
    plt.show()