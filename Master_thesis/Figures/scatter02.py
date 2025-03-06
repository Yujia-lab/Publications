#前测容量相关
from matplotlib import pyplot as plt 
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as sci #计算样本的偏度和峰度

df=pd.read_excel("颜色记忆容量前测_相关.xlsx",sheet_name="相关")
x=df["eprime"]
y=df["phone"]

def pearson(x,y):
    r,p=sci.pearsonr(x,y)
    

df.corr()
sns.jointplot(x,y, data = df,kind = 'reg')
plt.annotate(pearson,template='r:{val:.2f,p:{p:.3f}')
plt.show()