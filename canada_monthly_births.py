# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:13:44 2021

@author: Farzaneh
"""
del()
import os
import sys
assert sys.version_info >= (3, 5)
# Scikit-Learn â‰¥0.20 is required
import sklearn
assert sklearn.__version__ >= "0.20"    
from scipy.io import loadmat
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import hankel
import math
from scipy.fftpack import fft
from scipy.fftpack import ifft
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
import datetime as dt
import time
# Where to save the figures
PROJECT_ROOT_DIR = "."
CHAPTER_ID = "Trend"
IMAGES_PATH = os.path.join(PROJECT_ROOT_DIR, "images", CHAPTER_ID)
os.makedirs(IMAGES_PATH, exist_ok=True)

def save_fig(fig_id, tight_layout=True, fig_extension="png", resolution=300):
    path = os.path.join(IMAGES_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)
    
def antidiag(hank_red,p1,p2,lent):
    vector = np.zeros(lent)
    for k in range (2):
        if k==0:
            for j in range (p2):
                 vector[j] = hank_red[0,j]
        if k==1:
            for i in range (p1):
                vector[p2+i-1] = hank_red[i,p2-1] 
    return (vector)
#mydata = loadmat('C:\\Users\\Farzaneh\\Documents\\MATLAB\\DATA\\test_linear.mat')
mydata =pd.read_excel('C:\\Users\\Farzaneh\\Documents\\DATA\\Canada_birth.xlsx')
#data=mydata[0:,1]
#date=mydata[0:,0]
plt.figure(figsize=(6,4))

mydata["Data"].plot()

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import hankel
import math
from scipy.fftpack import fft
from scipy.fftpack import ifft
from sklearn.decomposition import PCA
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
date=mydata["Year"]
data=mydata["Data"]
win=85
m=340
rank=3
nwin=math.floor(m/win)
mid1=math.floor((win+1)/2)
mid2=math.floor((win+1)/2)
vec=np.zeros(win)
vector=np.zeros(win*nwin)
#window=math.floor(data.shape/win)
for j in range(nwin):
    for i in range(win):
        vec[i]=data[i+j*win]
#    verts=data[i*win:win+i*win,0]
    khank=hankel(vec)
    kkhankr=khank[0:mid1,0:mid2]
    U, s, Vt = np.linalg.svd(kkhankr)
    S = np.zeros(U.shape)
    S[:rank, :rank] = np.diag(s[:rank])
    hank_red = U.dot(S).dot(Vt)#(U@S@Vt) 
    vectorr= antidiag(hank_red,mid2,mid1,win)
    for i in range(win):
        vector[i+j*win]=vectorr[i]
   
plt.figure(figsize=(10,3))
xfmt = md.DateFormatter('%Y-%m')
ax=plt.gca()
plt.xticks( rotation=25 )
ax.xaxis.set_major_formatter(xfmt)
plt.plot(date,data ,'b',linewidth=2)
plt.plot(date,vector, 'r',linewidth=2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=5)
ax.tick_params(which='minor', length=5)
plt.xlabel('Date')
plt.ylabel('Number of live birth per month')
save_fig("Canada_birth")
plt.show()
