#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:03:22 2023

@author: lilingyu
"""

## Set random seed
import random
random.seed(1) 
import numpy as np
np.random.seed(1)


################################################################################
################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import statsmodels.api as sm
from scipy.stats import f
import scipy.stats as stat
import math
from sklearn import metrics
from sklearn.metrics import auc 
import scipy as sp
import os

## returns the list minimum index
def indexofMin(arr):
	minindex = 0
	currentindex = 1
	while currentindex < len(arr):
		if arr[currentindex] < arr[minindex]:
			minindex = currentindex
		currentindex += 1
	return minindex

# sorting 
def sor_g(gene_1,gene_2):
    com = []
    for j in range(len(gene_1)):
        tup = (gene_1[j],gene_2[j])
        com.append(tup)
    sort_1 = sorted(com, key=lambda x:x[0])#Sort by gene_1
    sort_2 = sorted(com, key=lambda x:x[1])#Sort by gene_2
    a1 = []
    b1 = []
    for k in range(len(gene_1)):
        a1.append(sort_1[k][0])
        b1.append(sort_1[k][1])   
    return a1,b1
#Lagged regression
def granger(gene_1,gene_2):
	com = []
	for j in range(len(samples1)):
		tup = (gene_1[j],gene_2[j])
		com.append(tup)
	sort_1 = sorted(com, key=lambda x:x[0])
	# granger test based on gene sorting, 
	a1 = []
	b1 = []
	for k in range(len(samples1)):
		a1.append(sort_1[k][0])
		b1.append(sort_1[k][1])
	a1_cha = []
	b1_cha = []
	for k in range(len(samples1)):
		a1_cha.append(300*a1[k]/90001)
		#a1_cha.append(3*a1[k]/90001)
		b1_cha.append(300*b1[k]/90001)
		#b1_cha.append(3*b1[k]/90001)
	# Lag phase 1
	a1_t_1 = np.array(a1_cha[0:len(samples1)-1])
	a1_t1_1= np.array(a1_cha[1:len(samples1)])
	b1_t_1 = np.array(b1_cha[0:len(samples1)-1])
	b1_t1_1 = np.array(b1_cha[1:len(samples1)])   
	b1_t1_t2_1= sm.add_constant(b1_t1_1)
	a1_t1_t2_b1_1 = np.vstack((a1_t1_1,b1_t1_1)).T
	a1_t1_t2_b1_1 = sm.add_constant(a1_t1_t2_b1_1)
	 #F-test
	model = LinearRegression()
	model.fit(a1_t1_t2_b1_1,b1_t_1)
	b1_rssu_1 = np.sum((model.predict(a1_t1_t2_b1_1)-b1_t_1)**2)
	model.fit(b1_t1_t2_1,b1_t_1)
	b1_rssr_1 = np.sum((model.predict(b1_t1_t2_1)-b1_t_1)**2)
	f1_1 = ((b1_rssr_1-b1_rssu_1)/1)/(b1_rssu_1/(len(samples1)-3))
	p1_1 = f.sf(f1_1,1,len(samples1)-3) 
	 #Lag phase 2
	a1_t = np.array(a1_cha[0:len(samples1)-2])
	a1_t1 = np.array(a1_cha[1:len(samples1)-1])
	a1_t2 = np.array(a1_cha[2:len(samples1)])
	b1_t = np.array(b1_cha[0:len(samples1)-2])
	b1_t1 = np.array(b1_cha[1:len(samples1)-1])
	b1_t2 = np.array(b1_cha[2:len(samples1)])		
	a1_t1_t2 =  np.vstack((a1_t1,a1_t2)).T
	a1_t1_t2 = sm.add_constant(a1_t1_t2)		
	b1_t1_t2 = np.vstack((b1_t1,b1_t2)).T
	b1_t1_t2 = sm.add_constant(b1_t1_t2)		
	a1_t1_t2_b1 = np.vstack((a1_t1,a1_t2,b1_t1,b1_t2)).T
	a1_t1_t2_b1 = sm.add_constant(a1_t1_t2_b1)
	#F-test
	model.fit(a1_t1_t2_b1,b1_t)
	b1_rssu_2 = np.sum((model.predict(a1_t1_t2_b1)-b1_t)**2)
	model.fit(b1_t1_t2,b1_t)
	b1_rssr_2 = np.sum((model.predict(b1_t1_t2)-b1_t)**2) 
	f1_2 = ((b1_rssr_2-b1_rssu_2)/2)/(b1_rssu_2/(len(samples1)-2-2-1))
	p1_2= f.sf(f1_2,2,len(samples1)-2-2-1)
	#Lag phase 3
	a1_t_3 = np.array(a1_cha[0:len(samples1)-3])
   
	a1_t1_3= np.array(a1_cha[1:len(samples1)-2])
	a1_t1_t3_3=np.array(a1_cha[2:len(samples1)-1])
	a1_t2_3= np.array(a1_cha[3:len(samples1)])
	
	b1_t_3= np.array(b1_cha[0:len(samples1)-3])
	
	b1_t1_3= np.array(b1_cha[1:len(samples1)-2])
	b1_t1_t3_3=np.array(b1_cha[2:len(samples1)-1])
	b1_t2_3= np.array(b1_cha[3:len(samples1)])  
	  
	a1_t1_t2_3=  np.vstack((a1_t1_3,a1_t2_3,a1_t1_t3_3)).T
	a1_t1_t2_3= sm.add_constant(a1_t1_t2_3) 
	   
	b1_t1_t2_3= np.vstack((b1_t1_3,b1_t2_3,b1_t1_t3_3)).T
	b1_t1_t2_3= sm.add_constant(b1_t1_t2_3)		
	a1_t1_t2_b1_3 = np.vstack((a1_t1_3,a1_t2_3,a1_t1_t3_3,b1_t1_3,b1_t2_3,b1_t1_t3_3)).T
	a1_t1_t2_b1_3 = sm.add_constant(a1_t1_t2_b1_3)
	#F-test
	model = LinearRegression()
	model.fit(a1_t1_t2_b1_3,b1_t_3)
	b1_rssu_3 = np.sum((model.predict(a1_t1_t2_b1_3)-b1_t_3)**2)
	model.fit(b1_t1_t2_3,b1_t_3)
	b1_rssr_3 = np.sum((model.predict(b1_t1_t2_3)-b1_t_3)**2)
	f1_3 = ((b1_rssr_3-b1_rssu_3)/3)/(b1_rssu_3/(len(samples1)-7))
	p1_3 = f.sf(f1_3,3,len(samples1)-7)
	p=[p1_1,p1_2,p1_3]
	p_1=min(p)
	return p_1

def LassoRegression(degree, alpha):
	return Pipeline([
		("poly", PolynomialFeatures(degree=degree)),
		("std_scaler", StandardScaler()),
		("lasso_reg", Lasso(alpha=alpha))])

#Lasso regression:[degree:Highest power exponent; regularization] 
#[parameters:different types of data for different training values]
def lasso_regress(gene_1,gene_2,degree,alpha):
	le_g1=len(gene_1)#长度
	co=max(gene_1)-min(gene_1)#极差
	st=co/le_g1 #步长 
	t=np.linspace(min(gene_1), max(gene_1), num=le_g1)
	#t= np.arange(min(gene_1), max(gene_1), st) #根据极差生成即将映射的y值   
	a1,b1=sor_g(gene_1,gene_2)#排序后的g1和g2
	a1_1=np.array(a1)#np下的g1
	a1_1= a1_1.reshape(-1, 1)
	b1_1=np.array(b1)
	b1_1=b1_1.reshape(-1,1)
	lasso1_reg = LassoRegression(degree,alpha) 
	lasso1_reg.fit(a1_1, t)
	gene2=lasso1_reg.predict(b1_1)
	gene1=lasso1_reg.predict(a1_1)
	return gene1,gene2



################################################################################
################################################################################

## Read Expression data
pathname = "/home/lly/"
Boolean = "Data"
## Mac
# pathname = "/Users/lilingyu/E/PhD/"
# Boolean = "Rshell/Data"
## Mac
# os.chdir('/home/lly/R/Boolean/Data/DREAM3result_GNIPLR/RData/scData')
os.chdir(str(pathname)+'R/Boolean/'+str(Boolean)+'/DREAM3result_GNIPLR/RData/')


## list
fileEcoli = ([10,1,4,45], [50,1,12,576], [100,1,29,1708], 
             [10,2,3,48], [50,2,6,585], [100,2,11,1996])

fileYeast = ([10,1,1,53], [50,1,10,415], [100,1,19,1497], 
             [10,2,4,59], [50,2,20,557], [100,2,4,1700],
             [10,3,2,40], [50,3,16,417], [100,3,42,1718])
## array
fileEcoli0 = np.array(fileEcoli)
fileYeast0 = np.array(fileYeast)
## test
fileEcoli0[1,1]


## Set parameters by user
################################################################################
################################################################################
## "Ecoli"   "Yeast"   fileEcoli0    fileYeast
dataname = "Ecoli"
file = fileEcoli0
datafile = 'scData'

np.shape(file)[0]
################################################################################
################################################################################


# i = 0
for i in range(np.shape(file)[0]):

    genenum = file[i,0]
    filenum = file[i,1]
    k = file[i,2]
    geneselect = file[i,3]
    
    
    
    ax_nl=0.04
    n_d4 =str(datafile)+'/'+str(datafile)+'_'+str(dataname)+str(filenum)+'_Node'+str(genenum)+'.txt'
    train_s=pd.read_csv(str(n_d4),sep='\t')
    train_s.T.to_csv(str(datafile)+'/'+str(dataname)+'/y1_10_d2_t.txt', sep='\t')
    
    
    
    ## ReadGround-truth GRN  
    os.chdir(str(pathname)+'R/Boolean/Data/DREAM/DREAM3 in silico challenge/')
    # os.chdir('/Users/lilingyu/E/PhD/R/Boolean/Data/DREAM/DREAM3 in silico challenge/')
    ## Link file
    fw=open('Size'+str(genenum)+'/DREAM3 gold standards/'
            +'DREAM3GoldStandard_InSilicoSize'
            +str(genenum)+'_'+str(dataname)+str(filenum)+'.txt')
    
    rels=[]
    for p in fw:
    	t1=p.split()
    	rels.append(t1)
    fw.close()  
    fla=0
    cancer1 = {}
    samples1 = []
    
    
    ## load new data 
    os.chdir(str(pathname)+'R/Boolean/'+str(Boolean)+'/DREAM3result_GNIPLR/RData/')
    file1 = open(str(datafile)+'/'+str(dataname)+'/y1_10_d2_t.txt')
    for p in file1:
    	fla+=1
    	t1 = p.split()
    	n1=len(t1)
    	samples1=t1[1:]
    	cancer1[t1[0]]=[t1[i] for i in range(1,n1)]
    if '0' in cancer1:
    	del cancer1['0']
    if 'wt' in cancer1:
    	del cancer1['wt']
    if 'strain' in cancer1:
    	del cancer1['strain']
    if 'Time' in cancer1:
    	del cancer1['Time']
    file1.close()
    le=len(cancer1)
    for sub in cancer1:
    	cancer1[sub]=[float(x) for x in cancer1[sub]] 
        
        
    ## calculate and return the prediction
    fw=open(str(datafile)+'/'+str(dataname)+'/e1_100_d2_p1.txt','w')
    for i in range(len(rels)):
    	entry_11=rels[i][0]
    	entry_21=rels[i][1]
    	gene_1,gene_2=lasso_regress(cancer1[entry_11],cancer1[entry_21],30,ax_nl) 
    	p=granger(gene_1,gene_2)
    	fw.write(str(entry_11)+'\t'+str(entry_21)+'\t'+str(p)+'\t'+str(rels[i][2])+'\n')
    fw.close()   
    
    
    ## Save results
    relust1=pd.read_csv(str(datafile)+'/'+str(dataname)+'/e1_100_d2_p1.txt',sep='\t',names = range(4))
    # bbb=relust1.ix[:,[0,1,2,3]]    ## 2023 LLY Add
    bbb=relust1.loc[:,[0,1,2,3]]
    bbb.to_csv(str(n_d4[:-4])+'_lasso.txt', sep='\t',index=False,header=False)
    
    
    
    
    
    
