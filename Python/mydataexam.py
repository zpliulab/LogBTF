# -*- coding: utf-8 -*-
"""
Created on Sat May 08 10:13:22 2021

@author: Utente
"""

"""
Created on Thu Oct 22 12:03:57 2020

@author: Utente
"""

## 2022.4.26
## Python 的 for循环，需要加：，并且需要严格的tab键
## 生成csv，直接打开没分隔符，只能用read.table，就没有问题
## l=0, 对应的gene_name这一列，所以l=1:2**9即可,，但是需要设置2**9+1
## 看一下文件夹下有多少csv文件 ls -l|grep "^-"| wc -l


import boolean2 as b2
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab
import pandas as pd


#p = 9
#l = 2**9-12
#l = 1

for l in range(1, (2**9+1)):
	list=[]
	while True:
	    try:
	        #insert name and extension file excel with binary values
	        #Table=raw_input("Please enter the excel file name: ")
	                
	        #get a dataframe from the excel sheet
	        df=pd.read_excel("dataAllList.xlsx")
	        
	        try:
	            #indicate the cell in which to look for attractors
	#            cellname="var_2"
	#		  cellname = df.columns[1]
	            
	            #tranforms the values of the columns to be processed from "1" and"0"
	            #to "True" and "False"
	            df['TF']=df[df.columns[l]].apply(lambda x : 'True' if x == 1 else 'False')
	        except IOError:
	            print("The cell name is incorrect")
	            
	    except IOError:
	         print("The file name is incorrect")
	         continue
	    break
	
	
	#create a new dataframe with only two columns,the name of the genes and 
	#their Boolean value "True" or "False"
	mod_df=df[['gene_name','TF']]
	
	#create a list with values in columns of the dataframe "mod_df"
	#interspersed with the symbol "="
	for i in mod_df.index:
	    list.append(mod_df["gene_name"][i] + " = " + mod_df["TF"][i]+"\n")
	
	#tranformation of the list format into a string required for processing
	model_definition1='  '.join(list)
	
	
	print(model_definition1)
	    
	
	
	# Update rule  
	model_definition2= '''
	x1 *= x1 or x2
	x2 *= (not x1) and (not x4) and (not x5) and (not x9)
	x3 *= (not x2) and (not x5) and (not x9)
	x4 *= (not x2) and x3
	x5 *= (not x2) and x3 and x5
	x6 *= x9
	x7 *= ((not x5) and (not x9)) or x6
	x8 *= (not x7) or (x7 and x8 and (x6 or x5 or x9))
	x9 *= (not x6) and (not x7)
	
	'''
	#union of the two defined strings
	model_definition = model_definition1 + model_definition2
	
	#Refers to the text containing in the model definition and the model update
	model = b2.Model(text=model_definition, mode='sync')
	model.initialize()
	#Number of iteration
	model.iterate(steps=1)
	
	for node in model.data:
	    # print node, model.data[node]
	     print(node, model.data[node])
	#Cheking for fixed states
	model.report_cycles()
	
	#Save the result
	model.save_states("Atta/attrdataexam_"+ str(l)+".csv")
	print(l)


