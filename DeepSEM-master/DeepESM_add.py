#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 14:21:56 2023

@author: lilingyu
"""


################################################################################
################################################################################
## my file -- made by LLY    --  scData
 python main.py --task celltype_GRN --data_file DREAM3_data/scData/scData_Ecoli1_Node10.csv --net_file DREAM3_data/label/Size10/label_Ecoli1_Node10.csv --setting new --alpha 0.1 --beta 0.01 --n_epochs 10  --save_name out

## my file -- made by LLY    --  bulkData
 python main.py --task celltype_GRN --data_file DREAM3_data/bulkData/bulkData_Ecoli1_Node10.csv --net_file DREAM3_data/label/Size10/label_Ecoli1_Node10.csv --setting new --alpha 0.1 --beta 0.01 --n_epochs 10  --save_name out
 

################################################################################
################################################################################
## Calculate EPR values
import pandas as pd
import os
os.chdir('/Users/lilingyu/E/PhD/Paper/Paper7/Bioinformatics/AddMethod/ncs曾坚阳/DeepSEM-master/')

output = pd.read_csv('out/GRN_inference_result.tsv',sep='\t')

## delete the rows with 'nan'
## reference https://www.cnpython.com/qa/1407463
# output = output[complete.cases(output), ]
## reference https://blog.csdn.net/weixin_39450145/article/details/120795055
output = output.dropna(axis=0, how='any')

output['EdgeWeight'] = abs(output['EdgeWeight'])
output = output.sort_values('EdgeWeight',ascending=False)

label = pd.read_csv('DREAM3_data/label/Size10/label_Ecoli1_Node10.csv')
TFs = set(label['Gene1'])
Genes = set(label['Gene1'])| set(label['Gene2'])
output = output[output['TF'].apply(lambda x: x in TFs)]
output = output[output['Target'].apply(lambda x: x in Genes)]
label_set = set(label['Gene1']+'|'+label['Gene2'])
output= output.iloc[:len(label_set)]
len(set(output['TF']+'|' +output['Target']) & label_set) / (len(label_set)**2/(len(TFs)*len(Genes)-len(TFs)))




################################################################################
################################################################################
## Calculate AUPR ratio values

from sklearn.metrics import average_precision_score
import numpy as np
import pandas as pd

output = pd.read_csv('out/GRN_inference_result.tsv',sep='\t')
output['EdgeWeight'] = abs(output['EdgeWeight'])
output = output.sort_values('EdgeWeight',ascending=False)
label = pd.read_csv('DREAM3_data/label/Size10/label_Ecoli1_Node10.csv')
TFs = set(label['Gene1'])
Genes = set(label['Gene1'])| set(label['Gene2'])
output = output[output['TF'].apply(lambda x: x in TFs)]
output = output[output['Target'].apply(lambda x: x in Genes)]
label_set = set(label['Gene1']+label['Gene2'])
preds,labels,randoms = [] ,[],[]
res_d = {}
l = []
p= []
for item in (output.to_dict('records')):
        res_d[item['TF']+item['Target']] = item['EdgeWeight']
for item in (set(label['Gene1'])):
        for item2 in  set(label['Gene1'])| set(label['Gene2']):
            if item+item2 in label_set:
                l.append(1)
            else:
                l.append(0)
            if item+ item2 in res_d:
                p.append(res_d[item+item2])
            else:
                p.append(-1)
# average_precision_score(l,p)/np.mean(l)

## output
average_precision_score(l,p)



################################################################################
################################################################################
## Ensemble DeepSEM result

res = []
for i in range(10):
    res.append(pd.read_csv('../../scGRN/Upload/GRN_inference_benchmark/cross_validation/500_STRING_hESC/rep_i.csv',sep='\t'))
res = pd.concat(res)
res['EdgeWeight'] = abs(res['EdgeWeight'])
res.groupby(['Gene1','Gene2']).mean()

