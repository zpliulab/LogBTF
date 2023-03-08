#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 19:41:16 2023

@author: lilingyu
"""


## 2023.2.11 The code is based on ExamTest1130.py established on 2022.12.3 by Lingyu Li
## 2023.2.14 Using the Melanoma Simulated Data form. 
##           /Users/lilingyu/E/PhD/R/Landscape/R/RegNetPlot.R


## Read Expression data
pathname = "/home/lly/"  

## Mac
# pathname = "/Users/lilingyu/E/PhD/"



import os
os.chdir(str(pathname)+'Python/SERGIO-master')



## Set random seed
import random
random.seed(1) 
import numpy as np
np.random.seed(1)


import numpy as np
import pandas as pd
from SERGIO.sergio import sergio


################################################################################
## Simulate Clean Data _ Steady-State Simulation


## simulation:
# number_genes：GRN 中存在的基因总数 
# number_bins：要模拟的不同细胞类型的总数 
# number_sc：要模拟的每种细胞类型的细胞总数
# noise_params：建议使⽤⼩值 (<0.5)

# sim = sergio(number_genes=100, number_bins = 9, number_sc = 300, noise_params = 1, 
#              decays=0.8, sampling_state=15, noise_type='dpd')


# =============================================================================
# ################################################################################
# ## Data -- genes 100; cell-type 2; cell-numb 10;
# number_genes=100
# number_bins = 2
# number_sc = 10
# =============================================================================


# =============================================================================
# ################################################################################
# ## Data -- genes 100; cell-type 2; cell-numb 25;
# number_genes=100
# number_bins = 2
# number_sc = 25
# =============================================================================


################################################################################
## Data -- genes 100; cell-type 2; cell-numb 50;
number_genes= 20
number_bins = 1
number_sc = 10

shape = 6.5
percentile = 52




## mycode
sim = sergio(number_genes, number_bins, number_sc, noise_params = 1, 
             decays=0.8, sampling_state=15, noise_type='dpd')


sim.build_graph(input_file_taregts ='data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/Interaction_cID_4_20node.txt', 
                input_file_regs='data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised/Regs_cID_4_'+str(number_bins)+'type_20node.txt', 
                shared_coop_state=2)


## 运⾏稳态模拟
sim.simulate()

## 在 steady_state 模拟调⽤ getExpressions ⽅法后获得⼲净的模拟表达式矩阵
## 返回⼀个 3d numpy 数组（#cell_types * #genes * #cells_per_type）
expr = sim.getExpressions()

## 转换为⼤⼩ 为 (#genes * #cells) 的⼆维矩阵
expr_clean = np.concatenate(expr, axis = 1)

## Save the result
os.chdir(str(pathname)+'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised')
np.savetxt('exper_clean_'+str(number_genes)+'_'+str(number_bins)+'_'+str(number_sc)+'.csv',
            expr_clean, delimiter = ',')



################################################################################
## Add Technical Noise _ Steady-State Simulations

"""
Add outlier genes
"""
expr_O = sim.outlier_effect(expr, outlier_prob = 0.01, mean = 0.8, scale = 1)

"""
Add Library Size Effect
"""
libFactor, expr_O_L = sim.lib_size_effect(expr_O, mean = 4.6, scale = 0.4)

"""
Add Dropouts
"""
# binary_ind = sim.dropout_indicator(expr_O_L, shape = 6.5, percentile = 82)
binary_ind = sim.dropout_indicator(expr_O_L, shape = shape, percentile = percentile)
## the number of 0
print(np.sum(binary_ind == 0))
expr_O_L_D = np.multiply(binary_ind, expr_O_L)

"""
Convert to UMI count
"""
count_matrix = sim.convert_to_UMIcounts(expr_O_L_D)

"""
Make a 2d gene expression matrix
"""
count_matrix = np.concatenate(count_matrix, axis = 1)
print(np.sum(count_matrix == 0))

## Save the result
os.chdir(str(pathname)+'Python/SERGIO-master/data_sets/De-noised_100G_9T_300cPerT_4_DS1LLYBioRevised')
np.savetxt('count_matrix'+str(number_genes)+'_'+str(number_bins)+'_'+str(number_sc)+'.csv',
            count_matrix, delimiter = ',')


################################################################################









