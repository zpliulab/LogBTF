# [LogBTF (Embedding logistic regression into Boolean threshold function to reconstruct Boolean threshold network model)](https://github.com/zpliulab/logBTF)

![Screenshot](Data/framework.jpg)

In this work, an **embedded Boolean threshold network model by aggregating logistic regression with Boolean threshold function (LogBTF)** for inferring gene regulatory networks from single-cell gene expression data was proposed. 


## LogBTFs
<!--START_SECTION:news-->
* **LogBTF**: A **embedded Boolean threshold network model (LogBTF)** is proposed to infer **gene regulatory networks (GRNs)**. 
* In the comparison study, we also proved the proposed **LogBTF method** results in better inference performance than one regression-based method **SINCERITIES**, two miscellaneous methods: **GENIE3** and **TIGRESS**, two correlation-based methods: **ARACNE** and **CLR**.
* If you have any questions about **LogBTF**, please directly contact the corresponding author [Prof. Zhi-Ping Liu](https://scholar.google.com/citations?user=zkBXb_kAAAAJ&hl=zh-CN&oi=ao) with the E-mail: zpliu@sdu.edu.cn
<!--END_SECTION:news-->


## Citation
Li, Lingyu, et al. "**LogBTF: Gene regulatory network inference using Boolean threshold networks model from single-cell gene expression data**." Submit to Briefings in Bioinformatics (https://academic.oup.com/bib). 


## Data
<!--START_SECTION:news-->
* **Supplementary Materials** file present the necessary **Additional files** contained in our work (**Table S1 and Table S2**).
* **Data** file give some necessary input/output files by the **R/Matlab/Python** codes. The subfile **DREAM3_RealData2** is the Matsumoto RNA-seq data, and the subfile **DREAM3_RealData16** is the Specific hHEP scRNA-seq data.
* Some of these input files only give the first few lines, limited by upload file size, but this does not affect the results of our work (**LogBTFs**).
* **Cytoscape** file, we give the inferred LMPP gene regulatory network. 
<!--END_SECTION:news-->


## R codes
The **serial number (1), (2), ..., (17)** represents the order in which the program runs in our work. All experiments are conducted on a workstation with two Xeon Gold 6226R CPUs and 256G of RAM

<!--START_SECTION:news-->
* (1) PermanceFunction.R  --  The function to evaluate inferred network performance.
* (2) mydataexamgather.R  --  The code for the Artificial data results. The data is generated by mydataexam.py code and saved in File Atta.
* (3) dataDREAMCLR.R  --  CLR method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (4) dataDREAMCLR_Yeast.R  --  CLR method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (5) dataDREAMGENIE3.R  --  GENIE3 method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (6) dataDREAMGENIE3_Yeast.R ---- GENIE3 method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (7) dataDREAMGlmRegAdj_Ecoli.R ---- LogBTFs method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (8) dataDREAMGlmRegAdj_Yeast.R -- LogBTFs method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (9) dataDREAMTARACNE.R – ARACNE method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (10) dataDREAMTARACNE_Yeast.R – ARACNE method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (11) dataDREAMTREGRESS.R -- TREGRESS method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (12) dataDREAMTREGRESS_Yeast.R -- TREGRESS method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (13)  dataDREAMGlmSINCERITIESadj_Ecoli.R  --  SINCERITIES method on DREAM3 dataset (Ecoli1 and Ecoli2).
* (14) dataDREAMGlmSINCERITIESadj_Yeast.R -- SINCERITIES method on DREAM3 dataset (Yeast1, Yeast2 and Yeast3).
* (15) Box_Group.R -- Visualize results of the DREAM3 dataset.
* (16) dataDREAMGlmRegAdjRealdata10fold.R, dataDREAMSINCERITIESRealdata(10fold).R, dataDREAMCLRRealdata.R, dataDREAMGENIE3Realdata.R, dataDREAMTARACNERealdata.R, dataDREAMTREGRESSRealdata.R  -- Six method on real dataset (Files Data2 and Data16).
* (17) LMPPGlmPenalty0803.R  -- LogBTFs method on LMPP dataset, compare with SINCERITIES.
<!--END_SECTION:news-->


## Python codes
<!--START_SECTION:news-->
* mydataexam.py – The generated data for the next time point of Artificial data.
<!--END_SECTION:news-->


## Matlab codes
<!--START_SECTION:news-->
* pseudo-data-order.m – Order the single-cells according to pseudo-time information.
<!--END_SECTION:news-->


## LogBTF (2022), Zhi-Ping Liu all rights reserved
This program package is supported by the copyright owners and coders "as is" and without warranty of any kind, express or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. In no event shall the copyright owner or contributor be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, without limitation, procurement of substitute goods or services; loss of use, data, or profits; or business interruption), regardless of the theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) for any use of the software, even if advised of the possibility of such damages.
