%% 2022.11.28 based on the following knowledge
% Demo for GRISLI
% Pierre-Cyril Aubin-Frankowski, 2018
% Conduct on Windows Matlab R2017b


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For simulated single-cell gene expression data

%% initical
clc; clear; close all; format long;

%% Please make sure to be located in GRISLI/.
addpath(genpath('D:\E\博士\Matlab_code\Boolean\GRISLI\DREAMsc'))
fprintf('This is a demo run of GRISLI intended to reproduce the results of the article.\n')

%% parameter
fileNUM = [1, 2, 3];
geneNUM = [10, 50, 100];

%% Ecoli
dataname = 'Ecoli';

for i = 1:2
    for j = 1:3
        [AUROC_GRISLI, AUPR_GRISLI, elapsedTime] = GRISLIfuction(fileNUM(i),geneNUM(j),dataname);
        ROCPRTEcoli((i-1)*3+j,1) = AUROC_GRISLI;
        ROCPRTEcoli((i-1)*3+j,2) = AUPR_GRISLI;
        ROCPRTEcoli((i-1)*3+j,8) = elapsedTime;
    end
        
end

%% save  row:1-10, 1-50, 1-100, 2-10, 2-50, 2-100
csvwrite(['D:\E\博士\Matlab_code\Boolean\GRISLI\DREAMscResult\ROCPRTEcoli.csv'],ROCPRTEcoli);

%% Yeast
dataname = 'Yeast';
for i = 1:3
    for j = 1:3
        
        [AUROC_GRISLI, AUPR_GRISLI, elapsedTime] = GRISLIfuction(fileNUM(i),geneNUM(j),dataname);
        ROCPRTYeast((i-1)*3+j,1) = AUROC_GRISLI;
        ROCPRTYeast((i-1)*3+j,2) = AUPR_GRISLI;
        ROCPRTYeast((i-1)*3+j,8) = elapsedTime;
    end
        
end
%% save
csvwrite(['D:\E\博士\Matlab_code\Boolean\GRISLI\DREAMscResult\ROCPRTYeast.csv'],ROCPRTYeast);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For real scRNA-seq data 


clear all; close all
addpath(genpath('D:\E\博士\Matlab_code\Boolean\GRISLI\'))
fprintf('This is a demo run of GRISLI intented to reproduce the results of the article.\n')

%% (Matsumoto (Treutlein et al., 2016) -- Retrieve the SCODE datasets
% filenumber=2;  %Choose the murine data (2 for SCODE Data2) or the human data 
% path_name=['SCODE-master/data' num2str(filenumber) '/']; 
% data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
% pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
% 
% A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature
% t_array=pseudotime_array;%Choose the time label (real or pseudo) 
% %(Advice: pseudo for SCODE Data2, real for SCODE Data3, pseudo for scvelo Data4)
% X=[t_array,data_matrix'];
% 
% [~,I]=sort(X(:,1));
% X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
% %chosen (real or pseudo) time in the first column and the gene expression in the others 

%% hHEP (Campet al., 2017) -- data 16  
tic

filenumber=16;
path_name=['SCODE-master/data' num2str(filenumber) 'old/']; 
data_matrix = dlmread([path_name 'data.txt'],'\t');
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature
t_array=pseudotime_array;%Choose the time label (real or pseudo) 
X=[t_array,data_matrix'];
[~,I]=sort(X(:,1));
X=X(I,:);

%% TESTING GRISLI  --  use R=L=10 for saving runtimes
Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The kernel we use
R=10;   
L_array=10;  
alpha_min=.4;
saveResultsBool=false;
nbtry=1;

%% Win 10 OK ， Linux nowork
saveFileName='AUROC_files/GRISLI_Data2_Gvelo.txt';
[A_app_array_ind,AUROC_GRISLI, AUPR_GRISLI, elapsedTime, TPR_array_area, FPR_array_area,...
    PPV_array_area,A_app_array_Rnk0]=Test_GRISLI_realdataLLY(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,false);%,V0,nbtry
ROCPlot(FPR_array_area,TPR_array_area,"GRISLI");

%% save
ROCPRT(1,1) = AUROC_GRISLI;
ROCPRT(1,2) = AUPR_GRISLI;
ROCPRT(1,8) = elapsedTime;

toc
disp(['运行时间: ',num2str(toc)])
%% ROCPRTdata5 and ROCPRTdata2 are both ROCPRTdata2
% csvwrite(['D:\E\博士\Matlab_code\Boolean\GRISLI\DREAMscResult\ROCPRTdata5.csv'],ROCPRT);

%% data 16 RL=10
csvwrite(['D:\E\博士\Matlab_code\Boolean\GRISLI\DREAMscResult\ROCPRTdata16RL10.csv'],ROCPRT);
