function [AUROC_GRISLI, AUPR_GRISLI, elapsedTime] = GRISLIfuction(filenumber,genenumber,dataname)

% Demo for GRISLI
% Pierre-Cyril Aubin-Frankowski, 2018
%% Retrieve È¡»Ø the SCODE datasets
path_name=[dataname num2str(filenumber) '-', num2str(genenumber) '/']; 
X = dlmread([path_name 'data.txt'],'\t'); %X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 
A = dlmread([path_name 'A.txt'],'\t'); %The real binary graph from the literature

%% TESTING GRISLI  
Alpha=@(Kx,Dt,sigx,sigt)exp(-(Kx.^2)/(2*sigx^2)).*exp(-(Dt.^2)/(2*sigt^2)).*(Dt.^2);%The kernel we use
R=100;
L_array=100;
alpha_min=.4;
saveResultsBool=false;
nbtry=1;

%% random seed
s = rng;
r = rand(1,1);
rng(s)

%% Win 10 good £¬ Linux error  perform GRISLI
saveFileName='AUROC_files/SimulatedData.txt';
[AUROC_GRISLI, AUPR_GRISLI, elapsedTime, TPR_array_area, FPR_array_area,...
    PPV_array_area,A_app_array_Rnk0]=Test_GRISLI_realdata(A,X,L_array,Alpha,R,alpha_min,saveResultsBool, saveFileName,false);%,V0,nbtry
% ROCPlot(FPR_array_area,TPR_array_area,"GRISLI");

end

