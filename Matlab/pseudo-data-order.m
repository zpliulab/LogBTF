
clear all
close all
addpath(genpath('./'))%Add the necessary path
fprintf('This is a demo run of GRISLI intented to reproduce the results of the article.\n')

%% Retrieve the SCODE datasets
filenumber=2;%Choose the murine data (2 for SCODE Data2) or the human data 
path_name=['SCODE-master/data' num2str(filenumber) '/']; 
data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature
t_array=pseudotime_array;%Choose the time label (real or pseudo) 
X=[t_array,data_matrix'];
[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 
