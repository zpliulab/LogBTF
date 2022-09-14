
clear all
close all
addpath(genpath('./'))%Add the necessary path
fprintf('This is a demo run of GRISLI intented to reproduce the results of the article.\n')

%% Retrieve 取回 the SCODE datasets
filenumber=2;%Choose the murine data (2 for SCODE Data2) or the human data 
%(3 for SCODE Data3) or the pancreatic data (4 for scvelo dataset)
path_name=['SCODE-master/data' num2str(filenumber) '/']; 
data_matrix = dlmread([path_name 'datamatrix.txt'],'\t');
% data_matrix = dlmread([path_name 'data' num2str(filenumber) 'zinb.txt'],'\t');
% dlmread: 从由filename指定的具有分割标志的ASCII⽂件中读取数值数据
pseudotime_array = dlmread([path_name 'pseudotime.txt'],'\t');
if filenumber ~=4
	realtime_array = dlmread([path_name 'realtime.txt'],'\t');
    % no real time is available for scvelo data4
end
if filenumber ==4
	Vvelocity_matrix= dlmread([path_name 'velocity.txt'],' ');
    % scevlo velocity is only available for scvelo data4
end

A = dlmread([path_name 'A.txt'],'\t');%The real binary graph from the literature

t_array=pseudotime_array;%Choose the time label (real or pseudo) 
size(t_array);% 373
%(Advice: pseudo for SCODE Data2, real for SCODE Data3, pseudo for scvelo Data4)
X=[t_array,data_matrix'];

[~,I]=sort(X(:,1));
X=X(I,:);%X is the spacetime matrix of the data (of size Cx(1+G)), with the
%chosen (real or pseudo) time in the first column and the gene expression in the others 
