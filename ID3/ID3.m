clc;clear;
%% ���ݴ���
[data,flag,active]=ID3_Handle();

%% ���������
tree=ID3_Structure(data,flag,active);

%% ���������
ID3_Printtree(tree);