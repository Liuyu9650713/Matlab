clc;clear;
%% 数据处理
[data,flag,active]=ID3_Handle();

%% 构造决策树
tree=ID3_Structure(data,flag,active);

%% 输出决策树
ID3_Printtree(tree);