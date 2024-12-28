%Aubin-Nische method 方法误差估计
%现在使用Aubin-Nische方法决定e在L2以及H1的收敛速度，该问题中m = 1,k = 2
clc;
clear;
close all;

%首先选择问题参数
k = 2%二次元素
m = 1
s_vals = 0:1%0<=s<=m
beta_vals = zeros(size(s_vals))%存储beta

%现在我们要假定一个参考解
%然后确定网格大小
%然后用一个for来迭代每个网格大小
%迭代的for里包括一个网格和元素矩阵来确定有限元解
%然后用解来确定L2和H1范数