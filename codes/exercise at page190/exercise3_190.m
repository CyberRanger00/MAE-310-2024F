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
