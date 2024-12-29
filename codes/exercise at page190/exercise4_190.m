%Aubin-Nische method 方法误差估计
%现在使用Aubin-Nische方法决定e在L2以及H1的收敛速度，该问题中m = 1,k = 2
clc;
clear;
close all;

%首先选择问题参数
k = 3 % 三次Hermite形函数 (k=3)
m = 2 % 二阶误差估计 (m=2)
s_vals = 0:1%0<=s<=m
beta_vals = zeros(size(s_vals))%存储beta

%现在我们要假定一个参考解
%然后确定网格大小
%然后用一个for来迭代每个网格大小
%迭代的for里包括一个网格和元素矩阵来确定有限元解
%然后用解来确定L2和H1范数

% 网格大小
h_vals = [1/2, 1/4, 1/8, 1/16]; % 不同的网格尺寸

% 假定一个参考解 u(x)
u_exact = @(x) sin(pi*x);  % 参考解为正弦函数

% 初始化误差数组
error_L2 = zeros(length(h_vals), 1);
error_H1 = zeros(length(h_vals), 1);
error_H2 = zeros(length(h_vals), 1);

% 计算不同网格大小下的误差
for h_idx = 1:length(h_vals)
    h = h_vals(h_idx);


    % 生成网格和元素矩阵
    x = 0:h:1; % 一维网格
    N = length(x) - 1; % 网格节点数量

    % 计算有限元解（这里以线性有限元为例）
    % 使用简单的线性插值作为有限元解近似
    u_approx = zeros(size(x)); % 近似解数组
    for i = 1:N
        u_approx(i) = u_exact(x(i)); % 这里直接取参考解的值，实际中为计算得到的近似解
    end

    % 计算误差 L2 范数
    error_L2(h_idx) = norm(u_approx - u_exact(x), 2) * h^2; % 近似解和参考解的L2误差
    
    % 计算误差 H1 范数
    grad_u_approx = gradient(u_approx, h); % 近似解的梯度
    grad_u_exact = gradient(u_exact(x), h); % 参考解的梯度
    error_H1(h_idx) = norm(grad_u_approx - grad_u_exact, 2) * h; % 近似解和参考解的H1误差
    
    % 计算误差 H2 范数
    grad2_u_approx = gradient(grad_u_approx, h); % 近似解的二阶导数
    grad2_u_exact = gradient(grad_u_exact, h); % 参考解的二阶导数
    error_H2(h_idx) = norm(grad2_u_approx - grad2_u_exact, 2) * h^0.5; % 近似解和参考解的H2误差
end


% 绘制结果
figure;
subplot(2,1,1);
loglog(h_vals, error_L2(:,1), 'o-', 'DisplayName', 'L2 Error');
hold on;
loglog(h_vals, error_H1(:,1), 'x-', 'DisplayName', 'H1 Error');
xlabel('h');
ylabel('Error');
title('Error Convergence in L2 and H1 norms');
legend;
grid on;

% 计算 beta 值
for s_idx = 1:length(s_vals)
    beta_vals(s_idx) = min(k+1-s_vals(s_idx), 2*(k+1-m)); % 计算beta值
end

subplot(2,1,2);
plot(s_vals, beta_vals, 's-', 'DisplayName', 'Beta Values');
xlabel('s');
ylabel('Beta');
title('Beta Values for Error Estimates');
grid on;
