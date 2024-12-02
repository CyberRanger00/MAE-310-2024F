%前面还是一样的设置
clear all; clc; clf;%清理内存以及图片之类的

%首先定义exact solution 和其导数
u_exact = @(x) x.^5;     %定义exact solution，来源是对外源力积分两次
du_exact = @(x) 5*x.^4;  %exact derivative 求导

% 定义外源力和边界条件
f = @(x) -20*x.^3; % 给的source term
g = 1.0;           % u = g at x = 1 (Dirichlet)
h_bc = 0.0;        % -du/dx = h at x = 0 (Neumann)

% 定义要迭代的polynomial degrees
p_array = [2, 3]; % 2: Quadratic, 3: Cubic
num_p = length(p_array);

%为了迭代方便，定义一下number of element
n_el_array = [2, 4, 6, 8, 10, 12, 14, 16];
num_cases = length(n_el_array);

% 定义公差
tolerances = [1e-10, 1e-8, 1e-6];
num_tol = length(tolerances);

% 初始化cell arrays以储存误差，mesh size以及公差结果
L2_errors = cell(num_p,1);
H1_errors = cell(num_p,1);
h_values = cell(num_p,1);
GMRES_results = cell(num_p,1); % 存储公差信息

%这里我们先循环每个polynomial degree
for p_idx = 1:num_p
    pp = p_array(p_idx);            
    n_en = pp + 1;                     
    L2_errors{p_idx} = zeros(num_cases,1);
    H1_errors{p_idx} = zeros(num_cases,1);
    h_values{p_idx} = zeros(num_cases,1);
    GMRES_results{p_idx} = struct();   % 初始化一个struct来存储GMRE tolerance的结果
    %这里的设置大差不差
