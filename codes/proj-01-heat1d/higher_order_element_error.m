clear all; clc; clf;%清理内存以及图片之类的

%首先定义exact solution 和其导数
u_exact = @(x) x.^5;     %定义exact solution，来源是对外源力积分两次
du_exact = @(x) 5*x.^4;  %exact derivative 求导

% 定义外源力和边界条件
f = @(x) -20*x.^3; % 给的source term
g = 1.0;           % u = g at x = 1 (Dirichlet)
h_bc = 0.0;        % -du/dx = h at x = 0 (Neumann)

%为了迭代方便，定义一下number of element
n_el_array = [2, 4, 6, 8, 10, 12, 14, 16];
num_cases = length(n_el_array);

%初始化数组以存储 errors 和 mesh sizes
L2_errors = zeros(num_cases,1);
H1_errors = zeros(num_cases,1);
h_values = zeros(num_cases,1);

% 定义要迭代的polynomial degrees 
p_array = [2, 3]; % 2: Quadratic, 3: Cubic
num_p = length(p_array);

