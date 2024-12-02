%calculate Relative L2 and H1 Errors
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



%现在我们要循环不同的number of element
for case_idx = 1:num_cases                %这是循环的开头
    n_el = n_el_array(case_idx);          % 当前 number of elements
    pp = 2;                                % 多项式次数 (二次)
    n_en = pp + 1;                         % 每个element的本地节点数量 (3)
    n_np = n_el * pp + 1;                  % 总节点数
    h = 1.0 / n_el;                        % Mesh size
    h_values(case_idx) = h;                % 存储 mesh size
    
    % Setup the mesh
    x_coor = linspace(0,1,n_np);           % 节点坐标
    
    % 定义 the IEN matrix
    IEN = zeros(n_el, n_en);
    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end