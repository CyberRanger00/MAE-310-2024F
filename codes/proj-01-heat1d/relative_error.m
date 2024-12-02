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

    % 设置 the ID array 
    ID = 1 : n_np;
    ID(end) = 0; %  Dirichlet condition
    
    % 设置Quadrature（使用10点Gauss quadrature）
    n_int = 10;
    [xi, weight] = Gauss(n_int, -1, 1);
    
    % 配置 the stiffness matrix 与 load vector
    n_eq = n_np - 1;                       % 方程个数（不含Dirichlet节点）
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq); % Sparse stiffness matrix
    F = zeros(n_eq, 1);                     % Load vector
    
    % 装配the stiffness matrix 与 load vector
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % Element stiffness matrix
        f_ele = zeros(n_en, 1);    % Element load vector
        
        x_ele = x_coor(IEN(ee,:)); % Coordinates of the current element's nodes

        % Quadrature loop
        for qua = 1 : n_int
            % 初始化
            dx_dxi = 0.0;
            x_l = 0.0;
            
            % 在 quadrature point 计算 shape functions 及其导数 
            for aa = 1 : n_en
                N_a = PolyShape(pp, aa, xi(qua), 0);    % Shape function的值
                dN_a_dxi = PolyShape(pp, aa, xi(qua), 1);% Shape functionn 的导数
                x_l    = x_l    + x_ele(aa) * N_a;      % Physical coordinate
                dx_dxi = dx_dxi + x_ele(aa) * dN_a_dxi; % Jacobian derivative
            end
            dxi_dx = 1.0 / dx_dxi; % Inverse of Jacobian
            
            % 组装 element load vector 与 stiffness matrix
            for aa = 1 : n_en
                N_a = PolyShape(pp, aa, xi(qua), 0);       % Shape function 的值
                dN_a_dxi = PolyShape(pp, aa, xi(qua), 1);  % Shape function 的导数
                F_elem_increment = weight(qua) * N_a * f(x_l) * dx_dxi;
                f_ele(aa) = f_ele(aa) + F_elem_increment;
                
                for bb = 1 : n_en
                    dN_b_dxi = PolyShape(pp, bb, xi(qua), 1); % Shape function 的导数
                    K_ele_increment = weight(qua) * dN_a_dxi * dN_b_dxi * dxi_dx;
                    k_ele(aa, bb) = k_ele(aa, bb) + K_ele_increment;
                end
            end
        end
        