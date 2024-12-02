%calculate Relative L2 and H1 Error
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
        
        % 装配成 global stiffness matrix 与 load vector
        for aa = 1 : n_en
            P = ID(IEN(ee,aa)); % Global DOF for node aa
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb)); % Global DOF for node bb
                    if(Q > 0)
                        % 确保P与Q在边界内
                        if P > n_eq || Q > n_eq
                            error(['Index exceeds matrix dimensions: P = ', num2str(P), ', Q = ', num2str(Q), ', n_eq = ', num2str(n_eq)]);
                        end
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        % Dirichlet boundary condition (u = g at Q=0)
                        F(P) = F(P) - k_ele(aa, bb) * g;
                    end
                end
            end
        end
    end
    
    
    % Apply Neumann boundary condition at x=0
    % Neumann BC： -du/dx = h_bc; since h_bc =0,这一步好像没有影响
    % 若 h_bc !=0, 就改一下
    % 本题中, h_bc=0, 因此无变化
    % 为了规整还是定义一下
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h_bc;
    
    % 解出 K * d = F
    d_temp = K \ F;
    
    % 构造完整的 displacement vector，包括 Dirichlet node
    disp_num = [d_temp; g];
    
    % 这是一个后处理步骤：计算 L2 and H1 errors
    L2_error = 0.0;
    H1_error = 0.0;
    
    for ee = 1 : n_el
        x_ele = x_coor(IEN(ee,:));    % Element 的坐标
        u_ele = disp_num(IEN(ee,:));   % Element 解的系数
        
        % 使用更高阶的quadrature进行误差积分
        [xi_err, weight_err] = Gauss(20, -1, 1); % 20-point quadrature
        
        for qua = 1 : length(weight_err)
            % 在 quadrature point计算 shape functions 及其导数 
            N = zeros(n_en,1);
            dN_dxi = zeros(n_en,1);
            for aa = 1 : n_en
                N(aa) = PolyShape(pp, aa, xi_err(qua), 0);
                dN_dxi(aa) = PolyShape(pp, aa, xi_err(qua), 1);
            end
            
            % 在 quadrature point 计算其物理坐标和解
            x_q = x_ele * N;               % (1 x 3) * (3 x 1) = scalar
            u_h = u_ele' * N;              % (1 x 3) * (3 x 1) = scalar
            du_h_dxi = u_ele' * dN_dxi;    % (1 x 3) * (3 x 1) = scalar
            
            % 计算x的导数
            dx_dxi = x_ele * dN_dxi;       % (1 x 3) * (3 x 1) = scalar
            du_h_dx = du_h_dxi / dx_dxi;
            
            % Exact solution , derivative
            u_ex = u_exact(x_q);
            du_ex = du_exact(x_q);
            
            % 计算误差
            L2_error = L2_error + weight_err(qua) * (u_h - u_ex)^2 * dx_dxi;
            H1_error = H1_error + weight_err(qua) * (du_h_dx - du_ex)^2 * dx_dxi;
        end
    end
    
    % 误差归一
    % 计算 exact solution 的norms
    L2_exact = 0.0;
    H1_exact = 0.0;
    
    for ee = 1 : n_el
        x_ele = x_coor(IEN(ee,:));    % Element coordinate
        
        % 使用更高阶的quadrature
        [xi_norm, weight_norm] = Gauss(20, -1, 1); % 20-point quadrature
        
        for qua = 1 : length(weight_norm)
            % 在 quadrature point计算 shape functions 及其 derivatives 
            N = zeros(n_en,1);
            dN_dxi = zeros(n_en,1);
            for aa = 1 : n_en
                N(aa) = PolyShape(pp, aa, xi_norm(qua), 0);
                dN_dxi(aa) = PolyShape(pp, aa, xi_norm(qua), 1);
            end
            
            % 计算坐标
            x_q = x_ele * N;               % (1 x 3) * (3 x 1) = scalar
            
            % Exact solution and derivative
            u_ex = u_exact(x_q);
            du_ex = du_exact(x_q);
            
            % Compute Jacobian
            dx_dxi = x_ele * dN_dxi;       % (1 x 3) * (3 x 1) = scalar
            
            % 累积exact norms
            L2_exact = L2_exact + weight_norm(qua) * u_ex^2 * dx_dxi;
            H1_exact = H1_exact + weight_norm(qua) * du_ex^2 * dx_dxi;
        end
    end
    
    % 计算 relative errors
    L2_errors(case_idx) = sqrt(L2_error) / sqrt(L2_exact);
    H1_errors(case_idx) = sqrt(H1_error) / sqrt(H1_exact);
    
    fprintf('n_el = %d, h = %.5f, L2_error = %.5e, H1_error = %.5e\n', ...
            n_el, h, L2_errors(case_idx), H1_errors(case_idx));
end


% 画图
figure;
loglog(h_values, L2_errors, '-o','LineWidth',2,'MarkerSize',8);
hold on;
loglog(h_values, H1_errors, '-s','LineWidth',2,'MarkerSize',8);
xlabel('Mesh size h');
ylabel('Relative Error');
title('Relative L_2 and H_1 Errors vs Mesh Size');
legend('L_2 Error','H_1 Error','Location','Best');
grid on;

% 得到斜率
% 拟合函数
p_L2 = polyfit(log(h_values), log(L2_errors), 1);
p_H1 = polyfit(log(h_values), log(H1_errors), 1);

% 展示斜率
fprintf('Slope of L2 error: %.2f\n', p_L2(1));
fprintf('Slope of H1 error: %.2f\n', p_H1(1));
