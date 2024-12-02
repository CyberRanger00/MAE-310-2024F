%前面还是一样的设置
clear all; clc; clf;%清理内存以及图片之类的

%首先定义exact solution 和其导数
u_exact = @(x) x.^5;     %定义exact solution，来源是对外源力积分两次
du_exact = @(x) 5*x.^4;  %exact derivative 求导

% 定义外源力和边界条件，顺便定义polynomial degree
f = @(x) -20*x.^3; % 给的source term
g = 1.0;           % u = g at x = 1 (Dirichlet)
h_bc = 0.0;        % -du/dx = h at x = 0 (Neumann)
pp = 3;

% 定义要迭代的polynomial degrees
p_array = [2, 3]; % 2: Quadratic, 3: Cubic
num_p = length(p_array);

%为了迭代方便，定义一下number of element
n_el_array = [2, 4, 6, 8, 10, 12, 14, 16];
num_cases = length(n_el_array);

% 定义公差
tolerances = [1e-10, 1e-8, 1e-6];
num_tol = length(tolerances);

% 定义要进行实验的quadrature points的个数
quad_points_array = [1, 2, 3, 4, 5, 6];
num_quad = length(quad_points_array);

% 初始化cell arrays以储存误差，mesh size以及公差结果，这里不再需要cell，因为不用循环polynomial degrees
L2_errors = zeros(num_quad, num_cases);
H1_errors = zeros(num_quad, num_cases);
h_values = zeros(num_quad, num_cases);
GMRES_iterations = zeros(num_quad, num_cases); % 存储公差信息


% 初始化存储
GMRES_iterations = zeros(num_quad, num_cases);

% 为了画图方便，定义一下颜色和标志
colors = {'b','g','r','c','m','k'};    
markers = {'o','s','^','d','p','h'};   


%这里我们先循环每个polynomial degree，好像用不着,不过留着也不影响
%得删除，变成循环不同的quadrature point数
for q_idx = 1:num_quad
    n_quad = quad_points_array(q_idx); % 当前quadrature points数
    fprintf('\n--- Quadrature Points: %d ---\n', n_quad);
    

    % 这里我们循环每个number of elements，直接用前面的代码
    % 对每个number of element循环
for     case_idx = 1:num_cases
        n_el = n_el_array(case_idx);      % 当前elements数量
        n_np = n_el * pp +1;              % 总节点数
        h = 1.0 / n_el;                    % Mesh size
        h_values(q_idx, case_idx) = h;    % 存储mesh size
      
    % 设置网格
    x_coor = linspace(0,1,n_np);           % 节点坐标
    
    % 定义 the IEN matrix
    IEN = zeros(n_el, pp + 1);
    for ee = 1 : n_el
        for aa = 1 : (pp + 1)
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end

    % 设置 the ID array,这里不用改
    ID = 1 : n_np;
    ID(end) = 0; %  Dirichlet condition
    
    % 设置Quadrature(这里是循环的quadrature point数，所以不用定义点数)
   [xi, weight] = Gauss(n_quad, -1, 1);
    
    % 配置 the stiffness matrix 与 load vector
    n_eq = n_np - 1;                       % 方程个数（不含Dirichlet节点）
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq); % Sparse stiffness matrix
    F = zeros(n_eq, 1);                     % Load vector
    
    % 装配the stiffness matrix 与 load vector，这个整个应该没问题，高阶也不影响计算过程
    for ee = 1 : n_el
        k_ele = zeros(pp + 1, pp + 1); % Element stiffness matrix
        f_ele = zeros(pp + 1, 1);    % Element load vector
        
        x_ele = x_coor(IEN(ee,:)); % Coordinates of the current element's nodes

        % Quadrature loop
        for qua = 1 : n_quad
            % 初始化
            dx_dxi = 0.0;
            x_l = 0.0;
            
            % 在 quadrature point 计算 shape functions 及其导数 
            for aa = 1 : (pp + 1)
                N_a = PolyShape(pp, aa, xi(qua), 0);    % Shape function的值
                dN_a_dxi = PolyShape(pp, aa, xi(qua), 1);% Shape functionn 的导数
                x_l    = x_l    + x_ele(aa) * N_a;      % Physical coordinate
                dx_dxi = dx_dxi + x_ele(aa) * dN_a_dxi; % Jacobian derivative
            end
            dxi_dx = 1.0 / dx_dxi; % Inverse of Jacobian
            
            % 组装 element load vector 与 stiffness matrix
            for aa = 1 : (pp + 1)
                N_a = PolyShape(pp, aa, xi(qua), 0);       % Shape function 的值
                dN_a_dxi = PolyShape(pp, aa, xi(qua), 1);  % Shape function 的导数
                F_elem_increment = weight(qua) * N_a * f(x_l) * dx_dxi;
                f_ele(aa) = f_ele(aa) + F_elem_increment;
                
                for bb = 1 : (pp + 1)
                    dN_b_dxi = PolyShape(pp, bb, xi(qua), 1); % Shape function 的导数
                    K_ele_increment = weight(qua) * dN_a_dxi * dN_b_dxi * dxi_dx;
                    k_ele(aa, bb) = k_ele(aa, bb) + K_ele_increment;
                end
            end
        end
        
        % 装配成 global stiffness matrix 与 load vector
        for aa = 1 : (pp + 1)
            P = ID(IEN(ee,aa)); % Global DOF for node aa
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : (pp + 1)
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
    
    %在x=0处应用Neumann BC
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h_bc;

    %Direct Solver(应用LU 分解的解法)
    tic; % 开始计时
        d_direct = K \ F; % 直接解
        time_direct = toc; % 结束计时
        
        % 包含Dirichlet node的完整位移向量
        disp_direct = [d_direct; g];

     %Iterative Solver
     % 设置GMRES 参数
        restart = min(50, n_eq); % 确保重启次数不超过方程式数
        maxit = 1000;            % 最大迭代次数
        tol = 1e-10;             % 初始公差
        
        % 解
        [d_gmres, flag, relres, iter, resvec] = gmres(K, F, restart, tol, maxit);

        % 存
        if isscalar(iter)
            total_iter = iter;
        else
            total_iter = sum(iter);
        end
        GMRES_iterations(q_idx, case_idx) = total_iter;
        
        %比较解
        % 构造包含Dirichlet node的完整位移向量
        disp_gmres = [d_gmres; g];


        % 计算两种解法的差异
        diff = disp_gmres - disp_direct;
        error_diff = norm(diff, inf); % 最大差范数
        fprintf('  Elements: %d, GMRES Iterations: %d, Max Diff: %.2e\n', ...
                n_el, total_iter, error_diff);
        
        %计算误差，直接用前面的代码改
            % 这是一个后处理步骤：计算 L2 and H1 errors
    L2_error = 0.0;
    H1_error = 0.0;
    
    for ee = 1 : n_el
        x_ele = x_coor(IEN(ee,:));    % Element 的坐标
         u_ele = disp_direct(IEN(ee,:));   % Element 解的系数
        
        % 使用更高阶的quadrature进行误差积分
        [xi_err, weight_err] = Gauss(20, -1, 1); % 20-point quadrature
        
        for qua = 1 : length(weight_err)
            % 在 quadrature point计算 shape functions 及其导数 
            N = zeros(pp + 1,1);
            dN_dxi = zeros(pp + 1,1);
            for aa = 1 : pp + 1
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
            N = zeros(pp + 1,1);
            dN_dxi = zeros(pp + 1,1);
            for aa = 1 : pp + 1
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
    
    % 前面的部分都只是数值计算，完全不需要更改，只需要把数据类型变换一下以兼容两个polynomial degrees就好
    % 计算 relative errors，这里需要把数据类型变更一下，不然跑不动
        L2_errors(q_idx, case_idx) = sqrt(L2_error) / sqrt(L2_exact);
        H1_errors(q_idx, case_idx) = sqrt(H1_error) / sqrt(H1_exact);

        fprintf('    Relative Errors: L2_error = %.5e, H1_error = %.5e\n', ...
                L2_errors(q_idx, case_idx), H1_errors(q_idx, case_idx));
end
end

