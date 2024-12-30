function  [u, mesh_size] = fem_solver(x, y, elements, f)

    %现在解这个有限元问题
    num_nodes = length(x);               % Total number of nodes
    num_elements = size(elements, 1);    % Total number of elements
    K = sparse(num_nodes, num_nodes);    % Stiffness matrix
    F = zeros(num_nodes, 1);             % Load vecto

    %循环元素组装总体的K和F
    for e = 1:num_elements
        nodes = elements(e, : );          %当前元素的node index
        xe = x(nodes);                   %element nodes的x坐标
        ye = y(nodes);                   %element nodes的y坐标

       %计算每个元素的stiffness matrix 和 load vector
       %首先计算三角形区域
        A = abs(det([1, xe(1), ye(1); 1, xe(2), ye(2); 1, xe(3), ye(3)])) / 2;
    
        %基函数的梯度（对于线性三角形单元，是常数）
        b = [ye(2) - ye(3), ye(3) - ye(1), ye(1) - ye(2)] / (2 * A);
        c = [xe(3) - xe(2), xe(1) - xe(3), xe(2) - xe(1)] / (2 * A);
    
    %元素的stiffness matrix
    Ke = A * (b' * b + c' * c);

    %元素的load vector (使用中点积分法)
    x_mid = mean(x);
    y_mid = mean(y);
    Fe = A / 3 * f(x_mid, y_mid) * ones(3, 1);

    %随便跑了两个算例感觉不太对劲，这里似乎需要保证节点在边界内，如果超出就报错让我知道
    if any(nodes > num_nodes)
        error('node index exceeds the number of global nodes.')
    end

    %将每个元素累加到global stiffness matrix和load vector里
       for i = 1:3 % 三角形，三个nodes
              for j = 1:3
                  K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + Ke(i, j);
              end
               F(nodes(i)) = F(nodes(i)) + Fe(i); % 加上去
        end
    end
     
    %Dirchlet BC（边界上u = 0）
    boundary_nodes = find(x == 0 | x == max(x) | y == 0 | y == max(y));
    K(boundary_nodes, :) = 0;                     % 将边界节点对应的行置零
    K(boundary_nodes, boundary_nodes) = speye(length(boundary_nodes)); % 设置边界节点的对角线元素 to 1
    F(boundary_nodes) = 0;                        % 边界节点处的载荷值设为零

    %解Ku = F
    u = K \ F;

    num_elements = size(elements, 1);
    total_area = 0;

    for e = 1:num_elements
        nodes = elements(e, :);
        xe = x(nodes);
        ye = y(nodes);
        total_area = total_area + polyarea(xe, ye);
    end

    mesh_size = sqrt(total_area / num_elements); % 平均element size

end
