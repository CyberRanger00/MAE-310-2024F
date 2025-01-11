% 定义边界条件
function [dispBC, forceBC] = defineBoundaryConditions(coords, traction,radius, sideLength)
    %coords: 节点坐标矩阵
    %traction: 外加荷载大小
    %radius: 半径参数
    %sideLength: 结构边长
    
    %设置容差值tol用于判断节点位置
    tol = 1e-6;
    %找出边界节点
    xSymNodes = find(abs(coords(:, 1)) < tol);  % x=0上的对称边界节点
    ySymNodes = find(abs(coords(:, 2)) < tol);  % y=0上的对称边界节点
    rightNodes = find(abs(coords(:, 1) - sideLength/2) < tol);  % 右边界节点
    topNodes = find(abs(coords(:, 2) - sideLength/2) < tol);    % 上边界节点
    %错误检查
    if isempty(xSymNodes) || isempty(ySymNodes) || isempty(rightNodes) || isempty(topNodes)
        error('Boundary nodes not found.');
    end
    %定义位移边界条件
    %x=0边界上的水平位移约束(u=0)
    %y=0边界上的竖直位移约束(v=0)
    dispBC = [xSymNodes, ones(length(xSymNodes), 1), zeros(length(xSymNodes), 1);
              ySymNodes, 2*ones(length(ySymNodes), 1), zeros(length(ySymNodes), 1)];
    
    %施加力边界条件
    %右边界上的水平拉力
    %边界上的竖直拉力
    forceBC = [rightNodes, ones(length(rightNodes), 1), traction*ones(length(rightNodes), 1);
               topNodes, 2*ones(length(topNodes), 1), traction*ones(length(topNodes), 1)];
end
