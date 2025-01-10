% 定义边界条件
function [dispBC, forceBC] = defineBoundaryConditions(coords, traction, radius, sideLength)
    %coords: 节点坐标矩阵
    %traction: 外加荷载大小
    %radius: 半径参数
    %sideLength: 结构边长
    
    %设置容差值tol用于判断节点位置
    tol = 1e-6;
    %找出边界节点
    leftNodes = find(abs(coords(:, 1)) < tol);
    bottomNodes = find(abs(coords(:, 2)) < tol);
    rightNodes = find(abs(coords(:, 1) - sideLength) < tol);
    %错误检查
    if isempty(leftNodes) || isempty(bottomNodes) || isempty(rightNodes)
        error('Boundary nodes not found.');
    end
    %定义位移边界条件
    dispBC = [leftNodes, ones(length(leftNodes), 1), zeros(length(leftNodes), 1);
    bottomNodes, 2 * ones(length(bottomNodes), 1), zeros(length(bottomNodes), 1)];
    %定义力边界条件:
    forceBC = [rightNodes, ones(length(rightNodes), 1), traction * ones(length(rightNodes), 1)];
end
