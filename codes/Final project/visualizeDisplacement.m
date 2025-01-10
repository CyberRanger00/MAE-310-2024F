% Visualize 位移
function visualizeDisplacement(coords, ien, U)
%coords: 节点坐标矩阵
%ien: 单元节点连接关系矩阵
%U: 节点位移
    %检查位移向量长度是否为偶数
    if mod(length(U), 2) ~= 0
        error('Length of U must be even.');
    end
    numNodes = size(coords, 1);
    %检查位移向量大小是否与节点数匹配(每个节点有x,y两个自由度)
    if length(U) ~= 2 * numNodes
        error('Size of U (%d) does not match 2 * numNodes (%d).', length(U), 2 * numNodes);
    end
    scale = 10;% 放大变形使其更容易观察
    displacements = reshape(U, 2, [])'; % 将位移向量重组为nx2矩阵
    displacedCoords = coords(:, 1:2) + scale * displacements; % 计算变形后的坐标
    triplot(ien, coords(:, 1), coords(:, 2), 'b'); % 绘制原始形状
    hold on;
    triplot(ien, displacedCoords(:, 1), displacedCoords(:, 2), 'r');% 绘制变形后形状

    %添加图形标注
    title('Displacement Visualization');
    legend('Original', 'Deformed');
    hold off;
end
