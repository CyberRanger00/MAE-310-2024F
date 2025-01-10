% 装配刚度矩阵
function K = assembleStiffness(coords, ien, E, nu, planeStress)
    numNodes = size(coords, 1);% 获取节点总数
    K = zeros(2 * numNodes);% 创建2n×2n的零矩阵（n为节点数）
                            % 每个节点有x和y两个自由度
    for e = 1:size(ien, 1)  % 遍历每个单元
        nodes = ien(e, :);  % 获取当前单元的节点编号
        if any(nodes <= 0) || any(nodes > numNodes)
            warning('Invalid element connectivity for element %d. Skipping.', e);
            continue;
        end
        %这里先留两个函数
        Ke = elementStiffness(coords(nodes, :), E, nu, planeStress); % 计算单元刚度矩阵
    end
end
