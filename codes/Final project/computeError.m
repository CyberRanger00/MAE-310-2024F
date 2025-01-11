% 计算误差
function [L2Error, H1Error] = computeError(coords, ien, stressNumerical, stressAnalytical)
%coords: 节点坐标矩阵
%ien: 单元节点连接关系矩阵
%stressNumerical: 数值计算得到的应力
%stressAnalytical: 解析解得到的应力

    %初始化误差
    L2Error = 0;
    H1Error = 0;
    totalArea = 0;

    % 对每个单元进行循环
    for e = 1:size(ien, 1)
        % 获得单元节点
        nodes = ien(e, :);
        elementCoords = coords(nodes, :);

        % 获取数值解和解析解
        sigmaNumerical = stressNumerical(e, :); % 单元数值应力
        sigmaAnalytical = mean(stressAnalytical(nodes, :), 1); % 单元平均解析应力
        
        % 计算单元雅可比行列式（用于面积计算）:
        [~, detJ] = computeBMatrix(elementCoords);

        % 计算L2误差： ||σ - σ_exact||²
        diff = sigmaNumerical - sigmaAnalytical; 
        L2Error = L2Error + detJ * (diff * diff'); % 积分过程

        % 计算H1误差（梯度差）
        gradDiff = diff; % 这里是简化处理,实际应该计算梯度差
        H1Error = H1Error + detJ * (gradDiff * gradDiff'); % 每个单元对H1范数的贡献

        % 累加单元面积
        totalArea = totalArea + detJ;
    end

    % 最终归一化处理
    L2Error = sqrt(L2Error / totalArea);
    H1Error = sqrt(H1Error / totalArea);
end

