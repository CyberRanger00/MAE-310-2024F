% 计算误差
function [L2Error, H1Error] = computeError(coords, ien, stressNumerical, stressAnalytical)
%coords: 节点坐标矩阵
%ien: 单元节点连接关系矩阵
%stressNumerical: 数值计算得到的应力
%stressAnalytical: 解析解得到的应力
    % 初始化误差
    L2Error = 0;
    H1Error = 0;
    gradError = 0;    % 梯度误差项
    totalArea = 0;    % 总面积
    
    % 遍历所有单元
    for e = 1:size(ien, 1)
        % 获取单元节点
        nodes = ien(e, :);
        elementCoords = coords(nodes, :);
        
        % 获取该单元的应力值
        sigmaNum = stressNumerical(e, :);         % [σxx, σyy, τxy]
        sigmaAnalytical = mean(stressAnalytical(nodes, :), 1);  % 取节点平均值作为解析值
        
        % 计算单元的B矩阵和雅可比行列式
        [B, detJ] = computeBMatrix(elementCoords);
        
        % 计算应力差值用于L2误差
        stressDiff = sigmaNum - sigmaAnalytical;
        
        % 计算该单元对L2误差的贡献
        L2Error = L2Error + detJ * (stressDiff * stressDiff');
        
        % 计算应力梯度用于H1误差
        % 计算数值解的梯度
        [gradNumXX, gradNumYY, gradNumXY] = computeStressGradients(elementCoords, sigmaNum, B);
        
        % 计算解析解的梯度
        [gradAnalXX, gradAnalYY, gradAnalXY] = computeAnalyticalStressGradients(...
            elementCoords, sigmaAnalytical, nodes);
        
        % 计算梯度差值
        gradDiffXX = gradNumXX - gradAnalXX;
        gradDiffYY = gradNumYY - gradAnalYY;
        gradDiffXY = gradNumXY - gradAnalXY;
        
        % 将梯度误差加入H1误差
        gradError = gradError + detJ * (gradDiffXX.^2 + gradDiffYY.^2 + gradDiffXY.^2);
        
        % 累加面积
        totalArea = totalArea + detJ;
    end
    
    % 归一化误差
    L2Error = sqrt(L2Error / totalArea);
    H1Error = sqrt((L2Error^2 + gradError) / totalArea);
end