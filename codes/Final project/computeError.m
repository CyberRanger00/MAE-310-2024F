function [L2Error, H1Error] = computeError(coords, ien, stressNumerical, stressAnalytical)
    % 初始化误差
    L2Error = 0;
    H1Error = 0;
    totalArea = 0;
    
    % 对每个单元进行循环
    for e = 1:size(ien, 1)
        % 获得单元节点
        nodes = ien(e, :);
        elementCoords = coords(nodes, :);
        
        % 先检查应力向量的维度
        numStressComponents = size(stressNumerical, 2);  % 获取实际的应力分量数
        
        % 获取数值解和解析解
        sigmaNumerical = stressNumerical(e, :);     % 获取当前单元的数值应力
        sigmaAnalytical = mean(stressAnalytical(nodes, :), 1);  % 平均节点值
        
        % 计算B矩阵和雅可比行列式
        [B, detJ] = computeBMatrix(elementCoords);  % B矩阵大小为 3x6
        
        % 计算L2误差
        diff = sigmaNumerical - sigmaAnalytical;
        L2Error = L2Error + detJ * (diff * diff');
        
        % 提取B矩阵的导数行
        dNdx = B(1,:);  % x方向导数 1x6
        dNdy = B(2,:);  % y方向导数 1x6
        
        % 为每个应力分量计算梯度
        gradNumerical = [];  % 动态大小，根据应力分量数确定
        gradAnalytical = [];
        
        % 对每个应力分量分别计算x和y方向的梯度
        for i = 1:numStressComponents
            % 创建形函数向量
            N = [sigmaNumerical(i); sigmaNumerical(i); sigmaNumerical(i); 0; 0; 0];
            NA = [sigmaAnalytical(i); sigmaAnalytical(i); sigmaAnalytical(i); 0; 0; 0];
            
            % 计算数值解梯度
            gradNumerical = [gradNumerical; 
                           dNdx * N;    % x方向导数
                           dNdy * N];   % y方向导数
                           
            % 计算解析解梯度
            gradAnalytical = [gradAnalytical;
                            dNdx * NA;   % x方向导数
                            dNdy * NA];  % y方向导数
        end
        
        % 计算梯度差
        gradDiff = gradNumerical - gradAnalytical;
        
        % 计算H1误差
        H1Error = H1Error + detJ * (diff * diff' + gradDiff' * gradDiff);
        
        % 累加单元面积
        totalArea = totalArea + detJ;
    end
    
    % 最终归一化处理
    L2Error = sqrt(L2Error / totalArea);
    H1Error = sqrt(H1Error / totalArea);
end