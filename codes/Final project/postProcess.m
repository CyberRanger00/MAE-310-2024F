function [stress, strain] = postProcess(coords, ien, U, E, nu, planeStress)
%coords: 节点坐标
%ien: 单元节点连接关系矩阵
%U: 整体位移向量
%E: 弹性模量
%nu: 泊松比
%planeStress: 平面应力状态标志(true/false)
    %初始化结果数组:
    stress = zeros(size(ien, 1), 1); % 存储每个单元的应力
    strain = zeros(size(ien, 1), 3);% 存储每个单元的应变(εx, εy, γxy)
    %对每个单元进行循环计算
    for e = 1:size(ien, 1)
        nodes = ien(e, :); % 获取当前单元的节点编号
        Ue = U([2 * nodes - 1, 2 * nodes]);% 提取单元节点位移
        %计算应变-位移矩阵B:
        [B, ~] = computeBMatrix(coords(nodes, :));
        %根据问题类型构建弹性矩阵C
        if planeStress
            %平面应力状态
            C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
        else
            %平面应变状态
            C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
        end
        strain(e, :) = B * Ue;% 计算应变
        stress(e) = C(1, :) * strain(e, :)';% 计算应力
    end
end
