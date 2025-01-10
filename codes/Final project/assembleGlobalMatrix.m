%装配总体刚度矩阵
function K = assembleGlobalMatrix(K, Ke, nodes)
    %自由度计算
    dof = reshape([2 * nodes - 1; 2 * nodes], 1, []); %每个节点有两个自由度
    %确保自由度数量与单元刚度矩阵的大小匹配
    %防止组装过程中出现维度不匹配的错误
    if length(dof) ~= size(Ke, 1)
        error('DOF size mismatch with element stiffness matrix.');
    end
    %将单元刚度矩阵Ke的元素加到全局刚度矩阵K的对应位置
    K(dof, dof) = K(dof, dof) + Ke;
end
