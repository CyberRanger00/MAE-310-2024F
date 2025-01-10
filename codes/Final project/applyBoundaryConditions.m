%施加边界条件
function [K_mod, F] = applyBoundaryConditions(K, dispBC, forceBC)
    %创建一个全零的力向量 F,其大小与刚度矩阵的行数相同
    F = zeros(size(K, 1), 1);
    %应用力边界条件
    for i = 1:size(forceBC, 1)
        F(2 * forceBC(i, 1) - (2 - forceBC(i, 2))) = forceBC(i, 3);
    end
    %应用位移边界条件
    for i = 1:size(dispBC, 1)
        node = dispBC(i, 1);
        dof = 2 * node - (2 - dispBC(i, 2));
        K(:, dof) = 0; K(dof, :) = 0; K(dof, dof) = 1;
        F(dof) = dispBC(i, 3);
    end
    %返回修改后的刚度矩阵 K_mod 和力向量 F
    K_mod = K;
end
