% 单元的刚度矩阵
function Ke = elementStiffness(coords, E, nu, planeStress)
    %构建本构矩阵
    if planeStress%平面应力条件
        C = E / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    else          %平面应变条件
        C = E / ((1 + nu) * (1 - 2 * nu)) * [1 - nu, nu, 0; nu, 1 - nu, 0; 0, 0, (1 - 2 * nu) / 2];
    end
    %计算应力-应变矩阵B
    [B, detJ] = computeBMatrix(coords);
    Ke = B' * C * B * detJ;
end
