function [gradXX, gradYY, gradXY] = computeStressGradients(coords, stress, B)
    % 使用形函数导数计算应力梯度
    
    % 计算单元尺寸
    x21 = coords(2,1) - coords(1,1);
    x31 = coords(3,1) - coords(1,1);
    y21 = coords(2,2) - coords(1,2);
    y31 = coords(3,2) - coords(1,2);
    
    % 计算单元面积
    area = 0.5 * (x21 * y31 - x31 * y21);
    
    % 计算形函数导数
    dNdx = [y31, -y21, y21-y31] / (2 * area);
    dNdy = [-x31, x21, x31-x21] / (2 * area);
    
    % 计算各应力分量的梯度
    gradXX = [dNdx * stress(1); dNdy * stress(1)];  % σxx的梯度
    gradYY = [dNdx * stress(2); dNdy * stress(2)];  % σyy的梯度
    gradXY = [dNdx * stress(3); dNdy * stress(3)];  % τxy的梯度
end