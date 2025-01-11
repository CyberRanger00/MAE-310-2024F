function [grad] = computeAnalyticalGradient(coords)
    % 计算单元中心点
    xc = mean(coords(:,1));
    yc = mean(coords(:,2));
    
    % 转换为极坐标
    r = sqrt(xc^2 + yc^2);
    theta = atan2(yc, xc);
    
    % 问题参数
    R = 0.5;  % 孔的半径
    Tx = 10e3;  % 应用的拉力
    
    % 计算对r的偏导数
    dSdr = -Tx/2 * (2*R^2/r^3) - Tx/2 * (-8*R^2/r^3 + 12*R^4/r^5) * cos(2*theta);
    
    % 计算对θ的偏导数
    dSdt = Tx * (1 - 4*R^2/r^2 + 3*R^4/r^4) * (-sin(2*theta));
    
    % 转换到笛卡尔坐标
    grad = [dSdr*cos(theta) - dSdt*sin(theta)/r;
           dSdr*sin(theta) + dSdt*cos(theta)/r];
end