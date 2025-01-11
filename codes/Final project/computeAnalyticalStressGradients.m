function [gradXX, gradYY, gradXY] = computeAnalyticalStressGradients(coords, stress, nodes)
    % 计算解析应力梯度
    % 注意：这个函数需要根据具体问题的解析解来实现
    
    % 计算单元中心点
    xc = mean(coords(:,1));
    yc = mean(coords(:,2));
    
    % 转换为极坐标
    r = sqrt(xc^2 + yc^2);
    theta = atan2(yc, xc);
    
    % 问题参数
    R = 0.5;      % 孔的半径
    Tx = 10e3;    % 施加的拉力
    
    % 计算对r的偏导数
    dSrrdr = -Tx/2 * (2*R^2/r^3) - Tx/2 * (-8*R^2/r^3 + 12*R^4/r^5) * cos(2*theta);
    dSttdr = -Tx/2 * (R^2/r^3) - Tx/2 * (-12*R^4/r^5) * cos(2*theta);
    dSrtdr = -Tx/2 * (-4*R^2/r^3 + 12*R^4/r^5) * sin(2*theta);
    
    % 计算对θ的偏导数
    dSrrdt = Tx * (1 - 4*R^2/r^2 + 3*R^4/r^4) * (-sin(2*theta));
    dSttdt = -Tx * (1 + 3*R^4/r^4) * sin(2*theta);
    dSrtdt = -Tx * (1 + 2*R^2/r^2 - 3*R^4/r^4) * cos(2*theta);
    
    % 转换回笛卡尔坐标系
    gradXX = [dSrrdr*cos(theta) - dSrtdr*sin(theta);
              (dSrrdt - r*dSrtdr)*sin(theta)/r + (r*dSrrdr - dSrtdt)*cos(theta)/r];
              
    gradYY = [dSrrdr*sin(theta) + dSrtdr*cos(theta);
              -(dSrrdt - r*dSrtdr)*cos(theta)/r + (r*dSrrdr - dSrtdt)*sin(theta)/r];
              
    gradXY = [dSrtdr;
              dSrtdt/r];
end