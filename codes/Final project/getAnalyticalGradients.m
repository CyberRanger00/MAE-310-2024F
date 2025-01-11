function [gradXX, gradYY, gradXY] = getAnalyticalGradients(coords)
  % 计算中心点坐标
    xc = mean(coords(:,1));
    yc = mean(coords(:,2));
    
    % 转换为极坐标
    r = sqrt(xc^2 + yc^2);
    theta = atan2(yc, xc);
    
    % 应力梯度（这里使用简化的梯度计算）
    if r > radius
        Tx = 10e3;  % 施加的拉力
        
        % 对r的偏导数
        dSdr = -Tx/2 * (2*radius^2/r^3) - ...
               Tx/2 * (-8*radius^2/r^3 + 12*radius^4/r^5) * cos(2*theta);
        
        % 对θ的偏导数
        dSdt = Tx * (1 - 4*radius^2/r^2 + 3*radius^4/r^4) * (-sin(2*theta));
        
        % 转换回笛卡尔坐标
        grad = [dSdr*cos(theta) - dSdt*sin(theta)/r;
               dSdr*sin(theta) + dSdt*cos(theta)/r];
    else
        grad = [0; 0];
    end
end