%计算带圆孔板在单向拉伸下的解析应力解
%这是无限板的解析解
%对于有限尺寸板体，在边界处会有一定误差
%在r=0处需要特别处理以避免除零错误
%应力值是在极坐标系下给出的
function stressAnalytical = computeAnalyticalStress(coords, radius, traction)
    x = coords(:, 1);
    y = coords(:, 2);
    r = sqrt(x.^2 + y.^2);% 计算径向距离
    theta = atan2(y, x);% 计算极角

    %计算径向应力
    sigma_rr = traction / 2 * (1 - radius^2 ./ r.^2) ...%均匀场项
             + traction / 2 * (1 - 4 * radius^2 ./ r.^2 + 3 * radius^4 ./ r.^4) .* cos(2 * theta);%扰动项
    %计算环向应力
    sigma_tt = traction / 2 * (1 + radius^2 ./ r.^2) ...%均匀场项
             - traction / 2 * (1 + 3 * radius^4 ./ r.^4) .* cos(2 * theta);%扰动项
    %计算剪切应力
    sigma_rt = -traction / 2 * (1 + 2 * radius^2 ./ r.^2 - 3 * radius^4 ./ r.^4) .* sin(2 * theta);
    
    %返回应力张量
    stressAnalytical = [sigma_rr, sigma_tt, sigma_rt];
end