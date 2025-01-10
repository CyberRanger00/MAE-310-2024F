% Visualize应力
function visualizeStress(coords, ien, stress)
    if size(coords, 2) < 3
        coords(:, 3) = 0; % 如果没有z坐标,添加0值
    end
    %检查应力向量的长度是否与单元数量匹配
    if size(stress, 1) ~= size(ien, 1)
        error('Size of stress (%d) does not match number of elements (%d).', size(stress, 1), size(ien, 1));
    end

    % 绘制应力分布
    trisurf(ien, coords(:, 1), coords(:, 2), coords(:, 3), stress, 'EdgeColor', 'none');
    colorbar;% 添加颜色条
    title('Stress Visualization (\sigma_{xx})');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(2); % 设置为俯视图
    axis equal;% 保持坐标轴等比例
end