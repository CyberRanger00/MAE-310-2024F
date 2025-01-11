% Visualize应变
%目前只显示εxx分量
%如需显示其他应变分量需要修改代码
function visualizeStrain(coords, ien, strain)
    if size(strain, 2) > 1
        strain = strain(:, 1);
    end
    trisurf(ien, coords(:, 1), coords(:, 2), zeros(size(coords, 1), 1), strain, 'EdgeColor', 'none');
    colorbar;
    title('Strain Distribution (ε_{xx})');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal;% 保持坐标轴等比例
end