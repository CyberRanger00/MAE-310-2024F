% Visualize Strain
function visualizeStrain(coords, ien, strain)
    % 首先检查strain的维度
    [~, numComponents] = size(strain);
    
    % 创建figure窗口
    figure('Position', [100, 100, 800, 600]);
    
    % 确定要显示的应变分量数量
    componentNames = {'ε_{xx}'};
    if numComponents > 1
        componentNames = {'ε_{xx}', 'ε_{yy}', 'γ_{xy}'};
    end
    
    % 根据应变分量的数量决定显示方式
    if numComponents == 1
        % 只有一个应变分量的情况
        trisurf(ien, coords(:,1), coords(:,2), zeros(size(coords,1),1), ...
               strain, 'EdgeColor', 'none');
        colorbar;
        title(['Strain Component: ' componentNames{1}]);
        xlabel('x');
        ylabel('y');
        view(2);
        axis equal;
    else
        % 多个应变分量的情况
        for i = 1:numComponents
            subplot(2, 2, i);
            trisurf(ien, coords(:,1), coords(:,2), zeros(size(coords,1),1), ...
                   strain(:,i), 'EdgeColor', 'none');
            colorbar;
            title(['Strain Component: ' componentNames{i}]);
            xlabel('x');
            ylabel('y');
            view(2);
            axis equal;
        end
        
        % 如果有足够的应变分量，添加等效应变
        if numComponents >= 3
            subplot(2, 2, 4);
            eqvStrain = sqrt(2/3 * ((strain(:,1) - strain(:,2)).^2 + ...
                         strain(:,2).^2 + strain(:,1).^2 + ...
                         3/2*strain(:,3).^2));
            trisurf(ien, coords(:,1), coords(:,2), zeros(size(coords,1),1), ...
                   eqvStrain, 'EdgeColor', 'none');
            colorbar;
            title('Equivalent Strain');
            xlabel('x');
            ylabel('y');
            view(2);
            axis equal;
        end
    end
    
    % 添加总标题
    sgtitle('Strain Distribution Analysis');
end