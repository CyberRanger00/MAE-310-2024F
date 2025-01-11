%平面弹性问题分析
clc;clear;close all;

%参数设置
E = 1e9;   %杨氏模量
nu = 0.3;  %泊松比
planeStress = true;   %平面应力还是平面应变
radius = 0.5;  %洞的半径
sideLength = 4;  %平板的边长
traction = 10e3;    %load traction

    % 网格收敛性分析参数
    numRefinements = 4; % 网格加密次数
    L2Errors = zeros(numRefinements, 1);% 存储L2范数误差
    H1Errors = zeros(numRefinements, 1);% 存储H1范数误差
    elementSizes = zeros(numRefinements, 1);% 存储网格尺寸

    for i = 1:numRefinements
        %网格加密循环
        meshSize = sideLength / (2^(i + 1)); % 网格尺寸逐次减半
        %生成和读取网格
        generateGeoFile('geometry.geo', radius, sideLength, meshSize);
        system(sprintf('gmsh geometry.geo -2 -format msh2 -o geometry.msh -clmax %.2f', meshSize));
        [coords, ien] = readGmsh('geometry.msh');

        %求解有限元问题
        K = assembleStiffness(coords, ien, E, nu, planeStress);
        [dispBC, forceBC] = defineBoundaryConditions(coords, traction, radius, sideLength);
        [K_mod, F] = applyBoundaryConditions(K, dispBC, forceBC);
        U = K_mod \ F;
        [stressNumerical, strainNumerical] = postProcess(coords, ien, U, E, nu, planeStress);

        %计算误差
        stressAnalytical = computeAnalyticalStress(coords, radius, traction);
        [L2Errors(i), H1Errors(i)] = computeError(coords, ien, stressNumerical, stressAnalytical);

        %记录单元大小
        elementSizes(i) = meshSize;
    end

% Visualization
figure;
visualizeDisplacement(coords, ien, U);
saveas(gcf, 'DisplacementVisualization.png');
figure;
visualizeStress(coords, ien, stressNumerical);
saveas(gcf, 'StressVisualization.png');
figure;
visualizeStrain(coords, ien, strainNumerical);
saveas(gcf, 'StrainVisualization.png');

disp(L2Errors);
disp(H1Errors);

%绘制收敛性曲线
figure;
loglog(elementSizes, L2Errors, '-o', 'DisplayName', 'L2 Error');
hold on;
loglog(elementSizes, H1Errors, '-o', 'DisplayName', 'H1 Error');
xlabel('Element Size');
ylabel('Error');
title('Convergence Analysis');
legend;
grid on;
saveas(gcf, 'Convergence Analysis.png');