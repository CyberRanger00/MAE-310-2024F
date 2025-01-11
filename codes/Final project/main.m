%平面弹性问题分析
clc;clear;close all;

%参数设置
E = 1e9;   %杨氏模量
nu = 0.3;  %泊松比
planeStress = true;   %平面应力还是平面应变
radius = 0.5;  %洞的半径
sideLength = 4;  %平板的边长
traction = 10e3;    %load traction

%主文件的生成以及阅读msh文件
meshSize = sideLength / (2^(i + 1)); % 减少mesh size
generateGeoFile('geometry.geo', radius, sideLength, meshSize);
system(sprintf('gmsh geometry.geo -2 -format msh2 -o geometry.msh -clmax %.2f', meshSize));
[coords, ien] = readGmsh('geometry.msh');
if isempty(coords) || isempty(ien)
    error('Mesh data could not be loaded. Check the Gmsh file and path.');
end

%现在就要解有限元问题了

%组装刚度矩阵
K = assembleStiffness(coords, ien, E, nu, planeStress);

%应用BC
[dispBC, forceBC] = defineBoundaryConditions(coords, traction, radius, sideLength);
[K_mod, F] = applyBoundaryConditions(K, dispBC, forceBC);

%解位移
U = K_mod \ F;

%后处理
[stressNumerical, strainNumerical] = postProcess(coords, ien, U, E, nu, planeStress);

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