% 用gmsh生成.geo文件
function generateGeoFile(filename, radius, sideLength, meshSize)
    fid = fopen(filename, 'w');
    fprintf(fid, 'SetFactory("OpenCASCADE");\n');%设置Gmsh
    %创建基本几何形状
    fprintf(fid, 'Rectangle(1) = {0, 0, 0, %.2f, %.2f};\n', sideLength / 2, sideLength / 2);
    fprintf(fid, 'Disk(2) = {0, 0, 0, %.2f};\n', radius);
    %从矩形中减去圆形,创建带孔板
    fprintf(fid, 'BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }\n');
    %设置网格参数
    fprintf(fid, 'Mesh.CharacteristicLengthMax = %.2f;\n', meshSize);
    fprintf(fid, 'Physical Surface("Plate") = {1};\n'); % 定义板的表面
    fprintf(fid, 'Physical Curve("Boundary") = {3, 2, 1};\n');% 定义边界
    fclose(fid);
end