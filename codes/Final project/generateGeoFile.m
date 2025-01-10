% 用gmsh生成.geo文件
function generateGeoFile(filename, radius, sideLength)
    fid = fopen(filename, 'w');
    fprintf(fid, 'SetFactory("OpenCASCADE");\n');
    fprintf(fid, 'Rectangle(1) = {0, 0, 0, %.2f, %.2f};\n', sideLength, sideLength);
    fprintf(fid, 'Disk(2) = {0, 0, 0, %.2f};\n', radius);
    fprintf(fid, 'BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }\n');
    fprintf(fid, 'Physical Surface("Plate") = {1};\n');
    fprintf(fid, 'Physical Curve("Boundary") = {4, 3, 2, 1};\n');
    fclose(fid);
end