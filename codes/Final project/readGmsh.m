%阅读gmesh生成的.msh文件
function [coords, ien] = readGmsh(filename)
    coords = [];
    ien = [];
    try
        fid = fopen(filename, 'r');
        isNodeSection = false;
        isElementSection = false;

        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, '$Nodes')
                isNodeSection = true;
                numNodes = str2double(fgetl(fid));
                coords = zeros(numNodes, 3);
                continue;
            end
            if contains(line, '$EndNodes')
                isNodeSection = false;
                continue;
            end
        end
    end
