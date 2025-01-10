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
            %阅读node
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
      
        %end
        %阅读element
            if contains(line, '$Elements')
                isElementSection = true;
                numElements = str2double(fgetl(fid));
                ien = zeros(numElements, 3);
                continue;
            end
            if contains(line, '$EndElements')
                isElementSection = false;
                continue;
            end
         %end
         %阅读
            if isNodeSection
                data = sscanf(line, '%d %f %f %f');% 解析节点数据行
                if numel(data) == 4
                    coords(data(1), :) = data(2:4);% 存储节点坐标 (x,y,z)
                end
            end
            if isElementSection
                data = sscanf(line, '%d %d %d %d %d %d %d %d'); % 解析单元数据行
                if numel(data) >= 5 && data(2) == 2 % 检查是否为二维单元(type=2)
                    ien(end + 1, :) = data(6:8); % 存储单元的三个节点编号
                end
            end

        ien = ien(any(ien, 2), :); % 移除空行
        fclose(fid); % 关闭文件
    catch
        fclose('all'); % 发生错误时关闭所有打开的文件
        warning('Failed to read mesh data.'); % 发出警告
    end
