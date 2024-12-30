function [x, y, elements] = generate_quad_mesh(L, H, Nx, Ny)
    % 生成quadrilateral形状的网格
    dx = L / Nx;
    dy = H / Ny;
    [x, y] = meshgrid(0:dx:L, 0:dy:H);
    x = x(:); y = y(:);
    elements = [];
    for i = 1:Ny
        for j = 1:Nx
            n1 = (i-1)*(Nx+1) + j;
            n2 = n1 + 1;
            n3 = n1 + Nx + 1;
            n4 = n3 + 1;
            elements = [elements; n1, n2, n4, n3];
        end
    end
end

function tri_elements = quad_to_tri(quad_elements)
    % 将quadrilateral网格变成两个三角形
    tri_elements = [];
    for i = 1:size(quad_elements, 1)
        quad = quad_elements(i, :);
        tri_elements = [tri_elements; quad(1), quad(2), quad(4)];
        tri_elements = [tri_elements; quad(1), quad(4), quad(3)];
    end
end