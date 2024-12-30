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
