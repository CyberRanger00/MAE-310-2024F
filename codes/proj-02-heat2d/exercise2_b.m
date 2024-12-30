%误差分析
clear;
clc;

%设置问题
L = 1;%定义域长度
H = 1;%定义域高度
Nx = 10;%x方向元素的数量
Ny = 10;%y方向元素的数量

u_exact = @(x, y) sin(pi*x).*sin(pi*y); % Exact solution
f = @(x, y) 2*pi^2 * sin(pi*x).*sin(pi*y); % RHS forcing function

%生成网格
[x, y, quad_elements] = generate_quad_mesh(L, H, Nx, Ny);

%将四边形变成三角形
tri_elements = quad_to_tri(quad_elements);

[u_quad, quad_mesh_size] = fem_solver(x, y, quad_elements, f);
[u_tri, tri_mesh_size] = fem_solver(x, y, tri_elements, f);

% 求导以计算Sobolev norms
grad_u_exact = @(x, y) [pi*cos(pi*x).*sin(pi*y), pi*sin(pi*x).*cos(pi*y)];


    display(u_quad);
    display(quad_mesh_size);
    display(u_tri);
    display(tri_mesh_size);

