%二维有限元求解（quadrilateral elment,可能会被分成两个三角形来求解）
clear;
clc;

%设置问题
L = 1;%定义域长度
H = 1;%定义域高度
Nx = 10;%x方向元素的数量
Ny = 10;%y方向元素的数量

%现在需要按照示例做一个生成网格的函数对吧
%网断了，push不上去
%好了

%生成网格
[x, y, quad_elements] = generate_quad_mesh(L, H, Nx, Ny);

%将四边形变成三角形
tri_elements = quad_to_tri(quad_elements);

%现在解这个有限元问题





