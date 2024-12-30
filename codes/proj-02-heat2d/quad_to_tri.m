function tri_elements = quad_to_tri(quad_elements)
    % 将quadrilateral网格变成两个三角形
    tri_elements = [];
    for i = 1:size(quad_elements, 1)
        quad = quad_elements(i, :);
        tri_elements = [tri_elements; quad(1), quad(2), quad(4)];
        tri_elements = [tri_elements; quad(1), quad(4), quad(3)];
    end
end