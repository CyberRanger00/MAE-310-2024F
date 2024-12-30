function [e0, e1] = compute_error(x, y, elements, uh, u_exact, grad_u_exact)
%计算L2与H1 error
e0 = 0; e1 = 0;
for e = 1:size(elements, 1)
        nodes = elements(e, :);
        xe = x(nodes); ye = y(nodes);
        uh_e = uh(nodes);

        % Quadrature points
        xq = mean(xe);
        yq = mean(ye);

        % quadrature points处的精确解和梯度 
        u_exact_q = u_exact(xq, yq);
        grad_u_exact_q = grad_u_exact(xq, yq);

         % 近似解和梯度
        uh_q = mean(uh_e);
        grad_uh_q = [0, 0]; %梯度的近似

        % L2 error
        e0 = e0 + (uh_q - u_exact_q)^2;

        % H1 error
        e1 = e1 + sum((grad_uh_q - grad_u_exact_q).^2);
        
        e0 = sqrt(e0);
        e1 = sqrt(e1);
end
