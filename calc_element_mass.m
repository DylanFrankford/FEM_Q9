function M_el = calc_element_mass(el_nodes, rho, thickness, gauss_points, weights)
    % 计算9节点四边形单元的一致质量矩阵
    M_el = zeros(18, 18);
    n_gauss = size(gauss_points, 1);
    
    for gp = 1:n_gauss
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        w = weights(gp);
        
        % 计算形函数
        [N, ~, ~] = calc_shape_functions(xi, eta);
        
        % 计算雅各比行列式
        [~, dN_dxi, dN_deta] = calc_shape_functions(xi, eta);
        dx_dxi = dN_dxi * el_nodes(:,1);
        dx_deta = dN_deta * el_nodes(:,1);
        dy_dxi = dN_dxi * el_nodes(:,2);
        dy_deta = dN_deta * el_nodes(:,2);
        det_J = det([dx_dxi, dx_deta; dy_dxi, dy_deta]);
        
        % 构建形函数矩阵（2×18）
        N_matrix = zeros(2, 18);
        for i = 1:9
            N_matrix(1, 2*i-1) = N(i);  % u方向
            N_matrix(2, 2*i)   = N(i);  % v方向
        end
        
        % 累加质量矩阵（ρ * N^T * N * dV）
        M_el = M_el + (N_matrix' * N_matrix) * rho * thickness * det_J * w;
    end
end