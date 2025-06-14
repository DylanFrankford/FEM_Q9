% 计算单元载荷向量（以体力为例）
function f_el = calc_element_load(el_nodes, material, body_force, gauss_points, weights)
    % 输入：单元节点坐标，材料属性，体力向量[fx; fy]，高斯点，权重
    % 输出：单元载荷向量(18x1)
    
    f_el = zeros(18, 1);
    n_gauss = size(gauss_points, 1);
    
    for gp = 1:n_gauss
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        w = weights(gp);
        
        % 计算形函数
        [N, ~, ~] = calc_shape_functions(xi, eta);
        
        % 构建形函数矩阵
        N_matrix = zeros(18, 1);
        for i = 1:9
            N_matrix(2*i-1) = N(i);
            N_matrix(2*i) = N(i);
        end
        
        % 计算雅各比行列式
        [~, dN_dxi, dN_deta] = calc_shape_functions(xi, eta);
        dx_dxi = dN_dxi * el_nodes(:, 1);
        dx_deta = dN_deta * el_nodes(:, 1);
        dy_dxi = dN_dxi * el_nodes(:, 2);
        dy_deta = dN_deta * el_nodes(:, 2);
        det_J = det([dx_dxi, dx_deta; dy_dxi, dy_deta]);
        
        % 累加载荷向量（N^T*body_force*|J|*权重）
        f_el = f_el + N_matrix * body_force' * det_J * w;
    end
end

% 计算9节点单元的形函数及其偏导数
function [N, dN_dxi, dN_deta] = calc_shape_functions(xi, eta)
    % 输入：局部坐标(xi, eta)
    % 输出：形函数N(1x9)，偏导数dN_dxi(1x9), dN_deta(1x9)
    
    % 角节点形函数及偏导（1-4号）
N(3) = 0.25 * (xi + xi^2) * (eta + eta^2);
N(4) = 0.25 * (-xi + xi^2) * (eta + eta^2);
N(1) = 0.25 * (-xi + xi^2) * (-eta + eta^2);
N(2) = 0.25 * (xi + xi^2) * (-eta + eta^2);

% 对xi求导（计算A'(xi)·B(eta)）
dN_dxi(3) = 0.25 * (1 + 2*xi) * (eta + eta^2);  % d(xi+xi²)/dxi = 1+2xi
dN_dxi(4) = 0.25 * (-1 + 2*xi) * (eta + eta^2); % d(-xi+xi²)/dxi = -1+2xi
dN_dxi(1) = 0.25 * (-1 + 2*xi) * (-eta + eta^2);
dN_dxi(2) = 0.25 * (1 + 2*xi) * (-eta + eta^2);

% 对eta求导（计算A(xi)·B'(eta)）
dN_deta(3) = 0.25 * (xi + xi^2) * (1 + 2*eta);  % d(eta+eta²)/deta = 1+2eta
dN_deta(4) = 0.25 * (-xi + xi^2) * (1 + 2*eta);
dN_deta(1) = 0.25 * (-xi + xi^2) * (-1 + 2*eta); % d(-eta+eta²)/deta = -1+2eta
dN_deta(2) = 0.25 * (xi + xi^2) * (-1 + 2*eta);
    
% 边中间节点形函数及偏导（5-8号）
N(5) = 0.5 * (1 - xi^2) * (eta^2 - eta);
N(6) = 0.5 * (1 - eta^2) * (xi^2 + xi);
N(7) = 0.5 * (1 - xi^2) * (eta^2 + eta);
N(8) = 0.5 * (xi^2 - xi) * (1 - eta^2);

% 对xi求导（提取关于xi的导数）
dN_dxi(5) = 0.5 * (-2*xi) * (eta^2 - eta);       % d(1-xi²)/dxi = -2xi
dN_dxi(6) = 0.5 * (2*xi + 1) * (1 - eta^2);     % d(xi²+xi)/dxi = 2xi+1
dN_dxi(7) = 0.5 * (-2*xi) * (eta^2 + eta);       % d(1-xi²)/dxi = -2xi
dN_dxi(8) = 0.5 * (2*xi - 1) * (1 - eta^2);     % d(xi²-xi)/dxi = 2xi-1

% 对eta求导（提取关于eta的导数）
dN_deta(5) = 0.5 * (1 - xi^2) * (2*eta - 1);     % d(eta²-eta)/deta = 2eta-1
dN_deta(6) = 0.5 * (-2*eta) * (xi^2 + xi);       % d(1-eta²)/deta = -2eta
dN_deta(7) = 0.5 * (1 - xi^2) * (2*eta + 1);     % d(eta²+eta)/deta = 2eta+1
dN_deta(8) = 0.5 * (xi^2 - xi) * (-2*eta);       % d(1-eta²)/deta = -2eta

    % 中心节点形函数及偏导（9号）
    N(9) = (1-xi^2) * (1-eta^2);
    dN_dxi(9) = -2*xi*(1-eta^2);
    dN_deta(9) = -2*eta*(1-xi^2);
end