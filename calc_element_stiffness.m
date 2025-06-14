% 计算单元刚度矩阵
function k_el = calc_element_stiffness(el_nodes, material, gauss_points, weights)
    % 输入：单元节点坐标(9x2)，材料矩阵，高斯点，权重
    % 输出：单元刚度矩阵(18x18)
    k_el = zeros(18, 18);
    %获取高斯积分点数量
    n_gauss = size(gauss_points, 1);
    
    for gp = 1:n_gauss
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        w = weights(gp);
        
        % 计算形函数偏导数
        [~, dN_dxi, dN_deta] = calc_shape_functions(xi, eta);
        
        % 计算雅各比矩阵
        dx_dxi = dN_dxi * el_nodes(:, 1);
        dx_deta = dN_deta * el_nodes(:, 1);
        dy_dxi = dN_dxi * el_nodes(:, 2);
        dy_deta = dN_deta * el_nodes(:, 2);
        J = [dx_dxi, dx_deta; dy_dxi, dy_deta];
       
        % 计算雅各比矩阵的逆和行列式
        J_inv = inv(J);
        det_J = det(J);
       
        if det_J <= 0
            disp(J);  % 打印矩阵J的值
            error('单元畸变：雅各比行列式非正！');
        end
        
        % 计算应变-位移矩阵B
        dN_dx = J_inv(1,1)*dN_dxi + J_inv(1,2)*dN_deta;
        dN_dy = J_inv(2,1)*dN_dxi + J_inv(2,2)*dN_deta;
        B = zeros(3, 18);
        for i = 1:9
            B(1, 2*i-1) = dN_dx(i);
            B(2, 2*i) = dN_dy(i);
            B(3, 2*i-1) = dN_dy(i);
            B(3, 2*i) = dN_dx(i);
        end
        
        % 累加刚度矩阵（B^T*D*B*|J|*权重）
        k_el = k_el + B' * material.D * B * det_J * w;
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


