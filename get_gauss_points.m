% 定义高斯积分点（3x3）
function [gauss_points, weights] = get_gauss_points()
    % 3x3高斯积分点和权重
    xi = [-sqrt(3/5), 0, sqrt(3/5)];
    eta = xi;
    w = [5/9, 8/9, 5/9];
    
    gauss_points = zeros(9, 2);
    weights = zeros(9, 1);
    idx = 0;
    for i = 1:3
        for j = 1:3
            idx = idx + 1;
            gauss_points(idx, :) = [xi(i), eta(j)];
            weights(idx) = w(i) * w(j);
        end
    end
end
% 计算9节点单元的形函数及其偏导数
function [N, dN_dxi, dN_deta] = calc_shape_functions(xi, eta)
    % 输入：局部坐标(xi, eta)
    % 输出：形函数N(1x9)，偏导数dN_dxi(1x9), dN_deta(1x9)
    
    N = zeros(1, 9);
    dN_dxi = zeros(1, 9);
    dN_deta = zeros(1, 9);
    
    % 角节点形函数及偏导（1-4号）
    xi_coords = [-1, 1, 1, -1];
    eta_coords = [1, 1, -1, -1];
    for i = 1:4
        xi_i = xi_coords(i);
        eta_i = eta_coords(i);% 角节点局部坐标
        N(i) = 1/4 * (1+xi_i*xi) * (1+eta_i*eta) * (xi_i*xi + eta_i*eta - 1);
        dN_dxi(i) = 1/4 * xi_i * (1+eta_i*eta) * (2*xi_i*xi + eta_i*eta);
        dN_deta(i) = 1/4 * eta_i * (1+xi_i*xi) * (xi_i*xi + 2*eta_i*eta);
    end
    
    % 边中间节点形函数及偏导（5-8号，5/6为xi方向边，7/8为eta方向边）
    % 边5（xi=-1, eta方向）
    N(5) = 1/2 * (1-xi^2) * (1+eta);
    dN_dxi(5) = -xi * (1+eta);
    dN_deta(5) = 1/2 * (1-xi^2);
    % 边6（xi=1, eta方向）
    N(6) = 1/2 * (1-xi^2) * (1-eta);
    dN_dxi(6) = -xi * (1-eta);
    dN_deta(6) = -1/2 * (1-xi^2);
    % 边7（eta=1, xi方向）
    N(7) = 1/2 * (1-eta^2) * (1+xi);
    dN_dxi(7) = 1/2 * (1-eta^2);
    dN_deta(7) = -eta * (1+xi);
    % 边8（eta=-1, xi方向）
    N(8) = 1/2 * (1-eta^2) * (1-xi);
    dN_dxi(8) = -1/2 * (1-eta^2);
    dN_deta(8) = -eta * (1-xi);
    
    % 中心节点形函数及偏导（9号）
    N(9) = (1-xi^2) * (1-eta^2);
    dN_dxi(9) = -2*xi*(1-eta^2);
    dN_deta(9) = -2*eta*(1-xi^2);
end


