%计算单元应力应变
function [element_stresses, element_strains] = calc_element_results(U, elements, nodes, material, gauss_points)
    % 输入：节点位移，单元连接表，节点坐标，材料属性，高斯点
    % 输出：单元应力(每个单元1个积分点结果)，单元应变
    
    n_elements = size(elements, 1);
    n_gauss = size(gauss_points, 1);
    element_stresses = zeros(n_elements, 3);
    element_strains = zeros(n_elements, 3);
    [~, weights] = get_gauss_points();
    
    for elem = 1:n_elements
        elem_node = elements(elem, :);
        el_node_coords = nodes(elem_node, :);
        
        % 修正：正确提取单元节点的位移向量
        dof_indices = [];
        for i = 1:length(elem_node)
            dof_indices = [dof_indices, 2*elem_node(i)-1, 2*elem_node(i)];
        end
        el_disp = U(dof_indices);  % 确保el_disp是18×1的向量
        
        % 取中心高斯点计算结果（简化处理，实际可平均）
        gp_center = ceil(n_gauss/2);
        xi = gauss_points(gp_center, 1);
        eta = gauss_points(gp_center, 2);
        
        % 计算形函数偏导数
        [~, dN_dxi, dN_deta] = calc_shape_functions(xi, eta);
        
        % 计算雅各比矩阵
        dx_dxi = dN_dxi * el_node_coords(:, 1);
        dx_deta = dN_deta * el_node_coords(:, 1);
        dy_dxi = dN_dxi * el_node_coords(:, 2);
        dy_deta = dN_deta * el_node_coords(:, 2);
        J = [dx_dxi, dx_deta; dy_dxi, dy_deta];
        J_inv = inv(J);
        
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
        
        % 检查矩阵维度，避免错误
        if size(B, 2) ~= size(el_disp, 1)
            error(['单元 ', num2str(elem), '：B矩阵列数(', num2str(size(B, 2)), ...
                  ')与el_disp行数(', num2str(size(el_disp, 1)), ')不匹配']);
        end
        
        % 计算应变和应力
        strain = B * el_disp;
        stress = material.D * strain;
        
        % 修正：使用elem作为索引，而不是之前错误的el
        element_strains(elem, :) = strain';
        element_stresses(elem, :) = stress';
    end
end
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


