% 9节点四边形单元求解悬臂梁静力学问题
clear;
close all;
clc;

% 1. 问题参数
length = 0.1;         % 梁长(m)
width = 0.02;        % 梁高(m)
n_el_x = 10;         % x方向单元数
n_el_y = 2;         % y方向单元数
E = 210e9;          % 弹性模量(Pa)
nu = 0.3;           % 泊松比
load=-300/0.01;       % 自由端载荷(N/m)，厚度0.01m

% 2. 前处理：生成网格
[nodes, elements, material] = generate_rectangle_mesh(length, width, n_el_x, n_el_y, E, nu);
[gauss_points, weights] = get_gauss_points();
n_nodes = size(nodes, 1);
n_elements = size(elements, 1);

% 3. 单元计算：刚度矩阵和载荷向量
k_elements = cell(n_elements, 1);
f_elements = cell(n_elements, 1);
body_force = [0; 0];  % 无体力

for el = 1:n_elements
    el_nodes = nodes(elements(el, :), :);
    % 计算单元刚度矩阵
    k_el = calc_element_stiffness(el_nodes, material, gauss_points, weights);
    k_elements{el} = k_el;
    % 计算单元载荷向量（仅考虑自由端集中力，在主程序中处理）
    f_el = calc_element_load(el_nodes, material, body_force, gauss_points, weights);
    f_elements{el} = f_el;
end

% 4. 整体组装
[K, F] = assemble_global_matrix(elements, nodes, k_elements, f_elements);

% 5. 施加边界条件：悬臂梁左端固定，右端下边缘施加载荷
% 确认左端节点存在
left_nodes = find(nodes(:, 1) == 0);
if isempty(left_nodes)
    error('未找到左端节点！');
else
    % 施加固定边界条件
    for i = 1:max(size(left_nodes))
        node_id = left_nodes(i);
        % 获取节点对应的自由度编号
        ux_dof = 2 * node_id - 1;  % 假设Ux为节点的第一个自由度
        uy_dof = 2 * node_id;      % 假设Uy为节点的第二个自由度
        
        % 使用dofXX格式设置位移约束
        disp_bcs.(sprintf('dof%d', ux_dof)) = 0;  % 固定Ux
        disp_bcs.(sprintf('dof%d', uy_dof)) = 0;  % 固定Uy
    end
end

% 确认右端边缘节点存在
right_bottom_node = find(nodes(:, 1) == 0.1 & nodes(:, 2) == 0.01);
% 获取节点对应的自由度编号（假设Fy对应Uy，即第二个自由度）
fy_dof = 2 * right_bottom_node;
force_loads.(sprintf('dof%d', fy_dof)) = load;


% 处理边界条件
[K_mod, F_mod, active_dofs] = apply_boundary_conditions(K, F, disp_bcs, force_loads);

% 6. 求解
U = solve_system(K_mod, F_mod, active_dofs);

% 7. 后处理：计算应力应变并可视化
[element_stresses, element_strains] = calc_element_results(U, elements, nodes, material, gauss_points);
plot_displacement(nodes, U, elements, 50);  % 绘制变形图

% 输出最大位移（自由端中点）
free_end_node = find(nodes(:,1)>=length-1e-6 & nodes(:,2)>=width/2-1e-6, 1);
max_disp = U(2*free_end_node);
fprintf('自由端最大挠度: %.8f m\n', max_disp);