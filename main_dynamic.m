% main_dynamic.m - 主程序文件（动态分析）
clear; close all; clc;

% 1. 问题参数
beam_length = 0.1;       % 梁长(m)
beam_width = 0.02;       % 梁高(m)
n_el_x = 10;            % x方向单元数
n_el_y = 2;             % y方向单元数
E = 210e9;              % 弹性模量(Pa)
nu = 0.3;               % 泊松比
rho = 7800;             % 材料密度(kg/m^3)
thickness = 0.01;       % 厚度(m)

% 2. 时间参数
T_total = 0.1;          % 总时间(s)
dt = 1e-4;              % 时间步长(s)
n_steps = T_total/dt;

% 3. Newmark参数
gamma = 0.5;            % Newmark参数
beta_nm = 0.25;         % Newmark参数

% 4. 前处理：生成网格
[nodes, elements, material] = generate_rectangle_mesh(beam_length, beam_width, n_el_x, n_el_y, E, nu);
[gauss_points, weights] = get_gauss_points();
n_nodes = size(nodes, 1);
n_elements = size(elements, 1);
n_dofs = 2 * n_nodes;

% 5. 单元计算：刚度矩阵和质量矩阵
K_elements = cell(n_elements, 1);
M_elements = cell(n_elements, 1);
body_force = [0; 0];    % 无体力

for el = 1:n_elements
    el_nodes = nodes(elements(el, :), :);
    % 计算单元刚度矩阵
    K_elements{el} = calc_element_stiffness(el_nodes, material, gauss_points, weights);
    % 计算单元质量矩阵
    M_elements{el} = calc_element_mass(el_nodes, rho, thickness, gauss_points, weights);
end

% 6. 整体组装
[K, M] = assemble_matrix_FEM(elements, nodes, K_elements, M_elements);

% 7. 施加边界条件：悬臂梁左端固定
left_nodes = find(abs(nodes(:,1)) < 1e-6); % x=0的节点
fixed_dofs = [];
for i = 1:length(left_nodes)
    fixed_dofs = [fixed_dofs, 2*left_nodes(i)-1, 2*left_nodes(i)]; % 固定x和y方向
end
all_dofs = 1:n_dofs;
free_dofs = setdiff(all_dofs, fixed_dofs);

% 8. 计算前两阶固有频率和Rayleigh阻尼
[V, D] = eigs(K(free_dofs, free_dofs), M(free_dofs, free_dofs), 2, 'smallestabs');
omega = sqrt(diag(D));  % 固有频率(rad/s)
freq_hz = omega/(2*pi); % 转换为Hz

% 设置0.2%阻尼比
xi = 0.002; % 0.2%阻尼比

% 计算Rayleigh阻尼系数
A = [1/(2*omega(1)), omega(1)/2;
     1/(2*omega(2)), omega(2)/2];
b = [xi; xi];
coeffs = A \ b;
alpha = coeffs(1);
beta = coeffs(2);

fprintf('前两阶固有频率: %.2f Hz 和 %.2f Hz\n', freq_hz(1), freq_hz(2));
fprintf('计算得到的Rayleigh阻尼系数: alpha=%.3e, beta=%.3e\n', alpha, beta);

% 计算阻尼矩阵
C = alpha * M + beta * K;

% 9. 验证实际阻尼比
for k = 1:2
    phi = V(:,k);
    modal_C = phi' * C(free_dofs, free_dofs) * phi;
    modal_M = phi' * M(free_dofs, free_dofs) * phi;
    modal_K = phi' * K(free_dofs, free_dofs) * phi;
    xi_actual = modal_C / (2 * sqrt(modal_K * modal_M));
    fprintf('模态%d: 理论阻尼比=%.4f%%, 实际阻尼比=%.4f%%\n', ...
            k, xi*100, xi_actual*100);
end

% 10. 初始化载荷条件（右端中点施加动态载荷）
right_mid_node = find(abs(nodes(:,1)-beam_length) < 1e-6 & abs(nodes(:,2)-beam_width/2) < 1e-6);
F = zeros(n_dofs, 1);
load_dof = 2*right_mid_node; % y方向自由度
load_magnitude = -300 / thickness; % 转换为集中力 (N)

% 11. 初始条件
U0 = zeros(n_dofs, 1);   % 初始位移
V0 = zeros(n_dofs, 1);   % 初始速度
A0 = zeros(n_dofs, 1);   % 初始加速度

% 12. 时间积分 (Newmark方法)
[U_history, t] = newmark_solver(M, C, K, F, load_dof, load_magnitude, U0, V0, A0, dt, T_total, ...
                               gamma, beta_nm, free_dofs);

% 13. 后处理
plot_dynamic_results(nodes, elements, U_history, t, beam_length, beam_width, load_magnitude);

