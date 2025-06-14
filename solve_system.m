function U = solve_system(K, F, active_dofs)
    % 输入：处理后的刚度矩阵，载荷向量，自由自由度索引
    % 输出：节点位移向量
    
    K_active = K(active_dofs, active_dofs);
    F_active = F(active_dofs);
    U_active = K_active \ F_active;  % MATLAB矩阵除法求解
    
    n_dofs = size(K, 1);
    U = zeros(n_dofs, 1);
    U(active_dofs) = U_active;
end
