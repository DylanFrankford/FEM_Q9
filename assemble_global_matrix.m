function [K, F] = assemble_global_matrix(elements, nodes, k_elements, f_elements)
    % 输入：单元连接表，节点坐标，单元刚度矩阵列表，单元载荷向量列表
    % 输出：总体刚度矩阵(稀疏)，总体载荷向量
    
    n_nodes = size(nodes, 1);
    n_dofs = 2 * n_nodes;
    n_elements = size(elements, 1);
    
    % 初始化稀疏刚度矩阵和载荷向量
    K = sparse(n_dofs, n_dofs);
    F = zeros(n_dofs, 1);
    
    for el = 1:n_elements
        el_nodes = elements(el, :);
        dof_indices = zeros(18, 1);
        for i = 1:9
            dof_indices(2*i-1) = 2*el_nodes(i) - 1;  % x自由度
            dof_indices(2*i) = 2*el_nodes(i);        % y自由度
        end
        
        % 组装刚度矩阵
        k_el = k_elements{el};
        for i = 1:18
            for j = 1:18
                K(dof_indices(i), dof_indices(j)) = K(dof_indices(i), dof_indices(j)) + k_el(i, j);
            end
        end
        
        % 组装载荷向量
        f_el = f_elements{el};
        for i = 1:18
            F(dof_indices(i)) = F(dof_indices(i)) + f_el(i);
        end
    end
end