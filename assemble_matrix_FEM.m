function [K, M] = assemble_matrix_FEM(elements, nodes, K_elements, M_elements)
    % 组装全局刚度矩阵和质量矩阵
    n_nodes = size(nodes, 1);
    n_dofs = 2 * n_nodes;
    n_elements = size(elements, 1);
    
    K = sparse(n_dofs, n_dofs);
    M = sparse(n_dofs, n_dofs);
    
    for el = 1:n_elements
        el_nodes = elements(el, :);
        dof_indices = zeros(18, 1);
        for i = 1:9
            dof_indices(2*i-1) = 2*el_nodes(i) - 1; % x自由度
            dof_indices(2*i)   = 2*el_nodes(i);     % y自由度
        end
        
        % 组装刚度矩阵
        K_el = K_elements{el};
        for i = 1:18
            for j = 1:18
                K(dof_indices(i), dof_indices(j)) = K(dof_indices(i), dof_indices(j)) + K_el(i,j);
            end
        end
        
        % 组装质量矩阵
        M_el = M_elements{el};
        for i = 1:18
            for j = 1:18
                M(dof_indices(i), dof_indices(j)) = M(dof_indices(i), dof_indices(j)) + M_el(i,j);
            end
        end
    end
end
