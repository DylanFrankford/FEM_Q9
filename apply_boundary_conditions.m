function [K_mod, F_mod, active_dofs] = apply_boundary_conditions(K, F, disp_bcs, force_loads, large_num)
    % 输入：总体刚度矩阵，载荷向量，位移约束，力载荷，大数因子
    % 输出：处理后的刚度矩阵、载荷向量，自由自由度索引
    
    if nargin < 5
        large_num = 1e15;  % 乘大数法的默认因子
    end
    
    n_dofs = size(K, 1);
    fixed_dofs = false(n_dofs, 1);
    active_dofs = 1:n_dofs;
    
   % 处理位移约束（乘大数法）
for bc = fieldnames(disp_bcs)'
    field_name = bc{1};
    % 提取字段名中的数字部分（假设格式为'dofXXX'）
    dof_str = regexp(field_name, 'dof(\d+)', 'tokens', 'once');
    if ~isempty(dof_str)
        dof = str2double(dof_str{1}); % 转换为数值
        val = disp_bcs.(field_name);
        fixed_dofs(dof) = true;
        K(dof, dof) = K(dof, dof) * large_num;
        F(dof) = K(dof, dof) * val;
    end
end
    for load = fieldnames(force_loads)'
    field_name = load{1};
    % 提取dof后面的所有数字字符
    dof_str = regexp(field_name, 'dof(\d+)', 'tokens', 'once');
        if ~isempty(dof_str)
            dof = str2double(dof_str{1}); % 转换为数值
            val = force_loads.(field_name);
            F(dof) = F(dof) + val;
        end
    end
    
    % 提取自由自由度
    active_dofs = find(~fixed_dofs);
    
    % 将修改后的矩阵赋值给输出参数
    K_mod = K;
    F_mod = F;
end