function plot_displacement(nodes, U, elements, scale_factor)
    if nargin < 4
        scale_factor = 1;  % 位移缩放因子，便于可视化
    end
    
    figure;
    hold on;
    
    % 确保U是n_nodes×2的矩阵或可重塑为该格式
    if size(U, 2) == 1  % 如果U是列向量（2n_nodes×1）
        U_matrix = reshape(U, 2, []).';  % 转换为n_nodes×2的矩阵
    else
        U_matrix = U;  % 假设U already是n_nodes×2的矩阵
    end
    
    % 计算变形后节点坐标
    deformed_nodes = nodes + scale_factor * U_matrix;
    
    % 绘制变形后网格
    for el = 1:size(elements, 1)
        el_nodes = elements(el, :);
        plot(deformed_nodes(el_nodes, 1), deformed_nodes(el_nodes, 2), 'b-', 'LineWidth', 0.1);
    end
    
    % 绘制原始网格（半透明红色虚线）
    for el = 1:size(elements, 1)
        el_nodes = elements(el, :);
        plot(nodes(el_nodes, 1), nodes(el_nodes, 2), ...
             'Color', [1 0 0 0.3], 'LineStyle', '--', 'LineWidth', 0.5);
    end
    
    axis equal;
    title('变形图（红色为原始网格，蓝色为变形后）');
    xlabel('X'); ylabel('Y');
    hold off;
end