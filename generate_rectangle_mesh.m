function [nodes, elements, material] = generate_rectangle_mesh(length, width, n_el_x, n_el_y, E, nu)
    % 输入：长度、宽度、x/y方向单元数、弹性模量、泊松比
    % 输出：节点坐标、单元连接表、材料属性
    
    % 计算节点数量（每个方向2n+1个节点，n为单元数）
    n_nodes_x = 2 * n_el_x + 1;
    n_nodes_y = 2 * n_el_y + 1;
    n_nodes = n_nodes_x * n_nodes_y;
    n_elements = n_el_x * n_el_y;
    
    % 初始化节点坐标（局部坐标到全局坐标映射）
    nodes = zeros(n_nodes, 2);
    for j = 0:n_nodes_y-1
        y = width * j / (n_nodes_y-1);
        for i = 0:n_nodes_x-1
            x = length * i / (n_nodes_x-1);
            node_idx = j * n_nodes_x + i + 1;
            nodes(node_idx, :) = [x, y];
        end
    end
    
    % 生成单元连接表（9节点/单元）
    elements = zeros(n_elements, 9);
    for el_y = 0:n_el_y-1
        for el_x = 0:n_el_x-1
            el_idx = el_x + el_y * n_el_x + 1;
            
            % 计算单元在节点数组中的基索引
            base_idx = el_y * 2 * n_nodes_x + el_x * 2 + 1;
            
            % 角节点（按逆时针顺序）
            elements(el_idx, 1) = base_idx;                         % 左下
            elements(el_idx, 2) = base_idx + 2;                     % 右下
            elements(el_idx, 3) = base_idx + 2 + 2 * n_nodes_x;     % 右上
            elements(el_idx, 4) = base_idx + 2 * n_nodes_x;         % 左上
            
            % 边中点
            elements(el_idx, 5) = base_idx + 1;                     % 下边中点
            elements(el_idx, 6) = base_idx + 2 + n_nodes_x;         % 右边中点
            elements(el_idx, 7) = base_idx + 2 * n_nodes_x + 1;     % 上边中点
            elements(el_idx, 8) = base_idx + n_nodes_x;             % 左边中点
            
            % 单元中心点
            elements(el_idx, 9) = base_idx + 1 + n_nodes_x;
        end
    end
    
    % 材料矩阵（平面应力）
    material.D = E/(1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    material.E = E;
    material.nu = nu;
    

        figure('Name', sprintf('矩形网格划分 (%dx%d个9节点四边形单元)', n_el_x, n_el_y), ...
               'Position', [100, 100, 800, 600]);
        hold on;
        grid on;
        axis equal;
        xlabel('X'); ylabel('Y');
        title(sprintf('矩形区域有限元网格划分 (%dx%d个单元)', n_el_x, n_el_y));
        
        % 绘制节点
        plot(nodes(:,1), nodes(:,2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        
        % 绘制单元边界
        for e = 1:n_elements
            % 提取当前单元的节点编号
            elem_nodes = elements(e, :);
            
            % 绘制四条边（包括中间节点）
            edges = [1 5 2; 2 6 3; 3 7 4; 4 8 1];
            for i = 1:4
                edge_nodes = edges(i, :);
                x = nodes(elem_nodes(edge_nodes), 1);
                y = nodes(elem_nodes(edge_nodes), 2);
                plot(x, y, 'b-', 'LineWidth', 1);
            end
            
            % 标记单元中心点
            center_node = elem_nodes(9);
            text(nodes(center_node,1), nodes(center_node,2), num2str(e), ...
                 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'r');
        end
        
        % 显示网格信息
        fprintf('矩形区域网格划分完成！\n');
        fprintf('总节点数: %d\n', n_nodes);
        fprintf('总单元数: %d\n', n_elements);
        fprintf('单元类型: 9节点四边形单元\n');
        fprintf('网格尺寸: %d×%d\n', n_el_x, n_el_y);
end