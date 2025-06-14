function plot_dynamic_results(nodes, elements, U_history, t, beam_length, beam_width, load_magnitude)
    % 绘制动态响应结果
    figure('Position', [100, 100, 800, 600]);
    
    % 1. 绘制自由端位移时程
    subplot(2,1,1);
    right_mid_node = find(abs(nodes(:,1)-beam_length) < 1e-6 & abs(nodes(:,2)-beam_width/2) < 1e-6);
    disp_dof = 2*right_mid_node; % y方向自由度
    plot(t, U_history(disp_dof,:), 'LineWidth', 1.5);
    xlabel('时间 (s)');
    ylabel('位移 (m)');
    title(['自由端中点位移时程 (载荷 = ', num2str(load_magnitude), ' N)']);
    grid on;
    
    % 2. 绘制变形动画
    subplot(2,1,2);
    axis equal;
    xlim([0, beam_length*1.2]);
    ylim([-beam_width*2, beam_width*2]);
    xlabel('X'); ylabel('Y');
    title('梁变形动画 (红色: 原始网格, 蓝色: 变形后)');
    hold on;
    
    % 标记载荷位置
    plot(nodes(right_mid_node,1), nodes(right_mid_node,2), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(nodes(right_mid_node,1), nodes(right_mid_node,2)-beam_width*0.5, ...
        ['F_y = ', num2str(load_magnitude), ' N'], 'Color', 'r', 'FontWeight', 'bold');
    
    % 绘制原始网格
    for el = 1:size(elements,1)
        el_nodes = elements(el,:);
        plot(nodes(el_nodes,1), nodes(el_nodes,2), 'r--', 'LineWidth', 0.5);
    end
    
    % 动画演示
    for i = 1:10:length(t)
        cla;
        % 绘制原始网格
        for el = 1:size(elements,1)
            el_nodes = elements(el,:);
            plot(nodes(el_nodes,1), nodes(el_nodes,2), 'r--', 'LineWidth', 0.5);
        end
        
        % 绘制变形网格
        U_matrix = reshape(U_history(:,i), 2, [])';
        deformed_nodes = nodes + 10 * U_matrix; % 位移放大10倍
        for el = 1:size(elements,1)
            el_nodes = elements(el,:);
            plot(deformed_nodes(el_nodes,1), deformed_nodes(el_nodes,2), 'b-', 'LineWidth', 1.5);
        end
        
        % 标记载荷点
        plot(nodes(right_mid_node,1), nodes(right_mid_node,2), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(nodes(right_mid_node,1), nodes(right_mid_node,2)-beam_width*0.5, ...
            ['F_y = ', num2str(load_magnitude), ' N'], 'Color', 'r', 'FontWeight', 'bold');
        
        title(['变形动画 (t = ', num2str(t(i), '%.3f'), ' s)']);
        drawnow;
    end
    hold off;
end