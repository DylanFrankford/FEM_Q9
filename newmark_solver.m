function [U_history, t] = newmark_solver(M, C, K, F, load_dof, load_magnitude, U0, V0, A0, dt, T_total, gamma, beta_nm, free_dofs)
    % Newmark时间积分求解器
    n_steps = ceil(T_total/dt);
    t = linspace(0, T_total, n_steps+1)';
    n_dofs = size(M, 1);
    U_history = zeros(n_dofs, n_steps+1);
    
    % 初始条件
    U = U0;
    V = V0;
    A = A0;
    U_history(:,1) = U;
    
    % 计算有效刚度矩阵
    K_eff = K(free_dofs, free_dofs) + (gamma/(beta_nm*dt))*C(free_dofs, free_dofs) + ...
            (1/(beta_nm*dt^2))*M(free_dofs, free_dofs);
    
    % 时间步进
    for i = 1:n_steps
        % 施加阶跃载荷
        F(load_dof) = load_magnitude;  % 恒定载荷
        
        % 计算预测值
        U_pred = U + dt*V + (0.5-beta_nm)*dt^2*A;
        V_pred = V + (1-gamma)*dt*A;
        
        % 计算有效载荷
        F_eff = F(free_dofs) - C(free_dofs, free_dofs)*V_pred(free_dofs) - ...
                K(free_dofs, free_dofs)*U_pred(free_dofs);
        F_eff = F_eff + M(free_dofs, free_dofs)*(U_pred(free_dofs)/(beta_nm*dt^2) + ...
                V_pred(free_dofs)/(beta_nm*dt));
        
        % 求解位移增量
        delta_U = K_eff \ F_eff;
        
        % 更新位移、速度和加速度
        U(free_dofs) = U_pred(free_dofs) + delta_U;
        V(free_dofs) = V_pred(free_dofs) + gamma/(beta_nm*dt)*delta_U;
        A(free_dofs) = A(free_dofs) + 1/(beta_nm*dt^2)*delta_U;
        
        % 存储结果
        U_history(:,i+1) = U;
    end
end
