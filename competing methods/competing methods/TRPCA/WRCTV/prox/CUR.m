function [X, U, V, Z] = CUR(M, r, d, num_observed)
    % CUR+ 算法实现（部分观测矩阵的低秩逼近）
    % 输入：
    %   M: 目标矩阵（n×m，可能部分观测）
    %   r: 目标秩
    %   d: 采样的行/列数（d >= r）
    %   num_observed: 随机观测的条目数（|Ω|）
    % 输出：
    %   X: 估计的低秩矩阵
    %   U: 左奇异向量（n×r）
    %   V: 右奇异向量（m×r）
    %   Z: 回归系数矩阵（r×r）

    [n, m] = size(M);
    
    %% 1. 随机采样d行和d列
    col_indices = randperm(m, d);
    row_indices = randperm(n, d);
    A = M(:, col_indices); % 随机列（n×d）
    B = M(row_indices, :); % 随机行（d×m）
    
    %% 2. 计算近似奇异向量 U 和 V
    [U, ~] = eigs(A * A', r);       % 前r个左奇异向量（n×r）
    [V, ~] = eigs(B' * B, r);         % 前r个右奇异向量（m×r）
    
    %% 3. 随机采样观测条目 Ω
    Omega = randperm(n * m, num_observed);
    [rows, cols] = ind2sub([n, m], Omega);
    M_Omega = M(Omega); % 观测值
    
    %% 4. 优化回归问题：min_Z ||R_Ω(M) - R_Ω(U_hat Z V_hat')||_F^2
    % 构造设计矩阵（每个观测条目对应一个线性方程）
    Z = zeros(r, r);         % 待优化变量 (r x r)
    Z_prev = Z;
    t = 1;                  % Nesterov动量参数
    L = 1;                  % Lipschitz常数估计 (需根据问题调整)
    max_iter = 100;
    tol = 1e-5;
    
    for iter = 1:max_iter
        % 计算当前预测矩阵
        M_pred = U * Z * V';
        
        % 计算梯度 (仅在Omega处计算误差)
        grad = zeros(r, r);
        for k = 1:length(Omega)
            [i, j] = ind2sub([n, m], Omega(k));
            grad = grad + (M_pred(i,j) - M_Omega(k)) * (U(i,:)' * V(j,:));
        end
        grad = grad / length(Omega);
        
        % Nesterov加速步骤
        t_next = (1 + sqrt(1 + 4 * t^2)) / 2;
        Z_next = Z - (1/L) * grad;
        Z = Z_next + ((t - 1)/t_next) * (Z_next - Z_prev);
        
        % 更新变量
        Z_prev = Z;
        t = t_next;
        
        % 检查收敛性
        if iter > 1 && norm(grad, 'fro') < tol
            break;
        end
    end
    
    %% 5. 重构低秩矩阵
    X = U * Z * V';
end