%卡尔曼滤波器
function [x_est,P_est] = kalman_filter(A, B, C, Q, R, measurements, x_init, P_init)
    % A: 状态转移矩阵
    % B: 控制输入矩阵
    % C: 观测矩阵
    % Q: 过程噪声协方差矩阵
    % R: 测量噪声方差
    % measurements: 观测值序列
    % x_init: 初始状态估计
    % P_init: 初始协方差估计

    %状态和协方差初始化
    x_est = x_init;
    P_est = P_init;

    %卡尔曼滤波
    % for t = 1:length(measurements)
        %预测
        x_pred = A * x_est;
        P_pred = A * P_est * A' + Q;
        % % P_pred = 0.64 * P_pred;
        %更新
        K = P_pred * C' / (C * P_pred * C' + R);
        x_est = x_pred + K * (measurements - C * x_pred);
        P_est = (eye(size(A)) - K * C) * P_pred;
    % end
end