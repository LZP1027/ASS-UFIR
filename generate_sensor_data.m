function [acc, mag] = generate_sensor_data(q_k, g_G, h_G, noise_std)
    % 计算传感器测量值
    % q_k: 当前姿态四元数
    % g_G: 全局坐标系中的重力向量
    % h_G: 全局坐标系中的地磁场向量
    % noise_std: 噪声标准差

    % 根据当前姿态计算加速度计测量值
    % 将全局坐标系中的重力向量转换到载体坐标系中
    q_k = q_k';
    g_L = quatrotate(quatconj(q_k), g_G');
    % 添加噪声
    % acc = g_L;
    acc = g_L + noise_std * randn(1, 3);
    
    % 根据当前姿态计算磁力计测量值
    % 将全局坐标系中的地磁场向量转换到载体坐标系中
    h_L = quatrotate(quatconj(q_k), h_G');
    % 添加噪声
    % mag = h_L;
    mag = h_L + noise_std * randn(1, 3);
end
