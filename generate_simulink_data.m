function [q_k, att, gyr, acc, mag] = generate_simulink_data(n, omega_x, omega_y, omega_z) 
% 设置仿真参数
% n:数据点数量
t = 0.01; % 时间步长
cov_g = 0.0012;
% 假设重力和地磁场向量在全局坐标系中为 [0; 0; 1]
% g_G = [0; 0; 1]; % 全局坐标系中的重力向量
% h_G = [0.6963; 0; 0.7177]; % 全局坐标系中的地磁场向量

g_G = [0; 0; 9.81]; % 全局坐标系中的重力向量
h_G = [0.35; 0; 0.35]; % 全局坐标系中的地磁场向量

% 生成姿态角速度数据（陀螺仪数据）
gyr = zeros(n, 3); % 初始化
att = zeros(3, n);
for i = 1:n       %%%test2

    omega_x = 0.01 * randi([-9,9],1,1) + 0.001 * randn; % 角速度（360度/s）
    omega_y = 0.01 * randi([-9,9],1,1) + 0.001 * randn;
    omega_z = 1 * randi([-9,9],1,1) + 0.001 * randn;


    gyr(i, 1) = omega_x; % 绕x轴旋转
    gyr(i, 2) = omega_y; % 绕y轴旋转
    gyr(i, 3) = omega_z; % 绕z轴旋转
end


% % % % for i = 1:n   %%%%%test1
% % % % 
% % % %     omega_x = 0.01 * randi([-9,9],1,1) + 0.001 * randn; % 角速度（360度/s）
% % % %     omega_y = 0.01 * randi([-9,9],1,1) + 0.001 * randn;
% % % %     omega_z = 0.01 * randi([-9,9],1,1) + 0.001 * randn;
% % % % 
% % % % 
% % % % 
% % % %     gyr(i, 1) = omega_x; % 绕x轴旋转
% % % %     gyr(i, 2) = omega_y; % 绕y轴旋转
% % % %     gyr(i, 3) = omega_z; % 绕z轴旋转
% % % % end


% 生成加速度计和磁力计数据
acc = zeros(n, 3); % 初始化加速度计数据
mag = zeros(n, 3); % 初始化磁力计数据

% 初始化姿态
q_k = zeros(4, n); % 四元数数组
q_k(:, 1) = [1; 0; 0; 0]; % 初始姿态为单位四元数

% 生成仿真数据
for k = 1:n-1
    % 根据陀螺仪数据更新姿态
    w = gyr(k, :);
    omega_w = [ 0  -w(1) -w(2) -w(3);
              w(1)    0   w(3) -w(2);
              w(2) -w(3)    0   w(1);
              w(3)  w(2) -w(1)    0];
    % % cov_k = [-q_k(2, k)   -q_k(3, k)   -q_k(4, k);
    % %           q_k(1, k)   -q_k(4, k)    q_k(3, k);
    % %           q_k(4, k)    q_k(1, k)   -q_k(2, k);
    % %          -q_k(3, k)    q_k(2, k)    q_k(1, k)];
    % % Q = 0.25 * t^2 * cov_k * cov_g * cov_k';
    % % % w_k = gauss_rnd(zeros(4,1),Q,4);
    q_k(:, k+1) = (eye(4) + 0.5 * omega_w * t) * q_k(:, k) + 0.001 * randn(4, 1);
    q_k(:, k+1) = normalize_vector(q_k(:, k+1));
    % % [qnb, att(:, k+1), Cnb] = attsyn(q_k(:, k+1));
    att(:, k+1) = q_att(q_k(:, k+1), 1);
    % 根据姿态生成加速度计和磁力计数据（示例中未添加噪声）
    % 生成传感器数据
    [acc(k, :), mag(k, :)] = generate_sensor_data(q_k(:, k), g_G, h_G, 0.01);


    % % if (k>=200 && k<=300)
    % %    acc(k, :) = acc(k, :) + 0.5 * randn(1, 3);   
    % % end
    % % if (k>=450 && k<=550)
    % %    acc(k, :) = acc(k, :) + randi([-1,1], 1, 1) * randn(1, 3);
    % % end
    % % if (k>=700 && k<=800)
    % %    acc(k, :) = acc(k, :) + 1.5 * randn(1, 3);
    % % end



    % %  if (k>=100 && k<=200)
    % %    mag(k, :) = mag(k, :) + 0.5*randn(1, 3);   
    % % end
    % if k == 20 || k == 100 || k == 1000
    %    acc(k, :) = acc(k, :) + randn(1, 3);   
    % end
end
end
