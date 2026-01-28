%A Linear Kalman Filter for MARG Orientation Estimation 
% Using the Algebraic Quaternion Algorithm
function [q_k, att_k, q_meas, att_meas, N_k] = UFIR_Delete(gyr, acc, mag, cov_g, cov_a, cov_m, n, t, N_max)
%定义L为为载体（传感器）坐标系，G为全局（地球）坐标系
%定义各测量向量和参考向量，均需归一化处理

%N_max = 50;
N_min = 5;
cov_acc = [cov_a^2    0       0;
              0    cov_a^2    0;
              0       0    cov_a^2] / (9.81^2);
cov_mag = [cov_m^2    0       0;
              0    cov_m^2    0;
              0       0    cov_m^2] / (0.52^2);
cov_u = [cov_acc     zeros(3, 3);
         zeros(3, 3)    cov_mag];
w = zeros(1, 3);

g_G = [0; 0; 9.81];

% h_G = [0.5547; 0; 0.8321];
% h_G = [-0.6036; 0.2182; -1.0868];
% h_G = [0.6963; 0; 0.7177]
h_G = [0.35; 0; 0.35];

q_meas = zeros(4, n);
q_k = zeros(4, n);
att_k = zeros(3, n);
att_meas = zeros(3, n);
q_k(:, 1) = [1; 0; 0; 0];


N_k = zeros(n,2);
% N = 12;
A_n = zeros(4, 4, N_max);
theta1_max_deg = 20;
theta2_max_deg = 20;
N1 = N_min;
N2 = N_min;
N_small = N_min;
N0 = N_min;
for K = 1:n-N0
    i = 1;
    if K <= N_max
        N = N_min;
    end
    for k = N0 + K - N : N0 + K - 1
        
        acc(k, :) = normalize_vector(acc(k, :));
        mag(k, :) = normalize_vector(mag(k, :));
        lambd_1 = sqrt((acc(k, 3) + 1) / 2);
        lambd_2 = sqrt((1 - acc(k, 3)) / 2);
        if acc(k, 3) >= 0
            q_acc = [lambd_1; -acc(k, 2) / 2 / lambd_1; acc(k, 1) / 2 / lambd_1; 0 ];
        else
            q_acc = [-acc(k, 2) / 2 / lambd_2; lambd_2; 0; acc(k, 1) / 2 / lambd_2 ];
        end
        l = (quaternion_to_rotation_matrix(q_acc)') * mag(k, :)';
        gamma = l(1)^2 + l(2)^2;
    
        if l(1) >= 0
            q_mag = [sqrt(gamma+l(1)*sqrt(gamma))/sqrt(2*gamma); 0; 0; l(2)/sqrt(2*(gamma+l(1)*sqrt(gamma)))];
        else
            q_mag = [l(2)/sqrt(2*(gamma-l(1)*sqrt(gamma))); 0; 0; sqrt(gamma-l(1)*sqrt(gamma))/sqrt(2*gamma)];
        end
        q_meas(:, k) = quatmultiply(q_acc', q_mag');
        % % % J = AQUA_J(acc(k, :), mag(k,:), l, q_acc, q_mag);
        % % % R = J * cov_u * J';
        % % % cov_k = [-q_k(2, k)   -q_k(3, k)   -q_k(4, k);
        % % %           q_k(1, k)   -q_k(4, k)    q_k(3, k);
        % % %           q_k(4, k)    q_k(1, k)   -q_k(2, k);
        % % %          -q_k(3, k)    q_k(2, k)    q_k(1, k)];
        % % % Q = 0.25 * t^2 * cov_k * cov_g * cov_k';
        w = gyr(k, :);
        omega_w = [ 0  -w(1) -w(2) -w(3);
                  w(1)    0   w(3) -w(2);
                  w(2) -w(3)    0   w(1);
                  w(3)  w(2) -w(1)    0];
        A = eye(4) + 0.5 * omega_w * t;
        A_n(:, :, i) = A;
        i = i + 1;
        if (q_meas(:, k)' * q_k(:, k)) < 0
            q_meas(:, k) = -q_meas(:, k);
        else
            q_meas(:, k) = q_meas(:, k);
        end
    end

    q_k(:, N0+K) = AdaptiveU_FIR_Filter(A_n, 0, eye(4), eye(4), N, N_small, q_meas(:, N0 + K - N : N0 + K - 1), 0, N0+K);
    q_k(:, N0+K) = normalize_vector(q_k(:, N0+K));
    % % [qnb, att_k(:, N+K), Cnb] = attsyn(q_k(:, N+K));
    att_k(:, N0+K) = q_att(q_k(:, N0+K), 1);
    
    % 计算当前时刻与前一时刻的四元数内积
    %归一化角度变化
    delta_pitch = att_k(1, N0+K) - att_k(1, N0+K-1);
    delta_pitch = delta_pitch * 180 / pi;
    pitch_ratio = min(delta_pitch / theta1_max_deg, 1);

    delta_yaw   = att_k(3, N0+K) - att_k(3, N0+K-1);
    delta_yaw = delta_yaw * 180 / pi;
    yaw_ratio   = min(delta_yaw   / theta2_max_deg, 1);


    % 基于theta动态调整N
    N1 = N_min + (N_max - N_min) * pitch_ratio;
    N1 = round(max(N_min, min(N_max, N1))); % 取整并限制在[min, max]
    N1 = 5;
    N_k(N_min+K, 1) = N1;

    % 基于theta动态调整N
    N2 = N_max - (N_max - N_min) * yaw_ratio;
    N2 = round(max(N_min, min(N_max, N2))); % 取整并限制在[min, max]
    N2 = N_max;
    N_k(N_min+K, 2) = N2;

    N = max(N1, N2);
    N_small = min(N1, N2);
    A_n = zeros(4, 4, N);

end

end
    


