 %A Linear Kalman Filter for MARG Orientation Estimation 
% Using the Algebraic Quaternion Algorithm
function [q_k, att_k, q_meas, att_meas, N_k] = UFIR_Separation(gyr, acc, mag, cov_g, cov_a, cov_m, n, t, N_max)
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
q_k1 = zeros(4, n);
q_k2 = zeros(4, n);
att_k = zeros(3, n);
att_meas = zeros(3, n);
q_k1(:, 1) = [1; 0; 0; 0];
q_k2(:, 1) = [1; 0; 0; 0];

N_k = zeros(n,2);
% N = 12;
A_n = zeros(4, 4, N_max);
theta1_max_deg = 20;
theta2_max_deg = 60;
N1 = N_min;
N2 = N_min;
N0 = N_min;
for K = 1:n-N0
    i = 1;
    if K <= N_max
        N1 = N_min;
    end
    for k = N0 + K - N1 : N0 + K - 1
        
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
        J = AQUA_J(acc(k, :), mag(k,:), l, q_acc, q_mag);
        R = J * cov_u * J';
        cov_k = [-q_k1(2, k)   -q_k1(3, k)   -q_k1(4, k);
                  q_k1(1, k)   -q_k1(4, k)    q_k1(3, k);
                  q_k1(4, k)    q_k1(1, k)   -q_k1(2, k);
                 -q_k1(3, k)    q_k1(2, k)    q_k1(1, k)];
        Q = 0.25 * t^2 * cov_k * cov_g * cov_k';
        w = gyr(k, :);
        omega_w = [ 0  -w(1) -w(2) -w(3);
                  w(1)    0   w(3) -w(2);
                  w(2) -w(3)    0   w(1);
                  w(3)  w(2) -w(1)    0];
        A = eye(4) + 0.5 * omega_w * t;
        A_n(:, :, i) = A;
        i = i + 1;
        if (q_meas(:, k)' * q_k1(:, k)) < 0
            q_meas(:, k) = -q_meas(:, k);
        else
            q_meas(:, k) = q_meas(:, k);
        end
        % [qnb, att_meas(:, k), Cnb] = attsyn(q_meas(:, k));
        % [q_k(:, k+1), P_init] = kalman_filter(A, 0, eye(4), Q, R, q_meas(:, k), q_k(:, k), P_init);
        % q_k(:, k+1) = normalize_vector(q_k(:, k+1));
        % [qnb, att_k(:, k+1), Cnb] = attsyn(q_k(:, k+1));
    end
    q_k1(:, N0+K) = U_FIR_Filter(A_n, 0, eye(4), eye(4), N1, q_meas(:, N0 + K - N1 : N0 + K - 1), 0, N0+K);
    q_k1(:, N0+K) = normalize_vector(q_k1(:, N0+K));
    % % [qnb, att_k(:, N+K), Cnb] = attsyn(q_k(:, N+K));
    att_k1(:, N0+K) = q_att(q_k1(:, N0+K), 1);
    
    % 计算当前时刻与前一时刻的四元数内积
    %归一化角度变化
    delta_pitch = att_k1(1, N0+K) - att_k1(1, N0+K-1);
    delta_pitch = delta_pitch * 180 / pi;
    pitch_ratio = min(delta_pitch / theta1_max_deg, 1);


    % 基于theta动态调整N
    N1 = N_min + (N_max - N_min) * pitch_ratio;
    N1 = round(max(N_min, min(N_max, N1))); % 取整并限制在[min, max]
    A_n = zeros(4, 4, N1);
    N_k(N_min+K, 1) = N1;

end



for K = 1:n-N0
    i = 1;
    if K <= N_max
        N2 = N_min;
    end
    for k = N0 + K - N2 : N0 + K - 1
    % % %     if (K == 1) || (k == K+N-1)
    % % %         if (abs(norm(acc(k, :), 2) - 9.81)) < 0.5
    % % %             acc(k, :) = acc(k, :);
    % % %         else
    % % %             acc(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * g_G;
    % % %         end
    % % %     % % %     % if (abs(norm(mag(k, :), 2) - 0.49)) < 0.2
    % % %     % % %     %     mag(k, :) = mag(k, :);
    % % %     % % %     % else
    % % %     % % %     %     mag(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * h_G;
    % % %     % % %     % end
    % % %     end


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
        J = AQUA_J(acc(k, :), mag(k,:), l, q_acc, q_mag);
        R = J * cov_u * J';
        cov_k = [-q_k2(2, k)   -q_k2(3, k)   -q_k2(4, k);
                  q_k2(1, k)   -q_k2(4, k)    q_k2(3, k);
                  q_k2(4, k)    q_k2(1, k)   -q_k2(2, k);
                 -q_k2(3, k)    q_k2(2, k)    q_k2(1, k)];
        Q = 0.25 * t^2 * cov_k * cov_g * cov_k';
        w = gyr(k, :);
        omega_w = [ 0  -w(1) -w(2) -w(3);
                  w(1)    0   w(3) -w(2);
                  w(2) -w(3)    0   w(1);
                  w(3)  w(2) -w(1)    0];
        A = eye(4) + 0.5 * omega_w * t;
        A_n(:, :, i) = A;
        i = i + 1;
        if (q_meas(:, k)' * q_k2(:, k)) < 0
            q_meas(:, k) = -q_meas(:, k);
        else
            q_meas(:, k) = q_meas(:, k);
        end
        % [qnb, att_meas(:, k), Cnb] = attsyn(q_meas(:, k));
        % [q_k(:, k+1), P_init] = kalman_filter(A, 0, eye(4), Q, R, q_meas(:, k), q_k(:, k), P_init);
        % q_k(:, k+1) = normalize_vector(q_k(:, k+1));
        % [qnb, att_k(:, k+1), Cnb] = attsyn(q_k(:, k+1));
    end
    q_k2(:, N0+K) = U_FIR_Filter(A_n, 0, eye(4), eye(4), N2, q_meas(:, N0 + K - N2 : N0 + K - 1), 0, N0+K);
    q_k2(:, N0+K) = normalize_vector(q_k2(:, N0+K));
    % % [qnb, att_k(:, N+K), Cnb] = attsyn(q_k(:, N+K));
    att_k2(:, N0+K) = q_att(q_k2(:, N0+K), 1);
    
    % % % % 计算当前时刻与前一时刻的四元数内积
    % % % dot_q = abs(dot(q_k(:, N0 + K), q_k(:, N0 + K - 1))); 
    % % % dot_q = min(max(dot_q, -1), 1); % 限制dot_q范围，避免acos出错
    % % % theta_rad = 2 * acos(dot_q);  % 姿态变化角度（弧度）
    % % % theta_deg = rad2deg(theta_rad);  % 转换为度数
    % % % 
    % % % % 基于theta动态调整N
    % % % theta_ratio = min(theta_deg / theta2_max_deg, 1);  % 限制为[0, 1]
    % % % N2 = N_max - (N_max - N_min) * theta_ratio;
    % % % N2 = round(max(N_min, min(N_max, N2))); % 取整并限制在[min, max]
    % % % A_n = zeros(4, 4, N2);
    % % % N_k(N_min+K, 2) = N2;


    %归一化角度变化
    delta_yaw   = att_k2(3, N0+K) - att_k2(3, N0+K-1);
    delta_yaw = delta_yaw * 180 / pi;
    yaw_ratio   = min(delta_yaw   / theta2_max_deg, 1);


    % 基于theta动态调整N
    N2 = N_max - (N_max - N_min) * yaw_ratio;
    N2 = round(max(N_min, min(N_max, N2))); % 取整并限制在[min, max]
    A_n = zeros(4, 4, N2);
    N_k(N_min+K, 2) = N2;




end

    q_k(1:3, :) = q_k1(1:3, :);
    q_k(4, :) = q_k2(4, :);
    att_k(1:2, :) = att_k1(1:2, :);
    att_k(3, :)  =  att_k2(3, :);
    


