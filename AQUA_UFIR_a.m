%A Linear Kalman Filter for MARG Orientation Estimation 
% Using the Algebraic Quaternion Algorithm
function [q_k, att_k, q_meas, att_meas] = AQUA_UFIR_a(gyr, acc, mag, cov_g, cov_a, cov_m, n, N, t)
%定义L为为载体（传感器）坐标系，G为全局（地球）坐标系
%定义各测量向量和参考向量，均需归一化处理


cov_acc = [cov_a^2    0       0;
              0    cov_a^2    0;
              0       0    cov_a^2] / (9.81^2);
cov_mag = [cov_m^2    0       0;
              0    cov_m^2    0;
              0       0    cov_m^2] / (0.52^2);
cov_u = [cov_acc     zeros(3, 3);
         zeros(3, 3)    cov_mag];
w = zeros(1, 3);
% a_L = [a_x; a_y; a_z]; 
g_G = [0; 0; 9.81];
% m_L = [m_x; m_y; m_z];
% h_G = [0.5547; 0; 0.8321];
% h_G = [-0.6036; 0.2182; -1.0868];
% h_G = [0.6963; 0; 0.7177]
h_G = [0.35; 0; 0.35];
% w_L = [w_x; w_y; w_z];
P_init = [0.1   0     0     0;
            0   0.1   0     0;
            0     0   0.1   0;
            0     0     0   0.1];
q_meas = zeros(4, n);
q_k = zeros(4, n);
att_k = zeros(3, n);
att_meas = zeros(3, n);
q_k(:, 1) = [1; 0; 0; 0];

% N = 12;
A_n = zeros(4, 4, N);

for K = 1:n-N
    i = 1;
    for k = K:K+N-1
        % % % if (K == 1) || (k == K+N-1)
        % % %     if (abs(norm(acc(k, :), 2) - 9.81)) < 0.1
        % % %         acc(k, :) = acc(k, :);
        % % %     else
        % % %         acc(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * g_G;
        % % %     end
        % % %     if (abs(norm(mag(k, :), 2) - 1.2621)) < 0.2
        % % %         mag(k, :) = mag(k, :);
        % % %     else
        % % %         mag(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * h_G;
        % % %     end
        % % % end
        % % % % % % if (K == 1) || (k == K+N-1)
        % % % % % %     if (abs(norm(acc(k, :), 2) - 9.81)) < 0.5
        % % % % % %         acc(k, :) = acc(k, :);
        % % % % % %     else
        % % % % % %         acc(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * g_G;
        % % % % % %     end
        % % % % % %     % if (abs(norm(mag(k, :), 2) - 0.49)) < 0.2
        % % % % % %     %     mag(k, :) = mag(k, :);
        % % % % % %     % else
        % % % % % %     %     mag(k, :) =  quaternion_to_rotation_matrix(q_k(:, k)) * h_G;
        % % % % % %     % end
        % % % % % % end


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
        %%%%%A_n(:, :, k) = A;
        A_n(:, :, i) = A;
        i = i + 1;
        if (q_meas(:, k)' * q_k(:, k)) < 0
            q_meas(:, k) = -q_meas(:, k);
        else
            q_meas(:, k) = q_meas(:, k);
        end
        % [qnb, att_meas(:, k), Cnb] = attsyn(q_meas(:, k));
        % [q_k(:, k+1), P_init] = kalman_filter(A, 0, eye(4), Q, R, q_meas(:, k), q_k(:, k), P_init);
        % q_k(:, k+1) = normalize_vector(q_k(:, k+1));
        % [qnb, att_k(:, k+1), Cnb] = attsyn(q_k(:, k+1));
    end
    q_k(:, N+K) = U_FIR_Filter(A_n, 0, eye(4), eye(4), N, q_meas(:, K:K+N-1), 0, N+K);
    q_k(:, N+K) = normalize_vector(q_k(:, N+K));
    % % [qnb, att_k(:, N+K), Cnb] = attsyn(q_k(:, N+K));
    att_k(:, N+K) = q_att(q_k(:, N+K), 1);
end
    
end


