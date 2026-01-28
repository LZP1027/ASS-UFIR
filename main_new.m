clearvars;
n = 1000;
N = 29;
N2 = 5;
N_max = N;
N0 = 5;
time = 0.01;
t_k = linspace(0, n*time, n);
% % t_UFIR = t_k(N-1:n-1);
% % t_UFIR_N = t_k(N0-1:n-1);
% % t_UFIR_q = t_k(N:n-1);
% % t_UFIR_N_q = t_k(N0:n-1);
t_UFIR = t_k(N:n);
t_UFIR_N = t_k(N0:n);
t_UFIR_q = t_k(N+1:n);
t_UFIR_N_q = t_k(N0+1:n);
cov_g = 0.0012; 
cov_a = 0.01;
cov_m = 0.01;
% cov_g = 0.00175; 
% cov_a = 0.0078;
% cov_m = 0.0050;

[q_true, att_true, gyr, acc, mag] = generate_simulink_data(n, 0, 0, 0);
att_true = att_true * 180 / pi;








%%N-5
tic;
[q_k_UFIR_5, att_k_UFIR_5, q_meas_2, att_meas] = AQUA_UFIR_a(gyr, acc, mag, cov_g, cov_a, cov_m, n, N2, time);
time_norm_5 = toc;
out_5 = rot_angle_error_metrics(q_k_UFIR_5(:, N2+1:n), q_true(:, N2+1:n));
q_rmse_UFIR_5 = calculate_rmse(q_k_UFIR_5(:, N2+1:n), q_true(:, N2+1:n));
att_k_UFIR_5 = att_k_UFIR_5 * 180 / pi;
att_rmse_UFIR_5 = calculate_rmse(att_k_UFIR_5(:, N2+1:n), att_true(:, N2+1:n));


%%N-29
tic;
[q_k_UFIR, att_k_UFIR, q_meas_2, att_meas] = AQUA_UFIR_a(gyr, acc, mag, cov_g, cov_a, cov_m, n, N, time);
time_norm = toc;
out_norm = rot_angle_error_metrics(q_k_UFIR(:, N+1:n), q_true(:, N+1:n));
q_rmse_UFIR = calculate_rmse(q_k_UFIR(:, N+1:n), q_true(:, N+1:n));
att_k_UFIR = att_k_UFIR * 180 / pi;
att_rmse_UFIR = calculate_rmse(att_k_UFIR(:, N+1:n), att_true(:, N+1:n));


% % % % %%AQUA-UFIR_N清除测量值
tic;
[q_k_del, att_k_del, q_meas_3, att_meas, N_k_3] = UFIR_Delete(gyr, acc, mag, cov_g, cov_a, cov_m, n, time, N_max);
time_del = toc;
out_del = rot_angle_error_metrics(q_k_del(:, N0+1:n), q_true(:, N0+1:n));
q_rmse_del = calculate_rmse(q_k_del(:, N0+1:n), q_true(:, N0+1:n));
att_k_del = att_k_del * 180 / pi;
att_k_del(3, N0+1:n) = yaw_calibration(att_k_del(3, N0+1:n), att_true(3, N0+1:n));
att_rmse_del = calculate_rmse(att_k_del(:, N0+1:n), att_true(:, N0+1:n));

% % % % %%AQUA-UFIR_N清除测量值
tic;
[q_k_QR, att_k_QR, q_meas_4, att_meas, N_k_4] = QR_UFIR_Delete(gyr, acc, mag, cov_g, cov_a, cov_m, n, time, 18);
time_QR = toc;
out_QR = rot_angle_error_metrics(q_k_QR(:, N0+1:n), q_true(:, N0+1:n));
q_rmse_QR = calculate_rmse(q_k_QR(:, N0+1:n), q_true(:, N0+1:n));
att_k_QR = att_k_QR * 180 / pi;
att_k_QR(3, N0+1:n) = yaw_calibration(att_k_QR(3, N0+1:n), att_true(3, N0+1:n));
att_rmse_QR = calculate_rmse(att_k_QR(:, N0+1:n), att_true(:, N0+1:n));




% 绘制四元数比较曲线
figure;

% 实部比较曲线
subplot(1, 4, 1);
plot(t_k, q_true(1, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR_q, q_k_UFIR(1, N+1:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N_q, q_k_UFIR_5(1, N2+1:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N_q, q_k_QR(1, N0+1:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('q0');
grid on;

% 虚部 i 比较曲线
subplot(1, 4, 2);
plot(t_k, q_true(2, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, q_k_UFIR(2, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_UFIR_5(2, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_QR(2, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('q1');
legend('Reference', 'UFIR-29', 'UFIR-5', 'SS-UFIR', 'Location', 'Best');
grid on;

% 虚部 j 比较曲线
subplot(1, 4, 3);
plot(t_k, q_true(3, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, q_k_UFIR(3, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_UFIR_5(3, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_QR(3, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('q2');
grid on;

% 虚部 k 比较曲线
subplot(1, 4, 4);
plot(t_k, q_true(4, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, q_k_UFIR(4, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_UFIR_5(4, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, q_k_QR(4, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('q3');
grid on;





% 绘制姿态角估计曲线与真值曲线
figure;

% 绘制俯仰角曲线
subplot(3, 1, 1);
plot(t_k, att_true(1, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, att_k_UFIR(1, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_UFIR_5(1, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_del(1, N0:n), '--', 'Color', [0, 0, 0], 'LineWidth', 2);
plot(t_UFIR_N, att_k_QR(1, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('Roll(deg)');
grid on;


% 绘制横滚角曲线
subplot(3, 1, 2);
plot(t_k, att_true(2, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, att_k_UFIR(2, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_UFIR_5(2, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_del(2, N0:n), '--', 'Color', [0, 0, 0], 'LineWidth', 2);
plot(t_UFIR_N, att_k_QR(2, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('Pitch(deg)');
grid on;
legend('Reference', 'UFIR-29', 'UFIR-5', 'SS-UFIR', 'ASS-UFIR', 'Location', 'Best');

% 绘制航向角曲线
subplot(3, 1, 3);
plot(t_k, att_true(3, 1:n), '-.', 'Color', [0.07,0.62,1.00], 'LineWidth', 3);
hold on; % 保持当前图形，以便绘制其他曲线
plot(t_UFIR, att_k_UFIR(3, N:n), '-', 'Color', [1, 0.5, 0], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_UFIR_5(3, N2:n), '--', 'Color', [1, 0.85, 0.25], 'LineWidth', 1.5);
plot(t_UFIR_N, att_k_del(3, N0:n), '--', 'Color', [0, 0, 0], 'LineWidth', 2);
plot(t_UFIR_N, att_k_QR(3, N0:n), '-', 'Color', [0.70,0.45,0.80], 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('Yaw(deg)');
grid on;






