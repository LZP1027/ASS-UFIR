clearvars;
n = 600;
N_min = 5;
N_max = 50;
N_list = N_min:N_max;
time = linspace(0, n/100, n);
t_k_1 = time(1:n);
t_k_2 = time(2:n);
cov_g = 0.0012; %0.01量级下误差几乎一样
cov_a = 0.01;
cov_m = 0.01;

rmse_set = zeros(3, N_max);

[q_true, att_true, gyr, acc, mag] = generate_simulink_data(n, 0, 0, 0); 
att_true = att_true * 180 / pi;
% % % [GCF_q_k, q_true, q_yk, att, att_y, att_error] = GCF(gyr, acc, mag, q_true, att_true, n);


%%%load("test0.mat", 'gyr', 'acc', 'mag', 'att_true');     %%%%%目前最好的图
load("test1-0.mat", 'gyr', 'acc', 'mag', 'att_true');

 

%%AQUA-UFIR_batch
for N = N_min:N_max
    fprintf('Running AQUA_UFIR_a: N = %d\n', N);
    tic;
    [q_k_UFIR, att_k_UFIR, q_meas, att_meas] = AQUA_UFIR_a(gyr, acc, mag, cov_g, cov_a, cov_m, n, N);
    time_ufir_batch = toc;                                   
    q_rmse_UFIR = calculate_rmse(q_k_UFIR(:, N:n), q_true(:, N:n));
    att_k_UFIR = att_k_UFIR * 180 / pi;
    att_rmse_UFIR = calculate_rmse(att_k_UFIR(:, N:n-1), att_true(:, N-1:n-2));
    rmse_set(:,N) = att_rmse_UFIR;
end




% ---------- RMSE 三子图（Roll / Pitch / Yaw），不统一 y 轴 ----------
figure;

axis_names = {'Roll','Pitch','Yaw'};
markers    = {'-o','-s','-^'};

for ax = 1:3
    subplot(3,1,ax);
    % 注意：用 rmse_set(ax, :) 而不是 rmse_set(ax, N_list)
    plot(N_list, rmse_set(ax, N_list), markers{ax}, 'LineWidth', 1.6, 'MarkerSize', 3); 
    hold on; grid on;
    xlabel('Horizon N'); ylabel([axis_names{ax} ' RMSE (deg)']);
    %title(sprintf('%s RMSE vs N (AQUA\\_UFIR\\_a)', axis_names{ax}));

    % 标出该轴最优 N（最小 RMSE）
    [min_val, min_idx] = min(rmse_set(ax, N_list));
    bestN_axis = N_list(min_idx);
    plot(bestN_axis, min_val, 'rp', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    xline(bestN_axis, '--r', sprintf('N=%d', bestN_axis), ...
          'LabelVerticalAlignment','bottom', 'LineWidth', 1);

    % 文本标注（点上方少量偏移）
    local_min = min(rmse_set(ax, :));
    local_max = max(rmse_set(ax, :));
    pad_local = 0.08 * (local_max - local_min + eps);
    y_text = min_val + 0.04 * (local_max - local_min + eps);
    text(bestN_axis, y_text, sprintf('min RMSE = %.4f°', min_val), ...
        'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');

    % 不统一 y 轴范围：每个子图按自身数据自适应并加少量 padding
    ylim([max(0, local_min - pad_local), local_max + pad_local]);

    hold off;
end

sgtitle('Attitude RMSE vs Horizon N (UFIR-batch)');







 

