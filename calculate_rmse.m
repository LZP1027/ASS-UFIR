function rmse = calculate_rmse(estimated_values, true_values)
    

    % 计算估计值与真实值的差异向量
    errors = estimated_values - true_values;
    
    % 计算均方误差
    squared_errors = errors .^ 2;
    mse = mean(squared_errors, 2);
    
    % 计算均方根误差
    rmse = sqrt(mse);


end
