function normalized_vector = normalize_vector(input_vector)
    % 计算向量的模
    input_vector = double(input_vector);
    vector_norm = norm(input_vector);
    
    %确保向量不是零向量
    % if vector_norm == 0
    %     error('input is a zero vector');
    % end
    % 对向量进行归一化处理
    normalized_vector = input_vector / vector_norm;
end
