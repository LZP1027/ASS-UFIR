function C = quaternion_to_rotation_matrix(q)
    % 根据四元数的四个分量定义四元数的虚部和实部
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % 计算旋转矩阵的元素
    C = [
        q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2);
        2*(q1*q2 + q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
        2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0^2 - q1^2 - q2^2 + q3^2;
    ];
end
