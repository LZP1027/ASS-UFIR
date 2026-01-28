function att = q_att(q,i)

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    att = zeros(3, 1);
    if i == 1 %按XYZ顺序旋转
       att(1) = atan2(2*(q0*q1 + q2*q3), 1-2*(q1^2 + q2^2));
       att(2) = asin(2*(q0*q2 - q3*q1));
       att(3) = atan2(2*(q0*q3 + q2*q1), 1-2*(q3^2 + q2^2));
    else if i == 2
       att(1) = atan2(2*(q0*q3 + q2*q1), 1-2*(q3^2 + q2^2));;
       att(2) = asin(2*(q0*q2 - q3*q1));
       att(3) = atan2(2*(q0*q1 + q2*q3), 1-2*(q1^2 + q2^2));    
    end

end