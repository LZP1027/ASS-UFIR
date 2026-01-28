function J = AQUA_J(a_L, m_L, l, q_acc, q_mag)
    gamma = l(1)^2 + l(2)^2;
    beta_1 = sqrt(gamma+l(1)*sqrt(gamma));
    beta_2 = sqrt(gamma-l(1)*sqrt(gamma));
if a_L(3) >= 0 
k = sqrt(1+a_L(3));
dq_f1 = [q_mag(1)     0       0       -q_mag(3)  q_acc(1)  -q_acc(2)  q_acc(3)     0    ;
            0     q_mag(1)  q_mag(4)      0      q_acc(2)   q_acc(1)     0      q_acc(3);
            0    -q_mag(4)  q_mag(1)      0      q_acc(3)       0     q_acc(1) -q_acc(2);
         q_mag(4)     0       0        q_mag(1)      0     -q_acc(3)  q_acc(2)  q_acc(1)]; 
  if l(1) > 0
df1_f2 = [0    0      1/k             0                      0             0;
          0  -2/k  a_L(2)/k^3         0                      0             0;
         2/k   0  -a_L(1)/k^3         0                      0             0;
          0    0       0              0                      0             0;
          0    0       0   l(2)^2/beta_1/gamma     l(1)*l(2)/beta_1/gamma  0;
          0    0       0              0                      0             0;
          0    0       0              0                      0             0;
          0    0       0  -l(2)*beta_1/gamma^1.5   l(1)*beta_1/gamma^1.5   0]/(2*sqrt(2));
  else 
df1_f2 = [0    0      1/k             0                      0             0;
          0  -2/k  a_L(2)/k^3         0                      0             0;
         2/k   0  -a_L(1)/k^3         0                      0             0;
          0    0       0              0                      0             0;
          0    0       0   l(2)*beta_2/gamma^1.5   l(1)*beta_2/gamma^1.5   0;
          0    0       0              0                      0             0;
          0    0       0              0                      0             0;
          0    0       0  -l(2)^2/beta_2/gamma     l(1)*l(2)/beta_2/gamma  0]/(2*sqrt(2));
  end
else 
k = sqrt(1-a_L(3));
dq_f1 = [q_mag(1)     0       0       -q_mag(3)  q_acc(1)  -q_acc(2)     0     -q_acc(4);
            0     q_mag(1)  q_mag(4)      0      q_acc(2)   q_acc(1) -q_acc(4)     0    ;
            0    -q_mag(4)  q_mag(1)      0          0      q_acc(4)  q_acc(1) -q_acc(2);
         q_mag(4)     0       0        q_mag(1)  q_acc(4)      0      q_acc(2)  q_acc(1)];
  if l(1) > 0
df1_f2 = [0  -2/k -a_L(2)/k^3         0                      0             0;
          0    0     -1/k             0                      0             0;
          0    0       0              0                      0             0;
         2/k   0   a_L(1)/k^3         0                      0             0;
          0    0       0   l(2)^2/beta_1/gamma    -l(1)*l(2)/beta_1/gamma  0;
          0    0       0              0                      0             0;
          0    0       0              0                      0             0;
          0    0       0  -l(2)*beta_1/gamma^1.5   l(1)*beta_1/gamma^1.5   0]/(2*sqrt(2));
  else
df1_f2 = [0  -2/k -a_L(2)/k^3         0                      0             0;
          0    0     -1/k             0                      0             0;
          0    0       0              0                      0             0;
         2/k   0   a_L(1)/k^3         0                      0             0;
          0    0       0   l(2)*beta_2/gamma^1.5  -l(1)*beta_2/gamma^1.5   0;
          0    0       0              0                      0             0;
          0    0       0              0                      0             0;
          0    0       0  -l(2)^2/beta_2/gamma     l(1)*l(2)/beta_2/gamma  0]/(2*sqrt(2));
     
  end
end

df2_u = [                  1                                       1                                               1                           0                    0             0  ;
                           1                                       1                                               1                           0                    0             0  ;
                           1                                       1                                               1                           0                    0             0  ;
         m_L(3)-(2*a_L(1)*m_L(1)+a_L(2)*m_L(2))/k^2      -a_L(1)*m_L(2)/k^2                  a_L(1)*(a_L(1)*m_L(1)+a_L(2)*m_L(2))/k^4    1-a_L(1)^2/k^2    -a_L(1)*a_L(2)/k^2  a_L(1);
                   -a_L(2)*m_L(1)/k^2           m_L(3)-(a_L(1)*m_L(1)+2*a_L(2)*m_L(2))/k^2   a_L(2)*(a_L(1)*m_L(1)+a_L(2)*m_L(2))/k^4  -a_L(1)*a_L(2)/k^2    1-a_L(2)^2/k^2    a_L(2);
                        -m_L(1)                                 -m_L(2)                                          m_L(3)                     -a_L(1)              -a_L(2)       a_L(3)];



J = dq_f1 * df1_f2 * df2_u;
end
