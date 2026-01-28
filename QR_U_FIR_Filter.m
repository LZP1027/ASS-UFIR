%%-------------state space equation----------
% x(k) = F * x(k-1) + E * u(k) + B * w(k)  
% y(k) = H * x(k) + v(k)

function x_UFIR = QR_U_FIR_Filter(F,E,B,H,N,N2,y,u,k)

dim_x = size(F,1);    

%% F_mk
F_mk = eye(4);
F_i = eye(4);
for i =1 : N-1  
    F_i = F_i * F(:, :, i);
    F_mk = [F_mk ; F_i];    
end

% %% S_mk
% S1_mk = F^0 * E;
% for i = 1 : N-1
%  S1_mk = [S1_mk ; F^i * E];
% end
% S_mk = S1_mk;
% for j = 2 : N      
%  S1_mk = zeros((j-1) * size(E,1) , size(E,2));
%    for i = 0 : N-j
%    S1_mk = [S1_mk ; F^i * E];
%    end  
%  S_mk = [S_mk  S1_mk];
% end
%S_mk = 0;

%% S_bar_mk
% S_bar_mk = S_mk(dim_x * (N-1) + 1:dim_x * N ,:);     % S_bar_mk
%S_bar_mk = 0;     % S_bar_mk

%% C_bar_mk
C_bar_mk = H;
for i = 1 : N-1
    C_bar_mk = blkdiag(C_bar_mk,H);
end

%去掉对应维数
dimto_del = [];
if N-N2 > 0
    for j = 0 : N-N2-1
        jdx = j *dim_x + (2:3);
        dimto_del = [dimto_del, jdx];
    end
    F_mk([dimto_del, jdx], :) = [];
    C_bar_mk([dimto_del, jdx], :) = [];
    C_bar_mk(:, [dimto_del, jdx]) = [];
end


%% H_mk  L_mk 
H_bar_mk = C_bar_mk * F_mk;
%%%%L_mk = C_bar_mk * S_mk;

%% H_mk  
H_mk = H_bar_mk *  (F_i)^(-1); 
%%%H_hat_mk = (H_mk' * H_mk)^(-1) * H_mk';


    
% % % m = k-N;
% % % Y = [y(:,m)];
% % % % U = [u(:,m)];
% % % U = 0;
% % % 
% % % for i = m+1 : k-1
% % %     Y = [Y ; y(:, i)];
% % %     % U = [U ; u(:,i)];
% % % end

Y = [y(:,1)];
% U = [u(:,m)];

for i = 2 : N
    Y = [Y ; y(:, i)];
    % U = [U ; u(:,i)];
end

%%%H_hat_mk = (H_mk' * H_mk)^(-1) * H_mk';

if N-N2 > 0
    Y(dimto_del, :) = []; 
end

% x_UFIR = H_hat_mk * Y + ( S_bar_mk - H_hat_mk * L_mk) * U;
%%%x_UFIR = H_hat_mk * Y;

%x_UFIR = (H_mk' * H_mk) \ (H_mk' * Y);

%%%%x_UFIR = (H_mk) \ (Y);

x_UFIR = Givens(H_mk, Y);
end