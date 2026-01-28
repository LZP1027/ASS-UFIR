function out = rot_angle_error_metrics(q_est, q_true, angle_unit)
%ROT_ANGLE_ERROR_METRICS  Global attitude error based on rotation angle.
%
% Inputs
%   q_est:  4xN estimated quaternions, format [q0;q1;q2;q3], q0 is scalar part
%   q_true: 4xN reference quaternions, same format
%   angle_unit: 'deg' or 'rad' (optional, default 'deg')
%
% Output (struct)
%   out.theta_q     : 1xK rotation-angle error from quaternion definition
%   out.theta_R     : 1xK rotation-angle error from rotation-matrix definition
%   out.rmse_q      : scalar RMSE of out.theta_q
%   out.rmse_R      : scalar RMSE of out.theta_R
%   out.meanabs_q   : scalar mean absolute of out.theta_q
%   out.meanabs_R   : scalar mean absolute of out.theta_R
%   out.valid_idx   : indices used (after removing NaN/Inf)
%
% Notes
%   - Quaternion sign ambiguity is handled in quaternion-angle metric via abs(qe0).
%   - Both inputs are normalized internally for robustness.
%   - If input contains NaN/Inf columns, they are ignored.

    if nargin < 3 || isempty(angle_unit)
        angle_unit = 'deg';
    end

    % Basic checks
    if size(q_est,1) ~= 4 || size(q_true,1) ~= 4
        error('q_est and q_true must be 4xN with proper quaternion format.');
    end

    N = min(size(q_est,2), size(q_true,2));
    q_est  = q_est(:,1:N);
    q_true = q_true(:,1:N);

    % Remove invalid columns
    valid = all(isfinite(q_est),1) & all(isfinite(q_true),1);
    q_est  = q_est(:,valid);
    q_true = q_true(:,valid);
    idx_valid = find(valid);

    if isempty(idx_valid)
        error('No valid samples after removing NaN/Inf columns.');
    end

    % Normalize quaternions (column-wise)
    q_est  = normalize_quat(q_est);
    q_true = normalize_quat(q_true);

    % ---------- 1) Quaternion-defined rotation angle ----------
    % qe = inv(q_true) ⊗ q_est
    q_true_inv = quat_inv(q_true);
    q_e = quat_mul(q_true_inv, q_est);

    % theta = 2*acos(|qe0|)
    c = abs(q_e(1,:));
    c = clamp(c, 0, 1);                 % due to numerical errors
    theta_q = 2 * acos(c);              % rad

    % ---------- 2) Rotation-matrix-defined rotation angle ----------
    K = size(q_est,2);
    theta_R = zeros(1,K);
    for k = 1:K
        R_true = quat_to_rotm(q_true(:,k));
        R_est  = quat_to_rotm(q_est(:,k));
        R_e = R_true.' * R_est;
        cos_th = (trace(R_e) - 1) / 2;
        cos_th = clamp(cos_th, -1, 1);
        theta_R(k) = acos(cos_th);      % rad
    end

    % Unit convert
    if strcmpi(angle_unit, 'deg')
        theta_q = theta_q * (180/pi);
        theta_R = theta_R * (180/pi);
    elseif ~strcmpi(angle_unit, 'rad')
        error('angle_unit must be deg or rad.');
    end

    % Global metrics (use RMSE as a global indicator)
    out.theta_q = theta_q;
    out.theta_R = theta_R;

    out.rmse_q = sqrt(mean(theta_q.^2));
    out.rmse_R = sqrt(mean(theta_R.^2));

    out.meanabs_q = mean(abs(theta_q));
    out.meanabs_R = mean(abs(theta_R));

    out.valid_idx = idx_valid;
end

% ========================= Helper functions =========================

function qn = normalize_quat(q)
    n = sqrt(sum(q.^2,1));
    n = max(n, eps);
    qn = q ./ n;
end

function qi = quat_inv(q)
% Inverse of quaternion (column-wise). Works for non-unit too.
% q = [w;x;y;z]
    w = q(1,:);
    v = q(2:4,:);
    n2 = w.^2 + sum(v.^2,1);
    n2 = max(n2, eps);
    qi = [w; -v] ./ n2;
end

function qc = quat_conj(q)
    qc = [q(1,:); -q(2:4,:)];
end

function q = quat_mul(q1, q2)
% Quaternion product, scalar-first, column-wise
% q = q1 ⊗ q2
% q1,q2 are 4xN or 4x1 and broadcastable by same N after indexing outside.
    w1 = q1(1,:); x1 = q1(2,:); y1 = q1(3,:); z1 = q1(4,:);
    w2 = q2(1,:); x2 = q2(2,:); y2 = q2(3,:); z2 = q2(4,:);

    w = w1.*w2 - x1.*x2 - y1.*y2 - z1.*z2;
    x = w1.*x2 + x1.*w2 + y1.*z2 - z1.*y2;
    y = w1.*y2 - x1.*z2 + y1.*w2 + z1.*x2;
    z = w1.*z2 + x1.*y2 - y1.*x2 + z1.*w2;

    q = [w; x; y; z];
end

function R = quat_to_rotm(q)
% Convert scalar-first unit quaternion to rotation matrix
% q = [w;x;y;z]
    w = q(1); x = q(2); y = q(3); z = q(4);

    % Ensure unit length (robust)
    n = sqrt(w*w + x*x + y*y + z*z);
    if n < eps
        R = eye(3);
        return;
    end
    w = w/n; x = x/n; y = y/n; z = z/n;

    R = [ ...
        1 - 2*(y*y + z*z),     2*(x*y - w*z),       2*(x*z + w*y); ...
        2*(x*y + w*z),         1 - 2*(x*x + z*z),   2*(y*z - w*x); ...
        2*(x*z - w*y),         2*(y*z + w*x),       1 - 2*(x*x + y*y) ...
    ];
end

function y = clamp(x, a, b)
    y = min(max(x, a), b);
end
