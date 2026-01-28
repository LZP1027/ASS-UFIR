function x = Givens(A, b)
% Solve Ax=b using Givens rotation based QR decomposition

[m, n] = size(A);

% In-place modification for QR decomposition (upper triangular R)
for col = 1 : n
for row = m:-1:(col+1)
% Compute c and s for Givens rotation
[c, s] = ComputeGivensRotation(A(col, col), A(row, col));

% Apply rotation to A
for k = col:n
a = A(col, k);
b_elem = A(row, k);
A(col, k) = c * a - s * b_elem;
A(row, k) = s * a + c * b_elem;
end

% Apply rotation to b
a = b(col);
b_elem = b(row);
b(col) = c * a - s * b_elem;
b(row) = s * a + c * b_elem;
end
end

% Back substitution to solve R*x = d
x = zeros(n, 1);
for row = n:-1:1
if abs(A(row, row)) < 1e-12
error('Matrix is singular or ill-conditioned.');
end
tmp = A(row, row+1:end) * x(row+1:end);
x(row) = (b(row) - tmp) / A(row, row);
end

end

function [c, s] = ComputeGivensRotation(aii, aji)
% Compute cosine (c) and sine (s) for Givens rotation
if abs(aji) < 1e-12
c = 1;
s = 0;
elseif abs(aji) > abs(aii)
aii_over_aji = aii / aji;
one_over_sqrt = 1 / sqrt(1 + aii_over_aji^2);
c = - aii_over_aji * one_over_sqrt;
s = one_over_sqrt;
else
aji_over_aii = aji / aii;
one_over_sqrt = 1 / sqrt(1 + aji_over_aii^2);
c = one_over_sqrt;
s = - aji_over_aii * one_over_sqrt;
end
end