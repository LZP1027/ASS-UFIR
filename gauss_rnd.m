%GAUSS_RND  Multivariate Gaussian random variables
% Syntax:
%   P = GAUSS_RND(M,S,N)
%
% In:
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%   N - Number of samples (optional, default 1)
%
% Out: 
%   X - DxN matrix of samples.
%   
% Description:
%   Draw N samples from multivariate Gaussian distribution
% 
%     X ~ N(M,S)

function X = gauss_rnd(M,S,N)

  if nargin < 3
    N = 1;
  end
  
  L = chol(S)';
  X = repmat(M,1,N) + L*randn(size(M,1),N);
  
  