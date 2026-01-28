function noise = gauss(sigma2,u)
%% info:
% 产生第i个元素是均值为ui，方差为sigma2i的i*1的高斯噪声列向量
% 输入参数：
% sigma2：协方差矩阵，对角阵
% u：均值，列向量，默认0
% 输出参数：
% noise：与维数参数匹配高斯噪声列向量
%% init
if nargin < 2
    u = 0;
end
n = size(u,1);
if n ~= size(sigma2,1)
    disp('均值与方差维数不匹配');
    return;
end
%% yield noise
noise = u + sqrt(diag(sigma2)) .* randn(n,1);
