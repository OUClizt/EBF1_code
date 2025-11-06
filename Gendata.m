clc,clear,close all

%% % Generate data for script gam_gau, lp and gamma 
I = imread('cameraman_big.jpg'); I = rgb2gray(I);
X = im2double(I);
X = imresize(X,0.5,'nearest');
% I = imread('house.png'); X = im2double(I); % house image
[m,n] = size(X);
C = dct2(X);

%% % Model parameters
tol = 1e-8; % global tolerance
ker_size = 15; 
ker_std = 1; % size and standard deviation of Gaussian kernel, controlling blur level
sigma = 0.1; % noise percentage level

%% % Coefficient pruning
cc = zeros(m,n);
idex = (abs(C) > 0.025);
cc(idex) = C(idex);
X = idct2(cc);
C = cc;
num_origin = numel(find(abs(C) < tol)) / (m*n);
fprintf('Proportion of zero elements after truncation: %f\n', num_origin);

%% % Convolution
p = fspecial('gaussian', ker_size, ker_std);
K = imfilter(X, p, 'symmetric');

%% % Add noise
rng('default'); rng(1,'twister');
e = randn(size(X)); e = e/norm(e,'fro');
noise = sigma * e * norm(X,'fro');
Y = K + noise; 
C_Y = dct2(Y);
noise_var = (sigma^2 * norm(X,'fro')^2) / (m*n); % compute noise variance


