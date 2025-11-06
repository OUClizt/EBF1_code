% The results of EBF with the half-Laplace hyperprior for restoring images degraded by different noise levels.
clc, clear, close all

%% % Generate data
I = imread('cameraman_big.jpg'); I = rgb2gray(I);
X = im2double(I);
X = imresize(X, 0.5, 'nearest');
% I = imread('house.png'); X = im2double(I); % house image
[m, n] = size(X);
C = dct2(X);

%% % Model parameters
tol = 1e-8; % global tolerance
ker_size = 15; ker_std = 0.5; % size and standard deviation of Gaussian kernel, controlling blur level

%% % Coefficient pruning
cc = zeros(m, n);
idex = (abs(C) > 0.025);
cc(idex) = C(idex);
X = idct2(cc);
C = cc;
num_origin = numel(find(abs(C) < tol)) / (m * n);
fprintf('Proportion of zero elements after truncation: %f\n', num_origin);

%% % Convolution
p = fspecial('gaussian', ker_size, ker_std);
K = imfilter(X, p, 'symmetric');

%% % Function input parameters
load('F_tF_small', 'F_tF');
IPT.F_tF = F_tF;
IPT.X = X;
IPT.ker = p;
IPT.alpha = 1; IPT.beta = 0.1;
IPT.MIN_ITER = 10; IPT.MAX_ITER = 200; IPT.RELTOL = 1e-5;
IPT.cut = 1e-16;

%% % Add noise
rng('default'); rng(1, 'twister');
e = randn(size(X)); e = e / norm(e, 'fro');
sigma_1 = 0.05;
sigma_2 = 0.1;
sigma_3 = 0.2;
noise_1 = sigma_1 * e * norm(X, 'fro');
noise_2 = sigma_2 * e * norm(X, 'fro');
noise_3 = sigma_3 * e * norm(X, 'fro');
Y_1 = K + noise_1;
Y_2 = K + noise_2;
Y_3 = K + noise_3;

%%
figure; imshow(X, [0, 1]); title('Original Image');
figure; imshow(Y_1, [0, 1]); title('Degraded Image (5% noise level)');
figure; imshow(Y_2, [0, 1]); title('Degraded Image (10% noise level)');
figure; imshow(Y_3, [0, 1]); title('Degraded Image (20% noise level)');

%% % Restoration under different noise levels
noise_var1 = (sigma_1^2 * norm(X, 'fro')^2) / (m * n); % compute noise variance
IPT.Y = Y_1; IPT.var = noise_var1;
[OPT_gam1, history_gam1] = solvefun_2d(IPT, 0, 'Gamma');
x_gamma1 = idct2(OPT_gam1.C);
figure; imshow(x_gamma1, [0, 1]); title(['Noise level = ', num2str(sigma_1)]);

noise_var2 = (sigma_2^2 * norm(X, 'fro')^2) / (m * n); % compute noise variance
IPT.Y = Y_2; IPT.var = noise_var2;
[OPT_gam2, history_gam2] = solvefun_2d(IPT, 0, 'Gamma');
x_gamma2 = idct2(OPT_gam2.C);
figure; imshow(x_gamma2, [0, 1]); title(['Noise level = ', num2str(sigma_2)]);

noise_var3 = (sigma_3^2 * norm(X, 'fro')^2) / (m * n); % compute noise variance
IPT.Y = Y_3; IPT.var = noise_var3;
[OPT_gam3, history_gam3] = solvefun_2d(IPT, 0, 'Gamma');
x_gamma3 = idct2(OPT_gam3.C);
figure; imshow(x_gamma3, [0, 1]); title(['Noise level = ', num2str(sigma_3)]);

%% % Report relative errors and sparsity
err_gam1 = OPT_gam1.relerr(end);
err_gam2 = OPT_gam2.relerr(end);
err_gam3 = OPT_gam3.relerr(end);

disp('5% noise level REL Err');
disp(err_gam1);
disp('10% noise level REL Err');
disp(err_gam2);
disp('20% noise level REL Err');
disp(err_gam3);

num_gam1 = numel(find(abs(OPT_gam1.C) < tol)) / (m * n);
num_gam2 = numel(find(abs(OPT_gam2.C) < tol)) / (m * n);
num_gam3 = numel(find(abs(OPT_gam3.C) < tol)) / (m * n);

disp('5% noise level sparsity');
disp(num_gam1);
disp('10% noise level sparsity');
disp(num_gam2);
disp('20% noise level sparsity');
disp(num_gam3);

%% % Visualization of coefficient maps
cmin = -16.5; cmax = 5;

figure; imagesc(log(abs(C))); axis image; colormap jet; colorbar; caxis([cmin cmax]); % unify color scale
title('log|C_{X}| (True Coefficients)'); axis image off;

figure; imagesc(log(abs(dct2(Y_1)))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Y}| (Degraded Coefficients, 5% noise)'); axis image off;
figure; imagesc(log(abs(dct2(Y_2)))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Y}| (Degraded Coefficients, 10% noise)'); axis image off;
figure; imagesc(log(abs(dct2(Y_3)))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Y}| (Degraded Coefficients, 20% noise)'); axis image off;

figure; imagesc(log(abs(OPT_gam1.C))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Gamma1}| (Restored Coefficients, 5% noise)'); axis image off;
figure; imagesc(log(abs(OPT_gam2.C))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Gamma2}| (Restored Coefficients, 10% noise)'); axis image off;
figure; imagesc(log(abs(OPT_gam3.C))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Gamma3}| (Restored Coefficients, 20% noise)'); axis image off;

