% The results of EBF with the half-Laplace hyperprior for restoring images degraded by different ill-posedness.
clc,clear,close all

%% % Generate data
I = imread('cameraman_big.jpg'); I = rgb2gray(I);
X = im2double(I);
X = imresize(X,0.5,'nearest');
% I = imread('house.png'); X = im2double(I); % house image
[m,n] = size(X);
C = dct2(X);

%% % Model parameters
tol = 1e-8; % global tolerance
ker_size = 15; % size of Gaussian kernel
sigma = 0.1; % noise level

%% % Coefficient pruning
cc = zeros(m,n);
idex = (abs(C) > 0.025);
cc(idex) = C(idex);
X = idct2(cc);
C = cc;
num_origin = numel(find(abs(C) < tol)) / (m*n);
fprintf('Proportion of zero elements after truncation: %f\n', num_origin);

%% % Add noise
rng('default'); rng(1,'twister');
e = randn(size(X)); e = e/norm(e,'fro');
noise = sigma * e * norm(X,'fro');
noise_var = (sigma^2 * norm(X,'fro')^2) / (m*n); % compute noise variance

%% % Input parameters for the function
IPT.X = X;
IPT.var = noise_var;
IPT.alpha = 1; IPT.beta = 0.1;
IPT.MIN_ITER = 10; IPT.MAX_ITER = 200; IPT.RELTOL=1e-5;
IPT.cut = 1e-16;

%% % Convolution
ker_std_1 = 0.5; p_1 = fspecial('gaussian', ker_size, ker_std_1);
ker_std_2 = 1.0; p_2 = fspecial('gaussian', ker_size, ker_std_2);
ker_std_3 = 1.5; p_3 = fspecial('gaussian', ker_size, ker_std_3);

K_1 = imfilter(X, p_1, 'symmetric');
K_2 = imfilter(X, p_2, 'symmetric');
K_3 = imfilter(X, p_3, 'symmetric');

Y_1 = K_1 + noise; 
Y_2 = K_2 + noise; 
Y_3 = K_3 + noise; 

%% % Display degraded images
figure;imshow(X,[0,1]);title('Original image');
figure;imshow(Y_1,[0,1]);title('\sigma_{ker}=0.5 degraded image');
figure;imshow(Y_2,[0,1]);title('\sigma_{ker}=1 degraded image');
figure;imshow(Y_3,[0,1]);title('\sigma_{ker}=1.5 degraded image');

%% % Restoration under different kernel standard deviations
load('F_tF_small','F_tF');IPT.F_tF = F_tF;
IPT.Y = Y_1; IPT.ker = p_1;
[OPT_gam1, history_gam1]=solvefun_2d(IPT,0,'Gamma');
x_gamma1 = idct2(OPT_gam1.C);
figure;imshow(x_gamma1,[0,1]);title(['\sigma_{ker} = ', num2str(ker_std_1)]);

load('F_tF_mid','F_tF');IPT.F_tF = F_tF;
IPT.Y = Y_2; IPT.ker = p_2;
[OPT_gam2, history_gam2]=solvefun_2d(IPT,0,'Gamma');
x_gamma2 = idct2(OPT_gam2.C);
figure;imshow(x_gamma2,[0,1]);title(['\sigma_{ker} = ', num2str(ker_std_2)]);

load('F_tF_big','F_tF');IPT.F_tF = F_tF;
IPT.Y = Y_3; IPT.ker = p_3;
[OPT_gam3, history_gam3]=solvefun_2d(IPT,0,'Gamma');
x_gamma3 = idct2(OPT_gam3.C);
figure;imshow(x_gamma3,[0,1]);title(['\sigma_{ker} = ', num2str(ker_std_3)]);

%% % Compute relative errors and sparsity
err_gam1 = OPT_gam1.relerr(end);
err_gam2 = OPT_gam2.relerr(end);
err_gam3 = OPT_gam3.relerr(end);

disp('ker_std=0.5 REL Err ');
disp(err_gam1);
disp('ker_std=1 REL Err ');
disp(err_gam2);
disp('ker_std=1.5 REL Err ');
disp(err_gam3);

num_gam1=numel(find(abs(OPT_gam1.C)<tol))/(m*n);
num_gam2=numel(find(abs(OPT_gam2.C)<tol))/(m*n);
num_gam3=numel(find(abs(OPT_gam3.C)<tol))/(m*n);

disp('ker_std=0.5 sparsity ');
disp(num_gam1);
disp('ker_std=1 sparsity ');
disp(num_gam2);
disp('ker_std=1.5 sparsity ');
disp(num_gam3);

%% % Display DCT coefficient maps
cmin = -16.5; cmax =  5;

figure;imagesc(log(abs(C))); axis image; colormap jet;colorbar;caxis([cmin cmax]); 
title('log|C_{X}| (True coefficients)');axis image off; 

figure; imagesc(log(abs(dct2(Y_1)))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{Y}| (\sigma_{ker}=0.5 degraded coefficients)');axis image off;
figure; imagesc(log(abs(dct2(Y_2)))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{Y}| (\sigma_{ker}=1 degraded coefficients)');axis image off;
figure; imagesc(log(abs(dct2(Y_3)))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{Y}| (\sigma_{ker}=1.5 degraded coefficients)');axis image off;

figure; imagesc(log(abs(OPT_gam1.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{gamma1}| (\sigma_{ker}=0.5 restored coefficients)');axis image off;
figure; imagesc(log(abs(OPT_gam2.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{gamma2}| (\sigma_{ker}=1 restored coefficients)');axis image off;
figure; imagesc(log(abs(OPT_gam3.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{gamma3}| (\sigma_{ker}=1.5 restored coefficients)');axis image off;
