%% Function input parameters
load('F_tF_mid','F_tF');
IPT.X = X; IPT.Y = Y; 
IPT.ker = p;
IPT.var = noise_var;
IPT.F_tF = F_tF;
IPT.MIN_ITER = 10; IPT.MAX_ITER = 100; IPT.RELTOL = 1e-5;
IPT.cut = 1e-16;

%%
figure; imshow(X,[0,1]); title('Original Image'); % savefig('origin cameraman') 
figure; imshow(Y,[0,1]); title('Degraded Image'); % savefig('degraded cameraman') 

%%
IPT.alpha = 1; IPT.beta = 0.1;
[OPT_gam, history_gam] = solvefun_2d(IPT,0,'Gamma');
x_gamma = idct2(OPT_gam.C);
figure; imshow(x_gamma,[0,1]); 
title(['Gamma alpha = ', num2str(IPT.alpha),' beta = ', num2str(IPT.beta)]);
 
%%
[OPT_no, history_no] = solvefun_2d(IPT,0,'SBL');
x_no = idct2(OPT_no.C);
figure; imshow(x_no,[0,1]); title('SBL');
 
%%
IPT.beta = 0.02;
[OPT_gau, history_gau] = solvefun_2d(IPT,0,'Gauss');
x_gau = idct2(OPT_gau.C);
imshow(x_gau,[0,1]); title(['Gauss, beta = ', num2str(IPT.beta)]);

%%
err_gam = OPT_gam.relerr(end);
err_gau = OPT_gau.relerr(end);
err_no = OPT_no.relerr(end);

disp('Gamma hyperprior REL Err');
disp(err_gam);
disp('Gauss hyperprior REL Err');
disp(err_gau);
disp('No hyperprior REL Err');
disp(err_no);

num_gam = numel(find(abs(OPT_gam.C) < tol)) / (m*n);
num_gau = numel(find(abs(OPT_gau.C) < tol)) / (m*n);
num_no = numel(find(abs(OPT_no.C) < tol)) / (m*n);

disp('Gamma hyperprior sparsity');
disp(num_gam);
disp('Gauss hyperprior sparsity');
disp(num_gau);
disp('No hyperprior sparsity');
disp(num_no);

%%
cmin = -16.5; cmax = 5;

figure; imagesc(log(abs(C))); axis image; colormap jet; colorbar; caxis([cmin cmax]); %unify color scale
title('log|C_{X}| (True Coefficients)'); 

figure; imagesc(log(abs(dct2(Y)))); axis image; colormap jet; colorbar;  caxis([cmin cmax]);
title('log|C_{Y}| (Degraded Coefficients)'); 

figure; imagesc(log(abs(OPT_no.C))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{SBL}| (SBL Restored Coefficients)');

figure; imagesc(log(abs(OPT_gau.C))); axis image; colormap jet; colorbar; caxis([cmin cmax]);
title('log|C_{Gauss}| (Gauss Restored Coefficients)'); 
