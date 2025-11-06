%% Function input parameters
load('F_tF_mid','F_tF');
IPT.X = X;IPT.Y = Y; 
IPT.ker = p;
IPT.var = noise_var;
IPT.F_tF = F_tF;
IPT.beta = 0.1;
IPT.MIN_ITER = 10; IPT.MAX_ITER = 100; IPT.RELTOL=1e-5;
IPT.cut = 1e-16;


%%
IPT.p = 0.5;
[OPT_l5, history_l5] = solvefun_2d(IPT,0,'lp');
x_l5 = idct2(OPT_l5.C);
imshow(x_l5,[0,1]);title(['General Gauss p = ', num2str(IPT.p),' beta = ', num2str(IPT.beta)])

%%
IPT.p = 0.75;
[OPT_l75, history_l75] = solvefun_2d(IPT,0,'lp');
x_l75 = idct2(OPT_l75.C);
imshow(x_l75,[0,1]);title(['General Gauss p = ', num2str(IPT.p),' beta = ', num2str(IPT.beta)])
 
%%
IPT.alpha = 1; IPT.beta = 0.1;
[OPT_l1, history_gam]=solvefun_2d(IPT,0,'Gamma');
x_gamma = idct2(OPT_l1.C);
figure;imshow(x_gamma,[0,1]);title(['General Gauss p = ', num2str(IPT.p),' beta = ', num2str(IPT.beta)]);

%%
err_l5 = OPT_l5.relerr(end);
err_l75 = OPT_l75.relerr(end);
err_l1 = OPT_l1.relerr(end);

disp('l5 hpyerprior REL Err ');
disp(err_l5);
disp('l75 hpyerprior REL Err ');
disp(err_l75);
disp('l1 hpyerprior REL Err ');
disp(err_l1);

num_l5=numel(find(abs(OPT_l5.C)<tol))/(m*n);
num_l75=numel(find(abs(OPT_l75.C)<tol))/(m*n);
num_l1=numel(find(abs(OPT_l1.C)<tol))/(m*n);
%%
cmin = -16.5;cmax = 5;

figure;imagesc(log(abs(C))); axis image; colormap jet;colorbar;caxis([cmin cmax]); % Í³Ò»ÑÕÉ«
title('log|C_{X}| (True Coefficients)');

figure; imagesc(log(abs(dct2(Y)))); axis image; colormap jet;colorbar;%caxis([cmin cmax]);
title('log|C_{Y}| (Degraded Coefficients)');
figure; imagesc(log(abs(OPT_l5.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{l0.5}| (l0.5 Restored Coefficients)');
figure; imagesc(log(abs(OPT_l75.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{l0.75}| (l0.75 Restored Coefficients)');
figure; imagesc(log(abs(OPT_l1.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);
title('log|C_{l1}| (l1 Restored Coefficients)');

