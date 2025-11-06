%% Function input parameters
load('F_tF_mid','F_tF');
IPT.X = X;IPT.Y = Y; 
IPT.ker = p;
IPT.var = noise_var;
IPT.F_tF = F_tF;
IPT.MIN_ITER = 10; IPT.MAX_ITER = 200; IPT.RELTOL=1e-5;
IPT.cut = 1e-16;

%%
IPT.alpha = 0.5; IPT.beta = 0.1;
[OPT_gam, history_gam]=solvefun_2d(IPT,0,'Gamma');
x_gamma = idct2(OPT_gam.C);
figure;imshow(x_gamma,[0,1]);title(['Gamma alpha = ', num2str(IPT.alpha),' beta = ', num2str(IPT.beta)]);

%%
IPT.alpha = 1.5; IPT.beta = 0.1;
[OPT_gam2, history_gam2]=solvefun_2d(IPT,0,'Gamma');
x_gamma2 = idct2(OPT_gam2.C);
figure;imshow(x_gamma2,[0,1]);title(['Gamma alpha = ', num2str(IPT.alpha),' beta = ', num2str(IPT.beta)]);

%%
IPT.alpha = 1; IPT.beta = 0.1;
[OPT_la, history_la]=solvefun_2d(IPT,0,'Gamma');
x_la = idct2(OPT_la.C);
figure;imshow(x_la,[0,1]);title(['Gamma alpha = ', num2str(IPT.alpha),' beta = ', num2str(IPT.beta)]);

%%
err_gam = OPT_gam.relerr(end);
err_gam2 = OPT_gam2.relerr(end);
err_la = OPT_la.relerr(end);
disp('gamma1 hpyerprior REL Err ');
disp(err_gam);
disp('gamma2 hpyerprior REL Err ');
disp(err_gam2);
disp('laplace hpyerprior REL Err ');
disp(err_la);

num_gam1=numel(find(abs(OPT_gam.C)<tol))/(m*n);
num_gam2=numel(find(abs(OPT_gam2.C)<tol))/(m*n);
num_la=numel(find(abs(OPT_la.C)<tol))/(m*n);

%%
cmin = -16.5;cmax =  5;

figure;imagesc(log(abs(C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);axis image off;
title('log|C_{X}| (True Coefficients)');
figure; imagesc(log(abs(dct2(Y)))); axis image; colormap jet;colorbar;caxis([cmin cmax]);axis image off;
title('log|C_{Y}| (Degraded Coefficients)'); 

figure; imagesc(log(abs(OPT_gam.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);axis image off;
title('log|C_{gamma1}| (\alpha=0.5 Restored Coefficients)');
figure; imagesc(log(abs(OPT_gam2.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);axis image off;
title('log|C_{gamma2}| (\alpha=1.5 Restored Coefficients)');
figure; imagesc(log(abs(OPT_la.C))); axis image; colormap jet;colorbar;caxis([cmin cmax]);axis image off;
title('log|C_{laplace}| (\alpha=1 Restored Coefficients)');
