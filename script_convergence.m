%Test the convergence results of the PALM alg.
%% Function input parameters
load('F_tF_mid','F_tF');
IPT.F_tF = F_tF;
IPT.X = X;IPT.Y = Y; 
IPT.ker = p;
IPT.var = noise_var;

IPT.MIN_ITER = 10; IPT.MAX_ITER = 200; IPT.RELTOL=1e-5;
IPT.cut = 1e-16;
 
%%
IPT.alpha = 1; 
IPT.beta = 10;[OPT_10, ~]=solvefun_2d(IPT,1,'Gamma');
IPT.beta = 1;[OPT_1, ~]=solvefun_2d(IPT,1,'Gamma');
IPT.beta = 0.1;[OPT_01, ~]=solvefun_2d(IPT,1,'Gamma');
IPT.beta = 0.01;[OPT_001, ~]=solvefun_2d(IPT,1,'Gamma');
%% Draw the convergence curve of the objective function
figure
semilogy(OPT_10.obj, 'LineWidth',2.0);hold on;
semilogy(OPT_1.obj,   'LineWidth',2.0);
semilogy(OPT_01.obj,  'LineWidth',2.0);
semilogy(OPT_001.obj, 'LineWidth',2.0);
hold off;
xlabel('Iteration','FontSize',20,'Interpreter','latex');
ylabel('Objective value','FontSize',20,'Interpreter','latex');

legend({'$\beta=10$', '$\beta=1$','$\beta=0.1$', '$\beta=0.01$'}, ...
       'Interpreter','latex','FontSize',18,'Location','northeast','Box','off');
set(gca,'FontSize',20,'LineWidth',1.2,'TickLabelInterpreter','latex');
axis tight; 
set(gca, 'LooseInset', get(gca, 'TightInset')); 
title('convergence comparison in \beta')