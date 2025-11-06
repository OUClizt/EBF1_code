%% solvefun_2d
% EBF solution for 2D image deconvolution problem
% INPUT: 
% IPT: X - original image matrix, Y - observed image matrix,
%        p - blur kernel matrix, var - noise variance,
%        F_tF - diagonal elements of the forward matrix F'F,
%        alpha, beta - Gamma hyperprior parameters,
%        p - lp hyperprior parameter,
%        cut - small value pruning threshold.
% QUIET - whether to suppress the iterative relative change output,
% METHOD - the used method, options: 'Gamma', 'SBL', 'Gauss', 'lp', 'invGamma'.
%
% OUTPUT: C, gamma, rel_err, obj
%         history - records relative change of C
function [OPT, history] = solvefun_2d(IPT, QUIET, METHOD)
    t_start = tic; 
    
   %% Data preprocessing
    MIN_ITER = IPT.MIN_ITER; 
    MAX_ITER = IPT.MAX_ITER; 
    RELTOL   = IPT.RELTOL;
    
    X = IPT.X; % original image
    Y = IPT.Y; % degraded image
    
    ker = IPT.ker; % blur kernel
    var = IPT.var; % noise variance
    A = IPT.F_tF/var;
    
    % Preprocessing for Gamma hyperprior parameters
    if ~isfield(IPT,'alpha') || isempty(IPT.alpha)
    IPT.alpha = 1;
    end
    if ~isfield(IPT,'beta') || isempty(IPT.beta)
    IPT.beta = 1;
    end
    alpha = IPT.alpha-1;
    beta = IPT.beta;
    beta2 = 2/beta;
    % Function handles
    F = @(x)imfilter(idct2(x),ker,'symmetric');
    Ft = @(x)dct2(imfilter(x,ker,'symmetric'));

    [m,n] = size(X);    
    gamma_new = abs(dct2(Y)); % initial value
    c_new = ones(m,n);    
    

    if ~QUIET
        fprintf('%3s\t%10s\t%10s\t%10s\n', ... 
            'iter', 'rel error', 'rel delta', 'rel tol');
    end

    if nargin < 3 || isempty(METHOD)
        METHOD = 'Gamma';
    end

    %% Iterate between the update steps until convergence or max iterations           
    Flag=true(m,n);
    temp = zeros(m,n);
    for counter = 1:MAX_ITER
         gamma_old = gamma_new;
         c_old = c_new;
         % Update x 
        Fty=Ft(Y); % F'y
        
        temp(Flag) = (1./(1./(gamma_old(Flag))+A(Flag))).*Fty(Flag);
        invSy = Y/var-F(temp)/(var^2);  % invS*y
        c_new = gamma_old.* Ft(invSy);

        
        % Update gamma
        gamma_active = gamma_old(Flag);
        c_active = c_new(Flag);
        grad = A(Flag)./(1+gamma_active.*A(Flag));
        
%-------Gamma hyperprior-----------         
        if strcmpi(METHOD, 'Gamma')
            if alpha>=0
                numerator = alpha + sqrt(alpha^2 + (grad + beta2) .* (c_active.^2));
                gamma_new_active = numerator ./ (grad + beta2);
            else
                gamma_new_active = sqrt((c_active.^2)./(grad+beta2-2*alpha./(gamma_active)));
            end
        
%-------no hyperprior-----------         
        elseif strcmpi(METHOD, 'SBL')
        numerator = sqrt(grad .* (c_active.^2));
        gamma_new_active = numerator ./ (grad+1e-8);
        
%-------Gauss hyperprior-----------         
        elseif strcmpi(METHOD, 'Gauss')
        gamma_new_active=zeros(size(c_active));
        for i = 1:length(c_active)
        df = @(x) (2/beta)*x^3+ grad(i)*x^2 -(c_active(i)^2);
        x0 = [0,1000];
        gamma_new_active(i) = fzero(df, x0);       
        end
        
%-------General Gaussian hyperprior (lp)-----------         
        elseif strcmpi(METHOD, 'lp')
        p = IPT.p;
        gamma_new_active = sqrt((c_active.^2)./(grad+beta2*p*(gamma_active).^(p-1)));

%-------inverse Gamma hyperprior-----------
        elseif strcmpi(METHOD, 'invGamma')
        numerator = -alpha+sqrt(alpha^2+grad.*(c_active.^2+2*beta));
        gamma_new_active = numerator ./ grad;
        end
    
        
        gamma_new = zeros(m,n);
        gamma_new(Flag) = gamma_new_active;
        
        
        Sigma = 1./(A+1./gamma_new);
        OPT.Sigma = Sigma;
        
        if counter > MIN_ITER/2  % Start pruning in the second half of iterations
        Flag(gamma_new < IPT.cut) = false; % Dynamic pruning: small gamma excluded from iteration
        gamma_new(~Flag) = 0;
        end
       % ==== Compute objective function value ====
        switch lower(METHOD)
            case 'sbl'
                H = 0;
            case 'gamma'
                H = (gamma_new/ beta) - alpha * log(gamma_new + 1e-12);
            case 'gauss'
                H = (gamma_new.^2) / (2*beta);
            case 'lp'
                p = IPT.p;
                H = (gamma_new.^p) / beta;
            case 'invgamma'
                H = beta ./ (gamma_new + 1e-12) + alpha .* log(gamma_new + 1e-12);
            otherwise
                H = 0;
        end

        % Avoid log(0)
        A_gamma = A .* gamma_new;
        term1 = 0.5 * sum(log(1+A_gamma(:))); 
        term2 = 0.5 / var * norm(F(c_new)-Y,'fro')^2;
        term3 = 0.5 * sum((c_new(:).^2) ./ (gamma_new(:) + 1e-12));
        term4 = sum(H(:));
        OPT.obj(counter) = term1 + term2 + term3 + term4;
        
        % Record convergence
        history.rel_error(counter) =  norm(c_old-c_new,'fro' )/norm(c_old,'fro' ) ; 
        
        OPT.relerr(counter) = norm(X-idct2(c_new),'fro')/norm(X,'fro');
        
 
        if ~QUIET
            fprintf('%3d\t%0.2e\t%0.2e\t%0.2e\n', ... 
                counter, OPT.relerr(counter),history.rel_error(counter), RELTOL);
        end
        
        % Termination condition 
        if ( history.rel_error(counter) < RELTOL && ... 
            counter > MIN_ITER )
             break;
        end
  
    end

    % Output the total running time 
    if ~QUIET
        toc(t_start);
    end
    OPT.C = c_new;
    OPT.gamma = gamma_new;   
end
