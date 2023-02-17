function [z,k,admm_optval] = solver_admml1(x0,z0,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL)
    n=size(A,2); %number of features
    m=size(A,1); %number of examples
    tic;
    x=x0;
    z=z0;
%     x = rand(n,1);
%     z = x;
    u = zeros(n,1);
    [L, U] = factor(A, rho);
    for k = 1:MAX_ITER
        % x-update
        q = Atb + rho*(z - u);
        if m >= n
           x = U \ (L \ q);
        else
           x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
        end
        % z-update
        zold = z;
        z = prox_l1(x + u, lambda*gamma);
        % u-update
        u = u + x - z;
        % diagnostics, reporting, termination checks
        admm_optval(k)   = objective(A, b, gamma, x, z);
        r_norm(k)   = norm(x - z);
        s_norm(k)   = norm(-rho*(z - zold));
        eps_pri(k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
        eps_dual(k) = sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
        if r_norm(k) < eps_pri(k) && s_norm(k) < eps_dual(k)
             break;
        end
    end
    x_admm = z;
    p_admm = admm_optval(end);
    admm_toc = toc;
    
end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if m >= n
       L = chol(A'*A + rho*speye(n), 'lower');
    else
       L = chol(speye(m) + 1/rho*(A*A'), 'lower');
    end
    L = sparse(L);
    U = sparse(L');
end

% function [ vX ] = prox_l1( vX, lambdaFactor )
% 
% % Soft Thresholding
% vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);
% end

function p = objective(A, b, gamma, x, z)
    %p = 0.5*sum_square(A*x - b) + gamma*norm(z,1);
    p = sum_square(A*x - b) + gamma*norm(z,1);
end