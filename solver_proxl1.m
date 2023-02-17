%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%% 
% This is an inexact Proximal-Gradient solver.
function [u_prox,F_prox_subopt,e1,e2,u0,k] = solver_proxl1(x0,A, b, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon0)
    %Some constant parameters
    tk = 1/L;
    beta = 0.5;
    %Cached computations for the gradient update step
    AtA = A'*A;
    Atb = A'*b;
%   minimize(sum_square(A * u - b)+gamma*l1_norm(u))
    f = @(u) sumsqr(A*u-b);
    objective = @(u) sumsqr(A*u-b) + gamma*l1_norm(u);
    
%     u0 = sprandn(n,N,0.05);
%     x0 = sprandn(n*N,1,0.05);
    x = x0;
    xprev = x;

    F_prox_subopt = zeros(MAX_ITER,1);
    e1 = zeros(n*N,MAX_ITER);
    e2 = zeros(MAX_ITER,1);
    for k = 1:MAX_ITER
        grad_x = AtA*x - Atb;
        %e1(:,k) = grad_x.*vpa(trandn(-delta,delta));
        %e2(k) = vpa(trandn(0,epsilon0));
        grad_x = grad_x; %+ e1(:,k);
        while 1
            %z = prox_l1(x - tk*grad_x, tk*1/lambda);
            z = prox_l1(x - tk*grad_x, lambda);%+sqrt(2*e2(k)/(L*n*N));
            if f(z) <= f(x) + grad_x'*(z - x) + (1/(2*tk))*sumsqr(z - x)
                break;
            end
            tk = beta*tk;
        end
        xprev = x;
        x = z;
        F_prox_subopt(k) = objective(x);
        F_prox_subopt(1) = objective(x0);
        if k > 1 && abs(F_prox_subopt(k) - F_prox_subopt(k-1)) < ABSTOL
            break;
        end
    end
    u0 = x0;
    u_prox = reshape(x,n,N);
end
%%Example: [u_prox,F_prox_subopt,e1,e2,u0,k] = lassompc_iprox(A_hat, b_hat, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,10^10*eps,10^10*eps)