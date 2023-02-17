%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%% 
% This is an inexact accelerated Proximal-Gradient solver
function [u_prox,F_prox_subopt,e1,e2,u0,k] = solver_axproxl0(A, b, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon0)
    % some constant parameters
    tk = 1/L;
    beta = 0.5;
    %Cached computations for gradient step
    AtA = A'*A;
    Atb = A'*b;
    %minimize(sum_square(A * u - b)+gamma*l1_norm(u))
    f = @(u) sumsqr(A*u-b);
    objective = @(u) sumsqr(A*u-b) + gamma*l1_norm(u);
    
    %Initialize;
    x0 = sprandn(n*N,1,0.5);
    x = x0;
    xprev = x;

    F_prox_subopt = zeros(MAX_ITER,1);
    e1 = zeros(n*N,MAX_ITER);
    e2 = zeros(MAX_ITER,1);
    for k = 1:MAX_ITER
        alpha = (k+1)/2;
        %alpha = trandn(1,k);
        zeta = (alpha-1)/alpha;
        y = (1+zeta)*x-zeta*xprev;
        grad_y = AtA*y - Atb;
        e1(:,k) = grad_y.*vpa(trandn(-delta,delta));
        e2(k) = vpa(trandn(0,epsilon0));
        grad_y = grad_y + e1(:,k);
        while 1
            z = prox_l0(y - tk*grad_y, lambda)+sqrt(2*e2(k)/(L*n*N));
            if f(z) <= f(y) + grad_y'*(z - y) + (1/(2*tk))*sumsqr(z - y)
                break;
            end
            tk = beta*tk;
        end
        xprev = y;
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