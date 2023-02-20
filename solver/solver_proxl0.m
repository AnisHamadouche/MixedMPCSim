%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%% 
% This is an inexact Proximal-Gradient solver.
function [u_prox,F_prox_subopt] = solver_proxl0(x0,grad_x0,A, b,AtA,Atb, n, N, gamma,L,MAX_ITER,ABSTOL,lambda)
    %Some constant parameters
    tk = 1/L;
    beta = 0.5;
%     x = cast(x0, 'like', T.x);
%     z=cast(x0, 'like', T.x);
%     grad_x=zeros(size(x0), 'like', T.x);
    x=x0;
    z=x0;
    grad_x = grad_x0;


    n=size(A,2); %number of features
    m=size(A,1); %number of examples
    tic;

%    f = @(u) 0.5*sumsqr(A*u-b)+gamma*l1_norm(u);
%     objective = @(u) sumsqr(A*u-b) + gamma*l1_norm(u);

%    F_prox_subopt = zeros(1,MAX_ITER,'like',T.x);
%    grad_x=zeros(size(x0), 'like', T.x);
    for k = 1:MAX_ITER
        grad_x(:) = AtA*x - Atb;
%         while 1
%        z = prox_l0(x - tk*grad_x, gamma);
         z = prox_l0(x - tk*grad_x, lambda*gamma);
%             if f(z) <= f(x) + grad_x'*(z - x) + (1/(2*tk))*sumsqr(z - x)
%                 break;
%             end
%             tk = beta*tk;
%         end
        xprev = x;
        x = z;
        F_prox_subopt(k) =  objective_l0(A, b, gamma, x, x);%f(x);
        F_prox_subopt(1) =  objective_l0(A, b, gamma, x0, x0);%f(z);
        if k > 1 && abs(F_prox_subopt(k) - F_prox_subopt(k-1)) < ABSTOL
            break;
        end
    end
    u0 = x0;
    u_prox = x;
end

function p = objective_l0(A, b, gamma, x, z)
    p = 0.5*sum_square(A*x - b) + gamma*sum(z(:)~=0);
end