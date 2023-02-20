function [z,admm_optval] = solver_admml1(x0,z0,u0,q,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL)
    n=size(A,2); %number of features
    m=size(A,1); %number of examples
    tic;
    
    %x0 = cast(x0, 'like', T.x);
    %z0 = cast(z0, 'like', T.x);
    x=x0;
    z=z0;
    u=u0;
    zold = z;
    %u = zeros(n,1, 'like', T.x);
    %q = zeros(n,1, 'like', T.x);

    [L, U] = factor(A, rho);
    U_inv = pinv(U);
    L_inv = pinv(L);
    admm_optval   = zeros(1,MAX_ITER);
    r_norm   = zeros(1,MAX_ITER);
    s_norm   = zeros(1,MAX_ITER);
    eps_pri  = zeros(1,MAX_ITER);
    eps_dual = zeros(1,MAX_ITER);    

    for k = 1:MAX_ITER
       % x-update
        q(:) = Atb + rho*(z - u);
        if m >= n
           %x = U \ (L \ q);
           x(:) = U_inv * (L_inv * q);
           %x_iter(:,k)=x;
        else
           %x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
           x(:) = q/rho - (A'*(U_inv*(L_inv*(A*q))))/rho^2;
           %x_iter(:,k)=x;
        end
        % z-update
        zold(:) = z;
        z(:) = prox_l1(x + u, lambda*gamma);
        %z_iter(:,k)=z;
        
        % u-update
        u(:) = u + x - z;
        %u_iter(:,k)=u;

        % diagnostics, reporting, termination checks
        admm_optval(k)   = objective(A, b, gamma, x, z);
        r_norm(k)   = sqrt(sumsqr(x - z));
        s_norm(k)   = sqrt(sumsqr(-rho*(z - zold)));
        eps_pri(k)  = sqrt(n)*ABSTOL + RELTOL*max(sqrt(sumsqr(x)), sqrt(sumsqr(-z)));
        eps_dual(k) = sqrt(n)*ABSTOL + RELTOL*sqrt(sumsqr(rho*u));
        if r_norm(k) < eps_pri(k) && s_norm(k) < eps_dual(k)
             break;
        end
    end
    x_admm = z;
    p_admm = admm_optval(end);
    admm_toc = toc;
end

function [L, U] = factor(A, rho)
    %A = double(A);
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( double(A'*A + rho*eye(n)), 'lower' );
    else            % if fat
       L = chol( double(eye(m) + 1/rho*(A*A')), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = L;
    U = L';
end

% function [ vX ] = prox_l1( vX, lambdaFactor )
% 
% % Soft Thresholding
% vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);
% end

function p = objective(A, b, gamma, x, z)
    p = 0.5*sum_square(A*x - b)+ gamma*l1_norm(z);
end