function [x_admm,k,admm_optval] = solver_xadmml1(x0,z0,mA,mb,paramLambda,paramL,SolverParams)
arguments
    x0 (:, 1) {mustBeFloat, mustBeReal}
    z0 (:, 1) {mustBeFloat, mustBeReal}
    mA (:, :) {mustBeFloat, mustBeReal}
    mb (:, 1) {mustBeFloat, mustBeReal}
    paramLambda (1, 1) {mustBeFloat, mustBeReal, mustBePositive} = 1
    paramL (:, :) {mustBeFloat, mustBeReal} = eye(size(mA,1)) %weight matrix
    SolverParams.paramRho (1, 1) {mustBeFloat, mustBeReal, mustBePositive} = 1 %1/paramLambda
    SolverParams.paramZeta (1, 1) {mustBeFloat, mustBeReal, mustBePositive} = 1 %in prox_g
    SolverParams.MAX_ITER (1, 1) {mustBeNumeric, mustBeReal, mustBePositive,mustBeInteger} = 1000
    SolverParams.ABSTOL (1, 1) {mustBeNumeric, mustBeReal, mustBePositive} = eps
    SolverParams.RELTOL (1, 1) {mustBeNumeric, mustBeReal, mustBePositive} = eps
    SolverParams.precision (1, 1) {mustBeNumeric, mustBeReal, mustBePositive} = eps
end
    %paramL = SolverParams.paramL;
    paramRho = SolverParams.paramRho;
    paramZeta = SolverParams.paramZeta;
    MAX_ITER = SolverParams.MAX_ITER;
    ABSTOL = SolverParams.ABSTOL;
    RELTOL = SolverParams.RELTOL;
    precision = SolverParams.precision;
    n = size(mA,2);
    A = 1;
    B = -1;
    tic;
    %Initialization
    x = x0;
    %x = zeros(n,1);
    X = zeros(n,MAX_ITER);
    z = z0;
    %z = zeros(n,1);
    Z = zeros(n,MAX_ITER);
    u = x-z;
    %u = zeros(n,1);
    U = zeros(n,MAX_ITER);
    
    X(:,1) = x;
    Z(:,1) = z;
    U(:,1) = u;
    admm_optval(1)   = objective(mA, mb, paramLambda, x, z);
%     [L, U] = factor(A, rho);
    for k = 2:MAX_ITER
        Gamma_1 = x-(A'*paramL*A*x + A'*paramL*(B*z+u))/paramLambda;
        %x-uodate
        [x,slvtol1] = prox_g(Gamma_1,mA,mb,paramZeta,precision);
        %z-update
        Gamma_2= z-(B'*paramL*B*z + B'*paramL*(A*x+u))/paramLambda;
        zold = z;
        %z = prox_l1(Gamma_2, paramLambda/ paramRho);
        z = prox_l1(Gamma_2,paramLambda/ paramRho);
        u = u + A * x + B * z;
%         % x-update
%         q = Atb + rho*(z - u);
%         if m >= n
%            x = U \ (L \ q);
%         else
%            x = lambda*(q - lambda*(A'*(U \ ( L \ (A*q) ))));
%         end
%         % z-update
%         zold = z;
%         z = prox_l1(x + u, paramLambda/paramZeta);
%         % u-update
%         u = u + x - z;
        %historical data
        X(:,k) = x;
        Z(:,k) = z;
        U(:,k) = u;
%        SlvTol1(k) = slvtol1;
%        SlvTol2(k) = slvtol2;
        % diagnostics, reporting, termination checks
        admm_optval(k)   = objective(mA, mb, paramLambda, x, z);
        r_norm(k)   = norm(x - z);
        s_norm(k)   = norm(-paramRho*(z - zold));
        eps_pri(k)  = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
        eps_dual(k) = sqrt(n)*ABSTOL + RELTOL*norm(paramRho*u);
        if r_norm(k) < eps_pri(k) && s_norm(k) < eps_dual(k)
             break;
        end
    end
    x_admm = z;
    p_admm = admm_optval(end);
    admm_toc = toc;
end

function [ vX ] = prox_l1( vX, lambdaFactor )

% Soft Thresholding
vX = max(vX - lambdaFactor, 0) + min(vX + lambdaFactor, 0);
end

function p = objective(mA, mb, lambda, x, z)
    p = 0.5*sum_square(mA*x - mb) + lambda*norm(z,1);
end

function [x,slvtol,cvxp] = prox_g(Gamma,mA,mb,paramZeta,precision)
    n = max(size(Gamma));
    cvx_begin
        %Return three precision values as follows:
        %solver prec. < solved prec. < Inaccurate/Solved prec.
        cvxp = cvx_precision(precision); 
        variable vx(n);
        %minimize(0.5 * sum_square(mA * vx - mb)+(0.5/paramZeta)*norm(vx- Gamma,2));
        minimize(sum_square(mA * vx - mb)+(1/paramZeta)*norm(vx- Gamma,2));
    cvx_end
    slvtol = cvx_slvtol; %Save the tolerance level the solver has achieved
    x = vx;
end

function [z,slvtol,cvxp] = xprox_l1(v,zeta,precision)
%   prox_l1(v,lambda) is the proximal operator of the l1 norm
%   with parameter lambda.
    n = max(size(v));
    cvx_begin
        %Return three precision values as follows:
        %solver prec. < solved prec. < Inaccurate/Solved prec.
        cvxp = cvx_precision(precision); 
        variable z(n);
%                minimize(norm(z,1)+(0.5/zeta)*norm(z- v,2));
        minimize(sum(abs(z))+(0.5/zeta)*norm(z- v,2));
    cvx_end
    slvtol = cvx_slvtol; %Save the tolerance level the solver has achieved
end