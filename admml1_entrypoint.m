function [z,k,admm_optval] = admml1_entrypoint(x0,z0,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL,dt)
    T = mixedTypes(dt);  
    %input
    A=cast(A, 'like', T.x);
    b=cast(b, 'like', T.x);
    AtA=cast(A'*A, 'like', T.x);
    Atb=cast(Atb, 'like', T.x);
    %parameters
    lambda=cast(lambda, 'like', T.param);
    rho = 1/lambda;
    gamma=cast(gamma, 'like', T.param);
    %beta=cast(beta, 'like', T.param);
    
    [z,k,admm_optval] = solver_admml1(x0,z0,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL,T);
    z = double(z);
    k = round(k);
    admm_optval = round(admm_optval);
end