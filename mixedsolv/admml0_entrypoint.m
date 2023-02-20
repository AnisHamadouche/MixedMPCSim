function [z,admm_optval] = admml0_entrypoint(x0,z0,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL,frmt)
    T = mixedTypes(frmt.data);  
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
    
    T = mixedTypes(frmt.dt);
    x0 = cast(x0, 'like', T.x);
    z0 = cast(z0, 'like', T.x);
    u0 = zeros(size(x0,1),1, 'like', T.x);
    q = zeros(size(x0,1),1, 'like', T.x);
    
    [z,admm_optval] = solver_admml0(x0,z0,u0,q,A, b, Atb, lambda, gamma, rho, MAX_ITER,ABSTOL, RELTOL);
    z = double(z);
    admm_optval = round(admm_optval);
end