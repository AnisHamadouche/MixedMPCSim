function [dU,F_prox_subopt] = pgdl1_entrypoint(x0,A, b, Atb, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,frmt)
    T = mixedTypes(frmt.data);  
    %input
    A=cast(A, 'like', T.x);
    b=cast(b, 'like', T.x);
    AtA=cast(A'*A, 'like', T.x);
    Atb=cast(Atb, 'like', T.x);
    %parameters
    lambda=cast(lambda, 'like', T.param);
    L = cast(L, 'like', T.param);
    gamma=cast(gamma, 'like', T.param);
    %beta=cast(beta, 'like', T.param);
    
    T = mixedTypes(frmt.dt); 
    x0 = cast(x0, 'like', T.x);
    grad_x0 =zeros(size(x0), 'like', T.x);
    
    [dU,F_prox_subopt] = solver_proxl1(x0, grad_x0, A, b,AtA,Atb, n, N, gamma,L,MAX_ITER,ABSTOL,lambda);
    dU = double(dU);
    F_prox_subopt = round(F_prox_subopt);
end