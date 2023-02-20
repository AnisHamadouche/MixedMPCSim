function p = objective_l0(A, b, gamma, x, z)
    %p = 0.5*sum_square(A*x - b) + gamma*norm(z,1);
    p = sum_square(A*x - b) + gamma*sum(z(:)~=0);
end