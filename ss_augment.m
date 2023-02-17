function [Ae,Be,Ce] = ss_augment(Ap,Bp,Cp,Dp)
%ss_augment() forms a discrete-time augmented state space description
%from a continuous-time state space model by the following transormations:
%
% x => [\delta x, y]'
% u => \delta u 
%
%See: (1.5) in https://drive.google.com/file/d/1Pfn9IzV24xJ07utEBiJOAV5xY8Kj_0yv/view?usp=sharing
%
%Inputs:
%   Ac : n1 x n1 matrix;
%   Bc : n1 x n_in matrix;
%   Cc : m1 x n1 matrix;
%   Dc : 1 x 1 matrox;
%   Delta_t : discretization sample time;
%
%Outputs:
%   Ae : n1+m1 x n1+m1 matrix
%   Be : n1+m1 x n_in matrix
%   Ce : m1 x n1+m1 matrix
    %[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
    [m1,n1]=size(Cp);
    [n1,n_in]=size(Bp);
    Ae=eye(n1+m1,n1+m1);
    Ae(1:n1,1:n1)=Ap;
    Ae(n1+1:n1+m1,1:n1)=Cp*Ap;
    Be=zeros(n1+m1,n_in);
    Be(1:n1,:)=Bp;
    Be(n1+1:n1+m1,:)=Cp*Bp;
    Ce=zeros(m1,n1+m1);
    Ce(:,n1+1:n1+m1)=eye(m1,m1);
end