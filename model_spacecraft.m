% clear all
% num=[0.1];
% denum=[1 0 0 0 0 0 0 -0.8];
% [Ap,Bp,Cp,Dp]=tf2ss(num,denum);
% [n1,n_in]=size(Bp);
% [m1,n1]=size(Cp);
% A=eye(n1+m1,n1+m1);
% A(1:n1,1:n1)=Ap;
% A(n1+1:n1+m1,1:n1)=Cp*Ap;
% B=zeros(n1+m1,n_in);
% B(1:n1,:)=Bp;
% B(n1+1:n1+m1,:)=Cp*Bp;
% C=zeros(m1,n1+m1);
% C(:,n1+1:n1+m1)=eye(m1,m1);
%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Yun Wu

%% Model of the SSETI/ESEO satellite (

% x=f(x,u)=[ \dot(\omega)^b_ob, \dot(\omega)_w, \dot(\mu), \dot(\varepsion)]^T

%, where \dot(\omega)^b_ob=[\omega_1,\omega_2,\omega_3]^T is the angular velocity
%of the bodyframe relative to the orbit frame, \omega_w is the angular velocity of
%the wheels about their spin axes (scalar in our case),and \varepsion=[\varepsion_1,
%\varepsion_2,\varepsion_3]^T, which together with \mu, make up the Euler parameters.

%The control input: is given as u =[\tau^T,\tau^T_w]^T=[\tau_1,\tau_2,\tau_3,\tau_w]^T.

%Initial conditions:b by choosing the equilibrium pointp equal to x^p=[0 0 0 0 1 0 0 0]^T,
%u^p=[0 0 0 0],

%Parameter space (dropping \mu): −[1,1,1,800,1,1,1]^T \leq x \leq [1,1,1,800,1,1,1]^T

%Constraints: u_max = -u_min = [0.0484, 0.0484, 0.0398, 0.0020]^T, |\omega_w| \leq 527

%Note: when deriving the controller, \mu is omitted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ref: Øyvind Hegrenæs, J. T. Gravdahl, and P. Tøndel, “Spacecraft
%attitude control using explicit model predictive control,” Automatica,
%vol. 41, no. 12, pp. 2107 – 2114, 2005. [Online]. Available:
%http://www.sciencedirect.com/science/article/pii/S0005109805002657
%% 
function [sys, H, F, x0, m, n, l, xmax, xmin, ymax, ymin] = model_spacecraft(N)
    %% Unconstrained Lasso MPC problem:
    % min sum_1^N ||x(k)^T*Q(k)*x(k)+u(k)^T*R(k)*u(k)||_2^2+lamda*||u(k)||_1
    % s.t. x(k+1) = A*x(k) + B*u(k)
    Ts = 0.1;
    w0 = 1;
    i11 = 4.250;
    i22 = 4.337;
    i33 = 3.664;
    Iw = 4.0e-5;
    k = i22 - Iw;
    kx = (i22 - i33)/i11;
    ky = (i11 - i33)/i22;
    kz = (i22 - i11)/i33;

    % example of a simple LTI system
    Ac = [  0               0               (1-kx)*w0           0                            -8*kx*w0^2           0                   0;          ...
            0               0               0                   0                             0                  -6*ky*i22*w0^2/k     0;          ...
            (kz -1)*w0      0               0                   0                             0                   0                  -2*kz*w0^2;  ...
            0               0               0                   0                             0                  6*ky*i22*w0^2/k      0;          ...
            0.5             0               0                   0                             0                   0                   0;          ...
            0               0.5             0                   0                             0                   0                   0;          ...
            0               0               0.5                 0                             0                   0                   0];    
   
   Bc = [  1/i11           0               0                   0;          ...
            0               1/k             0                   -1/k;       ... 
            0               0               1/i33               0;          ...
            0               -1/k            0                   i22/(k*Iw);   ...
            0               0               0                   0;          ...
            0               0               0                   0;
            0               0               0                   0];
    Cc = eye(7); 
    Dc = zeros(7,4);
    sys = ss(Ac,Bc,Cc,Dc,Ts);
    [A,B,C,D] = ssdata(sys);
    m = size(A,1);            % Dimension of state
    n = size(B,2);            % Dimension of input control
    l = size(C,1);            % Dimension of output
    Q  = diag([500,500,500,1.0e-7,1,1,1]);
    R  = diag([200,200,200,1]);
    x0 = [-0.05;0.15;-0.08;300;-25;60;90]; %[rad/s X 4,rad/s X 1,deg X 3]
    spt = [0; 0; 0; 0; 0; 0; 0; 0];
    x0_err = [0.0035;0.0052;0.0035;0.5;0.1;0.1;0.1]; % case II (noise)
    umax=[0.0484, 0.0484, 0.0398, 0.0020];
    umin=-umax;
    xmax = [1; 1; 1; 800; 1; 1; 1]; 
    xmin = -xmax;
    ymax =  [   1;          1;         1;           527;        1;      1;      1   ];
    ymin = -[   1;          1;         1;           527;        1;      1;      1   ];
    %T = 20;
    constraint = false;

    %% QP to Lasso convert
    PSI = zeros(size(A,1)*N,size(A,2));
    for i=1:N
        PSI((i-1)*size(A,1)+1:i*size(A,1),1:size(A,2)) = A^i;
    end 
    THETA = zeros(size(A,1)*N,size(B,2)*N);
    for i=1:N
        for j=1:i
            THETA((i-1)*size(A,1)+1:i*size(A,1), ...
                  (j-1)*size(B,2)+1:j*size(B,2)) = A^(i-j)*B;
        end
    end
    K_lqr = lqry(sys,Q,R);
    Q_hat = zeros(size(Q,1)*N, size(Q,2)*N);
    R_hat = zeros(size(R,1)*N, size(R,2)*N);
    for i=1:N
        Q_hat((i-1)*size(Q,1)+1:i*size(Q,1), ...
              (i-1)*size(Q,2)+1:i*size(Q,2)) = Q;
        R_hat((i-1)*size(R,1)+1:i*size(R,1), ...
              (i-1)*size(R,2)+1:i*size(R,2)) = R;
    end
    Q_bar = dlyap((A-B*K_lqr)', Q+K_lqr'*R*K_lqr);
    Q_hat((i-1)*m+1:i*m,(i-1)*m+1:i*m) = Q_bar;
    H = THETA'*Q_hat*THETA + R_hat;
    F = THETA'*Q_hat*PSI;
end