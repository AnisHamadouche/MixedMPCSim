%%PG solver for a normal application with MPC length of N = 10 and 300 max
%%iterations
%%
%Solver parameters
delta = 2.2E-1; %MAX tolerable rounding error (in gradient)
epsilon_0 = 1E-4; %MAX tolerable prox error
N = 4; %MPC horizon length
ABSTOL = eps; %Absolute Tolerance
MAX_ITER = 300; % gradient descent maximum iteration number
%Extra ADMM solver parameters
RELTOL = eps; %Relative Tolerance
lambda=1;
rho=1/lambda;
%%
%Q=Cp'*Cp;
%R=0.01*eye(size(B,2));
Q  = diag([500,500,500,1.0e-7,1,1,1]);
R  = diag([200,200,200,1]);
PSI = zeros(size(A,1)*N,size(A,2));
m = size(A,1);            % Dimension of state
n = size(B,2);            % Dimension of input control
l = size(C,1);            % Dimension of outpu
N_sim=200;
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
% Q_bar = dlyap((A-B*K_lqr)', Q+K_lqr'*R*K_lqr);
% Q_hat((i-1)*m+1:i*m,(i-1)*m+1:i*m) = Q_bar;
Q_hat = diag([diag(Q_hat)',diag(Q_hat)']);
H = THETA'*Q_hat*THETA + R_hat;
F = THETA'*Q_hat*PSI;
objective = @(u) sumsqr(A_hat*u-b_hat) + gamma*l1_norm(u);
%Initialization
xm = x0;
u=0;
y=Cp*x0;
%spt = zeros(7,N_sim+10);
spt=[zeros(7,120) 100*ones(7,120) 100*ones(7,200)]*0;
d=[zeros(4,50) 50*ones(4,20) zeros(4,200+170)];%disturbance
%Xf=[zeros(size(xm,1),1)]; %regulator
Xf=[zeros(size(xm,1),1);y-spt(:,1)]; %trajectory following
G = F * Xf;
A_hat = (1/sqrt(2)) * sqrtm(H);
b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
L = max(eigs(A_hat'*A_hat));
lambda = 1/L;
gamma_max = norm(A_hat' * b_hat,'inf');
gamma = 0.01*gamma_max;
dUp=sprand(N*n,1,0.01);
dU1=[];
u1=[];
y1=[];
Atb = A_hat'*b_hat;
THRESHOLD=eps;
for kk=1:N_sim
    if abs(y-spt(:,kk+1))>= THRESHOLD
        G = F * Xf;
        b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
        %gamma_max = norm(A_hat' * b_hat,'inf');
        %gamma = 0.001*gamma_max;
        %Run PG
        [dU,F_prox_subopt,e1,e2,u0,k]=solver_proxl0(dUp,A_hat, b_hat, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon_0);
    %     [dU,k,admm_optval] = solver_admml1(dUp,dUp,A_hat, b_hat, Atb, lambda, gamma, rho, MAX_ITER, ABSTOL, RELTOL);
    %     dUp=dU;
    %     dU = reshape(dU,n,N);
        u=u+dU(:,1);
        iter(kk)=k;
    else
        dU=zeros(4,1);
        u=0;
        iter(kk)=0;
    end
    dU1(:,kk)=dU(:,1);
    u1(:,kk)=u;
    y1(:,kk)=y;
    %%%%
    %plant simulation
    %%%%%%
    yp=y;
    xm_old=xm;
    xm=Ap*xm+Bp*(u+d(:,kk)); % calculate xm(k+1), %add input disturbance u+d(kk)
    y=Cp*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old; y-spt(:,kk+1)];
    %Xf=[xm-xm_old];
    up=u;
end
k=0:(N_sim-1);
figure(1)
subplot(2,1,1)
% hold;
plot(k,y1')
xlabel('Sampling instant');
ylabel('y');
legend('Output');
subplot(2,1,2)
% hold;
stairs(k,u1');
xlabel('Sampling instant');
ylabel('u');
legend('control');
figure(2)
% hold;
plot(k,dU1)
ylabel('\Delta u')
xlabel('Sampling Instant');

