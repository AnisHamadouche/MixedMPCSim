%%PG solver for a normal application with MPC length of N = 10 and 300 max
%%iterations
%%
clear all
clc

frmt.dt = 'fixed28'; % solver data format: 'single', 'fixed16', 'fixed12', 'fixed8'
frmt.data = 'double'; % data casting format: edit ./utils/mixedTypes
para.problem = 'l1'; % 'l1', 'l0'
para.solver = 'admm'; % 'pgd', 'admm'

addpath './utils'
addpath './models'
addpath './solver'
addpath './mixedsolv'

%% Read model and constraints
spacecraft_model
spacecraft_constraints

%% Solver parameters
N = 10; %MPC horizon length
MAX_ITER = 100; % gradient descent maximum iteration number
%Extra ADMM solver parameters
RELTOL = 1-2; %Relative Tolerance
ABSTOL = 1e-4; %Absolute Tolerance
lambda=1;
rho=1/lambda;
%paramL=eye(size(A_hat));
%%
%Q=Cp'*Cp;
%R=0.01*eye(size(B,2));
Q  = diag([500,500,500,1.0e-7,1,1,1]);
R  = diag([200,200,200,1]);
PSI = zeros(size(A,1)*N,size(A,2));
m = size(A,1);            % Dimension of state
n = size(B,2);            % Dimension of input control
l = size(C,1);            % Dimension of output
N_sim=100;
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
%sysaugmented = ss(A,B,C,Dc);
%K_lqr = lqr(sysaugmented,diag([diag(Q)',diag(Q)']),R);
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
%objective = @(u) sumsqr(A_hat*u-b_hat) + gamma*l1_norm(u);
% Initialization
xm = x0;
u=0;
y=Cp*x0;
y_init=y;

%% Read setpoint and disturbance signals
signals
%%

% Xf=[zeros(size(xm,1),1)]; %regulator
Xf=[zeros(size(xm,1),1);y-spt(:,1)]; %trajectory following
G = F * Xf;
A_hat = (1/sqrt(2)) * sqrtm(H);
b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
L = max(eigs(A_hat'*A_hat));
lambda = 1/L;
gamma_max = norm(A_hat' * b_hat,'inf');
gamma = 0.01*gamma_max;
%dUp=zeros(N*n,1);
dUp=100*rand(N*n,1);
dU1=[];
u1=[];
y1=[];
clear iter
clear k
optval = [];
for kk=1:N_sim
    G = F * Xf;
    b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
    %gamma_max = norm(A_hat' * b_hat,'inf');
    %gamma = 0.001*gamma_max;
    Atb = A_hat'*b_hat;
    %Run Solver
    if strcmp(para.problem, 'l1')
        if strcmp(para.solver, 'admm')
            [dU,admm_optval] = admml1_entrypoint(dUp,dUp,A_hat, b_hat, Atb, lambda, gamma, rho, MAX_ITER, ABSTOL, RELTOL, frmt);
            optval = [optval, admm_optval(end)];
        elseif strcmp(para.solver, 'pgd')
            [dU,F_prox_subopt] = pgdl1_entrypoint(dUp,A_hat, b_hat, Atb, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,frmt);
            optval = [optval, F_prox_subopt(end)];
        end
    elseif strcmp(para.problem, 'l0')
        if strcmp(para.solver, 'admm')
            [dU,admm_optval] = admml0_entrypoint(dUp,dUp,A_hat, b_hat, Atb, lambda, gamma, rho, MAX_ITER, ABSTOL, RELTOL, frmt);
            optval = [optval, admm_optval(end)];
        elseif strcmp(para.solver, 'pgd')
            [dU,F_prox_subopt] = pgdl0_entrypoint(dUp,A_hat, b_hat, Atb, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,frmt);
            optval = [optval, F_prox_subopt(end)];
        end
    end
    dUp=dU; %uncomment for admm
    dU = reshape(dU,n,N); %uncomment for admm
    u=u+dU(:,1);
%     iter(kk)=k;
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
    if all(abs(y-spt(:,kk+1)) < 0.02*abs(y_init-spt(:,kk+1)))
        t_f=kk;
        %N_sim=t_f;
        iter=[iter, zeros(1,N_sim-kk)];
        break
    end
end

k=0:(kk-1);
figure(1)
subplot(2,1,1)
str_titlle = sprintf('Solver: %s, Type: %s, Signedness: %s,\n Word length: %d, Fraction length: %d',upper(para.solver), mixedTypes(frmt.dt).x.DataTypeMode, mixedTypes(frmt.dt).x.Signedness, mixedTypes(frmt.dt).x.WordLength, mixedTypes(frmt.dt).x.FractionLength);
% hold;
%semilogy(k,abs(y1'))
plot(k,y1(5:end,:)')
xlabel('Sampling instant','FontSize',16);
ylabel('y','FontSize',16);
legend('$Roll$','$Pitch$','$Yaw$', ...
       'Interpreter','latex');
sgtitle(str_titlle);
subplot(2,1,2)
% hold;
stairs(k,u1');
xlabel('Sampling instant','FontSize',16);
ylabel('u','FontSize',16);
legend('$\tau_1$','$\tau_2$','$\tau_3$','$\tau_\omega$', ...
       'Interpreter','latex');
figure(2)
% hold;
plot(k,dU1)
title(str_titlle);
ylabel('$\Delta u$','FontSize',16)
xlabel('Sampling Instant','FontSize',16);
legend('$\tau_1$','$\tau_2$','$\tau_3$','$\tau_\omega$', ...
       'Interpreter','latex');
% PEAK(alg,:)=max(abs(y1(2:end)'));
% ITER(alg,:) = iter;
% ERROR(alg,:) = sum_square(y1'-spt(1:kk)');
% SumdU(alg,:) = sum_square(dU1');
% SumU(alg,:) = sum_square(u1');
% SettTime(alg,:)=kk;
% MAXU(alg,:)=max(u1');
% MAXdU=max(dU1');
% MINU(alg,:)=min(u1');
% MINdU=min(dU1');
% ERRORF(alg,:)=norm(y-spt(N_sim));
% alg=alg+1;