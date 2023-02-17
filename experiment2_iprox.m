%%PG solver for a normal application with MPC length of N = 10 and 300 max
%%iterations
clc
clear all

delta = 2.2E-1; %MAX tolerable rounding error (in gradient)
epsilon_0 = 1E-4; %MAX tolerable prox error
N = 30; %MPC horizon length
ABSTOL = eps; %Tolerance
MAX_ITER = 20; % gradient descent maximum iteration number
%digits(64);% for vpa(.)

[sys, H, F, x0, m, n, l, xmax, xmin, ymax, ymin] = model_spacecraft(N);

x_cvx = x0;
G = F * x_cvx;
%A_hat = (1/sqrt(2)) * H^(1/2);
A_hat = (1/sqrt(2)) * sqrtm(H);
%b_hat = (-1/sqrt(2)) * (H^(-1/2) * G);
b_hat = (-1/sqrt(2)) * (inv(sqrtm(H)) * G);
L = max(eigs(A_hat'*A_hat));
lambda = 1/L;
gamma_max = norm(A_hat' * b_hat,'inf');
gamma = 0.01*gamma_max;

objective = @(u) sumsqr(A_hat*u-b_hat) + gamma*l1_norm(u);

%Simulation
N_sim=400;

%Initialization
xm=x0;
u=0;
y=sys.C*x0;

for kk=1:N_sim
%     eta=-(Omega\Psi)*Xf;
%     deltau=L0'*eta;
%     beta=E*Xf+M_act*eta;
%     if (beta<y_bar_min)
%         M_act1=-M_act;
%         lambda=-inv(M_act1*inv(Omega)*M_act1')*(-y_bar_min+E*Xf+M_act1*(-eta));
%         eta=-inv(Omega)*(Psi*Xf+M_act1'*lambda);
%         deltau=L0'*eta;
%     end
%     if (beta>y_bar_max)
%         lambda=-inv(M_act*inv(Omega)*M_act')*(y_bar_max-E*Xf+M_act*(-eta));
%         eta=-inv(Omega)*(Psi*Xf+M_act'*lambda);
%         deltau=L0'*eta; 
%     end
%     if(deltau>deltau_max) 
%         deltau=deltau_max;
%     end
%     if (deltau<deltau_min)
%         deltau=deltau_min;
%     end
%     u=u+deltau;
%     if (u>u_max) 
%         deltau=u_max-up; 
%         u=u_max;
%     end
%     if (u<u_min) 
%         deltau=u_min-up; 
%         u=u_min; 
%     end
    deltau=solver_iproxl1(A_hat, b_hat, n, N, gamma,L,MAX_ITER,ABSTOL,lambda,delta,epsilon_0);
    deltau1(1,kk)=deltau;
    u1(1,kk)=u;
    y1(1,kk)=y;
    %%%%
    %plant simulation
    %%%%%%
    yp=y;
    xm_old=xm;
    xm=sys.A*xm+sys.B*(u+d(kk)); % calculate xm(k+1)
    y=sys.C*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;(y-sp(:,kk+1))];
    %y_bar_min=y_min-sp(1,kk+1);
    %y_bar_max=y_max-sp(1,kk+1);
    up=u;
end
dU1=deltau1;
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