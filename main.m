%MPC design parameters
% Q  = diag([500,500,500,1.0e-7,1,1,1]);
% R  = diag([200,200,200,1]);
a=0.2*ones(size(B,2),1);%tuning parameter 1
N=1*ones(size(B,2),1); %tuning parameter 2
Np=46;

[A,B,C] = ss_augment(A,B,C,D);

Q=C'*C;
R=0.1;

[Omega,Psi]=dmpclagr(A,B,a,N,Np,Q,R);


%Closed-loop performance evaluation
[A1,L0]=lagd(a(1),N(1));
K_mpc=L0'*(Omega\Psi);
Acl=A-B*K_mpc;
Pole_close=eig(Acl);
M_act=C*B*L0';
E=C*A;

%Constraints
u_min=-100;
u_max=100;
deltau_min=-100;
deltau_max=100;
y_min=-100;
y_max=100;

%Simulation
N_sim=400;

%Initialization
%xm=zeros(n1,1);
u=0;
%y=0;
y_delta_k=0;
y_delta_k_m1=0;
u_delta_k_m1=0;
%sp=[zeros(1,10) ones(1,N_sim)];
sp=[ones(1,120) -zeros(1,120) ones(1,200)]*0;
d=[ones(1,120) -zeros(1,120) ones(1,200)]*0;
%Xf=[y_delta_k; y_delta_k_m1;u_delta_k_m1;y-r(1,1)];
Xf=[xm;(y-sp(:,1))];
%y_bar_min=y_min-sp(1,1);
%y_bar_max=y_max-sp(1,1);