[A,B,C]=ss_augment(Ac,Bc,Cc,Dc);
sysaugmented = ss(A,B,C,Dc);
Q  = diag([500,500,500,1.0e-7,1,1,1]);
R  = diag([200,200,200,1]);
m = size(A,1);            % Dimension of state
n = size(B,2);            % Dimension of input control
l = size(C,1);            % Dimension of output
K_lqr = lqr(sysaugmented,diag([diag(Q)',diag(Q)']),R);
N_sim=200;
xm = x0;
x=xm;
u=0;
y=Cc*x0;
y_init=y;
%spt = zeros(7,N_sim+10);
spt=[zeros(7,120) ones(7,120) 100*ones(7,200+N_sim)]*0;
d=[zeros(4,50) 100*ones(4,20) zeros(4,200+170+N_sim)]*1;%disturbance
dUp=[];
Xf=[zeros(size(xm,1),1)];
dU1=[];
u1=[];
y1=[];
clear iter
for kk=1:N_sim
    dU = -K_lqr*Xf;
    dUp=dU;
    %dU = reshape(dU,n,N); %uncomment for admm
    u=u+dU(:,1);
    iter(kk)=k;
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
    %Xf=[xm-xm_old; y-spt(:,kk+1)];
    Xf=[xm-xm_old];
    up=u;
    if all(abs(y-spt(:,kk+1)) < 0.02*abs(y_init-spt(:,kk+1)))
        t_f=kk;
        %N_sim=t_f;
        iter=[iter, zeros(1,N_sim-kk)];
        break
    end
end
k=0:(N_sim-1);
figure(1)
subplot(2,1,1)
% hold;
plot(k,y1')
xlabel('Sampling instant');
ylabel('y');
legend('Output');