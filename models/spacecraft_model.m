clear all

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
[Ap,Bp,Cp,Dp] = ssdata(sys);
%Q=C'*C;
[A,B,C] = ss_augment(Ap,Bp,Cp,Dp); %trajectory following
%sys=ss(A,B,C,Dp); %trajectroy following
x0 = [-0.05;0.15;-0.08;300;-25;60;90]; %[rad/s X 4,rad/s X 1,deg X 3]
%spt = [0; 0; 0; 0; 0; 0; 0; 0];
x0_err = [0.0035;0.0052;0.0035;0.5;0.1;0.1;0.1]; % case II (noise)