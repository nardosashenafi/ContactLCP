% Inverted Pendulum
clc; clear all; close all

Vmax = 10;
T = 0.001;     % Sample period
theta_0 = 1*pi/180;  %Initial pendulum angle
x_0 = 0*0.1;   % Initial cart position
rg = 3.71;     % Gear ratio n2/n1
Kt = 0.00767;  % Motor torque constant Nm/A
rm = 6.35e-3;  % Radius of the powered wheel
M = 0.57;      % M = Mc2+Mw = 0.57 + 0.37
m = 0.165;      % Mass of the pendulum bar
g = 9.81;      % m/sec^2 acceleration of gravity
Lp = 0.6413/2; % length of pendulum rod is 0.6413               
R = 2.6;       % Armature resistance in Ohms
% Jp = 7.88e-3;  % Moment of inertia of the pendulum about its COM
Jp = m*Lp*Lp/3; %
r_enc = 0.01483;        % Encoder wheel radius in meters
Nenc = 4096;   % Number of counts per rev
Kenc = 2.275e-5; % Kenc = (2*pi/4096)*r_enc or r_enc = Kenc*4096/2/pi=0.0148
% r_enc = Kenc*4096/2/pi;
% Kenc = (2*pi/Nenc)*r_enc;
Fmax = Vmax*rg*Kt/R/rm;
Kemf = 1/R*(rg*Kt/rm)^2;

% Sample Frequency in Hz
fs = 1/T;
% Set Butterworth filter
fsc = (1/2)*fs; % Nyquist frequency = (1/2)*sample_frequency;
wn = 20/fsc;    % Low pass cutoff frequency is 20 Hz  
                % wn is the cutoff freq divided by the Nyquist freq
[bf,af] = butter(2,wn);  % Coefficients of the Butterworth filter
% xdot = filter(bf,af,xdot1); % Filtered speed
% acc = filter(bf,af,acc1);   % Filtered acceleration
% den = Jp*(M+m)+M*m*Lp*Lp;
kappa = 1/(Jp*(M+m)+M*m*Lp*Lp);
% beta_0 = m*g*Lp/(Jp*(M+m)+M*m*Lp*Lp);
% beta_2 = Jp/(Jp*(M+m)+M*m*Lp*Lp);

alpha3 = (Kt^2*rg^2/(R*rm^2))*kappa*(Jp+m*Lp^2);
alpha2 = -kappa*m*g*Lp*(M+m);
alpha1 = -(Kt^2*rg^2/(R*rm^2))*kappa*m*g*Lp;
alpha0 = 0;

a22 = -(rg^2*Kt^2/(rm^2*R))*kappa*(Jp+m*Lp^2);
a23 = -kappa*m*g*m*Lp^2;
a42 = (rg^2*Kt^2/(rm^2*R))*kappa*m*Lp;
a43 = kappa*m*g*Lp*(M + m);

b2 = (rg*Kt/(R*rm))*kappa*(Jp+m*Lp^2);
b4 = -(rg*Kt/(R*rm))*kappa*m*Lp;

A = [0 1 0 0; 0 a22 a23 0; 0 0 0 1; 0 a42 a43 0];
b = [0; b2; 0; b4];
C = [b A*b A^2*b A^3*b];

%LQR 
% Q = [.4 0 0 0; 0 0.45 0 0; 0 0 0 0; 0 0 0 0];
% r = .0002;
% [K_lqr,P,E]=lqr(A,b,Q,r);
% %K_lqr
% eig(A-b*K_lqr)

Ac = (inv(C)*A*C)'; % See text.
%Ac = [0 1 0 0; 0 0 1 0; 0 0 0 1; -alpha0 -alpha1 -alpha2 -alpha3];
bc = [0; 0; 0; 1];
Cc = [bc Ac*bc Ac*Ac*bc Ac*Ac*Ac*bc];
alpha0 = -Ac(4,1); alpha1 = -Ac(4,2); alpha2 = -Ac(4,3); alpha3 = -Ac(4,4);

% r1 = -(-16 + 14i);
% r2 = -(-16 - 14i);
% r3 = -(-1.4 + 1.2i);
% r4 = -(-1.4 - 1.2i);
r1 = 5; r2 = r1; r3 = r1; r4 = r1;

alphad3 = r1 + r2 + r3 + r4;
alphad2 = r1*r2 + r1*r3 + r1*r4 + r2*r3 + r2*r4 + r3*r4;
alphad1 = r1*r2*r3 + r1*r2*r4 + + r1*r3*r4 + r2*r3*r4;
alphad0 = r1*r2*r3*r4;

% Ac = [0 1 0 0; 0 0 1 0; 0 0 0 1; -alpha0 -alpha1 -alpha2 -alpha3];
% bc = [0; 0; 0; 1];
Kc = [alphad0-alpha0 alphad1-alpha1 alphad2-alpha2 alphad3-alpha3];
K = Kc*Cc*inv(C);
eig(Ac-bc*Kc);
eig(A-b*K)
K1 = K(1,1); K2 = K(1,2); K3 = K(1,3); K4 = K(1,4);

% Set the gains using Ackermann's Formula

% K_acker = [0 0 0 1]*inv(C)*(A^4 + alphad3*A^3 + alphad2*A^2 +...
%     alphad1*A + alphad0*eye(4));
% 

% format long 
% eig(A-b*K_acker);




