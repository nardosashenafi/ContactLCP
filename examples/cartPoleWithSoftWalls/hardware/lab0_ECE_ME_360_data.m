clc; close all; clear

% Encoder Parameters
r_enc = 0.01483;        % Encoder wheel radius in meters
Nenc = 4096;
Kenc = r_enc*2*pi/Nenc;  % Kenc = 2.275*10^-5; % meters/count
T = 0.001; % seconds
speed_error = Kenc/T;
% Motor Parameters for Calculating a0 and b0
R = 2.6; KT = 7.67e-3; Kb = KT; rm = 6.35e-3;
rg = 3.71; M = 0.57; f = 0; J1 = 3.9e-7; J2 = 0; Jenc = 0;
J = Jenc*rg^2 + J2 + J1*rg^2;
% Motor Model   G(s) =  b0/(s(s+a0))
a0 = (R*f/rm^2 + Kb*KT*(rg/rm)^2)/(R*J/rm^2 + R*M);  % a0 = 10.9846
b0 = (KT*rg/rm)/(R*J/rm^2 + R*M);                    % b0 = 2.4513

