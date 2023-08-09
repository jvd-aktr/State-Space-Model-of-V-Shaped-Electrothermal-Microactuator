%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Material Properties
rho=2330; %kg/m3
res=51e-6; %ohm.m
alpha=2.5e-6;
cp=732; %J/KgK
k=156; %w/mk
k_hat=k/(rho*cp);
E=160e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Geometrical
Wb=5e-6;
Lb=176e-6;
tb=20e-6;
ts=tb;
Ls=20e-6;
Ws=10e-6;
U=4;
L=(2*Lb)+Ls;
H=2;
theta=deg2rad(5.6);%deg
%V-Shaped formulas
B=(L^2)/(Wb^2*(cos(theta)^4)+(L^2*sin(theta)^2));
B1=Lb*B*sin(theta);
B2=(B*Lb)/(Wb*tb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steady State
A_=(2*Lb*Wb*tb)+(Ls*Ws*ts);
B_=((2*Lb*Ws*ts)+(Ls*Wb*tb));
C_=Ws*ts*Wb*ts;
% qdot_ss=(U^2)/(res*((B/C)*A)); 
% delT_w=((1/12)*((qdot_ss*L^2)/(rho*cp*k_hat)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derivative of system

 
A=(-k_hat*pi^2)/(L^2);
B=(pi^2/(12*rho*cp))*(1/(res*((B_/C_)*A_)));
C = 1;
D = 0;

figure(1);
sys = ss(A,B,C,D)
% initial(sys,100)
opt = stepDataOptions('StepAmplitude',9);
step(sys,opt)

figure(2);
Q = 0.1;
R = 1;
K = lqr(A,B,Q,R);
[numCL,denCL] = ss2tf(A-B*K,B,C,D);
CLtf = tf(numCL,denCL)
F = denCL(end)/numCL(end)
%K = lqr(A,B,Q,R);
sys_lqr = ss((A-B*K),F*B,C,D);
[u,t] = gensig("square",4e-3,6e-3);
% lsim(sys, 9*u, t, 22)
lsim(sys_lqr,90*u,t);

 
 
% t = 0:0.0001:1.4e-3;  % 201 points
% u = max(0,9*min(t,1));

% [u,t] = gensig("square",4e-3,6e-3);
% lsim(sys, 9*u, t, 22)
% 
C1 = pidtune(sys,'PI');
sys1 = feedback(sys*C1,1);

figure(3);
[u1,t1] = gensig("square",4e-3,6e-3);
lsim(sys1, 1e-6*u1, t1)
