%% Set parameter variable values
g = 9.81; % Gravity
cl = 0.4; % Standard cl value
cd = 0.05; % Standard cd value
p_air = 1.225; % in kg/m^3
m = 13.5; % in kg
kn = 500; % in N-m
A_ref = 0.5; % in m^2
speed_c = 9;  % in m/s
u1 = 1;
u2 = 1;
u3 = pi;
u4 = pi/2;
u5 = pi/2;
u6 = pi;

%% Linearized state-space equations about x = [speed_c 0 0 0 0 0 0 0]
%syms u1 u2 u3 u4 u5 u6

A_lin = [-1*cd*p_air*A_ref*speed_c/m 0 0 0 0 0 0 0;
    cl*p_air*A_ref*speed_c/m 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 -1;
    0 0 u6 0 0 0 0 0;
    0 0 0 0 0 0 0 u6;
    0 0 0 0 0 0 0 u4];
B_lin = [-1*kn*sin(u6)/m 1/m -1*kn*u1*cos(u3)/m 0 0 0;
    kn*sin(u6)/m 0 kn*u1*cos(u3)/m 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0;
    0 0 0 1 0 0;
    0 0 0 0 0 -1];

C = eye(8);

sys = ss(A_lin,B_lin,C,0);

%% Using Bryson's rule to set LQR Weights
x1_max = 15; % in m/s
x2_max = 4; % in m/s
x3_max = 10*pi/180; % in rad
x4_max = 10*pi/180; % in rad
x5_max = 10*pi/180; % in rad
x6_max = 10*pi/180; % in rad/s
x7_max = 10*pi/180; % in rad/s
x8_max = 10*pi/180; % in rad/s

u1_max = 5; % in /s
u2_max = 20; % in N
u3_max = 5*pi/180; % in rad
u4_max = 10*pi/180; % in rad
u5_max = 10*pi/180; % in rad
u6_max = 10*pi/180; % in rad


Q = [1/(x1_max)^2 0 0 0 0 0 0;
    0 1/(x2_max)^2 0 0 0 0 0;
    0 0 1/(x4_max)^2 0 0 0 0;
    0 0 0 1/(x5_max)^2 0 0 0;
    0 0 0 0 1/(x6_max)^2 0 0;
    0 0 0 0 0 1/(x7_max)^2 0;
    0 0 0 0 0 0 1/(x8_max)^2];

R = [1/(u1_max)^2 0 0 0 0 0;
   0 1/(u2_max)^2 0 0 0 0;
   0 0 1/(u3_max)^2 0 0 0;
   0 0 0 1/(u4_max)^2 0 0;
   0 0 0 0 1/(u5_max)^2 0;
   0 0 0 0 0 1/(u6_max)^2];

%R = eye(6);

%% Solve LQR
sysr = minreal(sys)
[K, S, e] = lqr(sysr, Q,R)
