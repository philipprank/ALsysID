clear clc;
h1max = 1.36; % m
h2max = 1.36; % m
h3max = 1.30; % m
h4max = 1.30; % m
hmin = 0.20; % m (same for all tanks)

qamax = 3.26/3600;
qbmax = 4/3600;
qmin = 0/3600;

a1 = 1.31e-4; % m2
a2 = 1.51e-4; % m2
a3 = 9.27e-5; % m2
a4 = 8.82e-5; % m2

S = 0.06; % m2

ga = 0.30; % always between 0 and 1
gb = 0.40; % always between 0 and 1

h1o = 0.65; % m
h2o = 0.66; % m
h3o = 0.65; % m
h4o = 0.66; % m

qao = 1.63/3600;
qbo = 2/3600;

g = 9.81; % m/s2

Ts = 25; % s (Sampling Time)

T1 = S*sqrt(2*h1o/g)/a1;% 1st Time Constant
T2 = S*sqrt(2*h2o/g)/a2;% 2nd Time Constant
T3 = S*sqrt(2*h3o/g)/a3;% 3rd Time Constant
T4 = S*sqrt(2*h4o/g)/a4;% 4th Time Constant