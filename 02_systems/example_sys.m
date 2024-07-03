function sys = example_sys(options)

if strcmp(options.sys,'sys1') == true
%% Example 1 (Second-Order-System)
%%====================%%
% Transfer Function:
% G(q) = (theta_1*q^(-2))/(1+theta_2*q^(-1)+theta_3*q^(-2))
%%====================%%
na = options.param.sys1.na;
nb = options.param.sys1.nb;
sys = idpoly(na,nb,1,1,1,[],options.Ts);

elseif strcmp(options.sys,'sys2') == true
%% Example 2 (Flexible Transmission System)
%%====================%%
% Taken from: "A Flexible Transmission System as a Benchmark for Robust
% Digital Control"; Landau 1995
% https://doi.org/10.1016/S0947-3580(95)70011-5
% Transfer Function (unloaded case):
% G(q) = (theta_1*q^(-3) + theta_2*q^(-4))/
%        (1 - theta_3*q^(-1) - theta_4*q^(-2) - theta_5*q^(-3) - theta_6*q^(-4))
%%====================%%
na = options.param.sys2.na;
nb = options.param.sys2.nb;  
sys = idpoly(na,nb,1,1,1,[],options.Ts);

elseif ex == 3
%% Example 3 (Coupled Water Tanks)
%%====================%%
% YET TO DO
% Transfer Function (MIMO) [2x2]:
% G11(q) =
% G12(q) =
% G22(q) =
%%====================%%
end
