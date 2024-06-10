function sim_options = sim_options
% Systems:
% 1. Arbitrary oe-model (idpoly object)
%
% 2. Example systems
%  i.   FIR-System (yet need to be implemented)
%  i.   Second Order System (SISO, order: na=2 & nb=2) -> 'sys1'
%  ii.  Flexible Transmission System (SISO, order: na=4 & nb=4) -> 'sys2'
%  iii. Coupled Water Tank (MIMO) -> 'sys3'
%
% Parameters:
% 1. Horizon length H (depends on noise level, good choice for example
%    systems -> H = 6-10)
% 2. Optimality Criteria
%  i.   A-optimal design -> 'a-opt'
%  ii.  D-optimal design -> 'd-opt'
%  iii. E-optimal design -> 'e-opt'
% All designed inputs gave better results for the consdiered examples,
% however the best choice was D-optimal design (under appropriate
% sim_options)
% 3. Number of Monte-Carlo Simuations (MC)
% 4. Observations (T)

%% Load System Parameters
load sys_param.mat sys_param;

sim_options = struct();
sim_options.param = sys_param;

%% General Options
sim_options.MC = 200;
sim_options.T = 400;
sim_options.Ts = 1;
% don't save or save simulated data (0 or 1)
sim_options.svdata = 0;

%% System Options
% 'sys1', 'sys2', 'sys3', 'sys4' or idpoly with A(q),C(q),D(q) = 0
sim_options.sys = 'sys1';
sim_options.var_noise = 0.01;
sim_options.amp_con = 1;

%% Input Options
% 'prbs', 'wnoise', 'chirp' 'all'
sim_options.method.rand = 'prbs';
%'act_rcdhz', 'act_freq', 'act_cvrlx', 'all'
sim_options.method.active = 'act_rdhz';
% 'a-opt', 'd-opt', 'e-opt'
sim_options.opt = 'd-opt';
% 'alg1', 'alg2', 'alg3', 'all'
sim_options.alg = 'alg1';
% 'true', 'recursive'
sim_options.estimation = 'recursive';
sim_options.init_T = 60;
% relevant only for 'true' and thereby created plots
sim_options.interval = 10; % time steps, where parameters get evaluated
% integer or intervall for comparison, e.g. [4 10] or [4 2 10]
% (needs to be sufficiently large)
sim_options.H = 10;


end