 %%============================================%%
%%%%   Optimal Input Design Algorithm for   %%%% 
%%%%      Partially Observable Systems      %%%%
%%============================================%%

tic
%% Load Options
options = sim_options;

%% Generate True System
sys_sim = example_sys(options);

%% Simulation
e = sqrt(options.var_noise)*randn(options.T,options.MC); % white noise (used for simulation)

[sys_id_rcdhz,data_id_rcdhz] = opt_input_rcdhz(sys_sim,options,e);
[sys_id_rand,data_id_rand] = input_rand(sys_sim,options,e);

sys_id = struct('rand',sys_id_rand,'rcdhz',sys_id_rcdhz);
data_id = struct('rand',data_id_rand,'rcdhz',data_id_rcdhz);

%% Data % Plots (Selection)
id_struct = save_data(options,sys_id,data_id); 
plots_main(options,sys_id,sys_sim,data_id)

toc