function [data_id,sys_id] = opt_input_freq(sys,options,e)
%%==================================%%
%%%     ACTIVE INPUT DESIGN        %%%
%%%     (FREUENCY DOMAIN)          %%%

% Optimal input design by optimizing offline
% some cost function of the Fisher information
% in the frequency domain, resluting in
% periodic input/output signals
%%==================================%%

%% Initialization
MC = options.MC;
T = options.T;
interval = options.interval;


end