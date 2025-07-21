function y_next = sim_system(sys,u,y,e)

%% Simulation
% (MIMO) ARX Model
% y_{i} = -a_{1}*y_{i-1}-...-y_{p}*y_{i-p}+b_{1}*u_{i-1}+...+b_{q}*u_{i-q}
% + noise
theta = [sys.B(2:end) -sys.A(2:end)];
y_next = theta*[u; y] + e;

end