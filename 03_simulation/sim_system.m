function y_f = sim_system(sys,u,y,e)

%% Simulation
theta = [sys.B(2:end) -sys.A(2:end)];
y_f = theta*[u; y] + e;

end