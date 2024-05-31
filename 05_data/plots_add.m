function plots_add(data,options)
%% Additional Plots generated with Stored Data Files

%% Appearance & Initialization
layout_options;

MC = options.MC;
T = options.T;
interval = options.interval;
colors = layout_options;

%% Plot 1
f3 = figure;
f3.Position = [200 100 1000 500];
t3 = tiledlayout(2,4);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';

ax1 = nexttile([1 2]);
plot(1:T-1,mean(data3{1,1}.cost_400_wnoise_det(1:end-1,1),2), 1:T-1,mean(data3{1,2}.cost_400_prbs_det(1:end-1,1),2),...
    1:T-1,mean(data3{1,3}.cost_400_act_det(1:end-1,1),2), 'LineWidth',1.5)
title('D-optimality', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$J[I(\theta)]$', 'Interpreter','latex')
legend('$u_{wnoise}$','$u_{prbs}$','$u_{act}$',...
    'Location','northwest', 'Interpreter','latex')

nexttile([1,2])
plot(1:T-1,mean(data4{1,1}.cost_400_wnoise_eig(1:end-1,1),2), 1:T-1,mean(data4{1,2}.cost_400_prbs_eig(1:end-1,1),2),...
    1:T-1,mean(data4{1,3}.cost_400_act_eig(1:end-1,1),2), 'LineWidth',1.5)
title('E-optimality', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$J[I(\theta)]$', 'Interpreter','latex')
legend('$u_{wnoise}$','$u_{prbs}$','$u_{act}$',...
    'Location','northwest', 'Interpreter','latex')

nexttile(6,[1,2])
plot(1:T-1,mean(data5{1,1}.cost_400_wnoise_tr(1:end-1,1),2), 1:T-1,mean(data5{1,2}.cost_400_prbs_tr(1:end-1,1),2),...
    1:T-1,mean(data5{1,3}.cost_400_act_tr(1:end-1,1),2), 'LineWidth',1.5)
title('A-optimality', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$J[I(\theta)]$', 'Interpreter','latex')
legend('$u_{wnoise}$','$u_{prbs}$','$u_{act}$',...
    'Location','northwest', 'Interpreter','latex')

end