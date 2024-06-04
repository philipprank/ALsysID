function plots_main(options,sys_id,sys_sim,data_id)

%% Appearance & Initialization
layout_options;

MC = options.MC;
T = options.T;
interval = options.interval;
method = options.method.rand;
layout = layout_options;
colors = layout.colors;

theta_sim = [sys_sim.B(2:end) sys_sim.A(2:end)];
theta_rand = sys_id.rand;
theta_rcdhz = sys_id.rcdhz;

avg_rand = mean(theta_rand,3);
avg_rcdhz = mean(theta_rcdhz,3);


%% Plot 1
f1 = figure;
f1.Position = [200 100 1000 500];

t1 = tiledlayout(3,1);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

ax1 = nexttile;
set(ax1,'xticklabel',[])
plot(1:MC,theta_sim(2)*ones(MC,1), 1:MC,squeeze(theta_rand(end,3,:)), 1:MC,squeeze(theta_rcdhz(end,3,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{1} = ' num2str(theta_sim(2)) ' \quad \hat{\theta}_{1,method} = ' num2str(round(mean(theta_rand(end,3,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(end,3,:)),3)) ') ' ...
    '\quad \hat{\theta}_{1,act} = ' num2str(round(mean(theta_rcdhz(end,3,:)),3)) ' \; (\sigma = ' num2str(round(std(theta_rcdhz(end,3,:)),3)) ')$'])
ylabel('$\theta_{1}$')
legend('$\theta$','$\hat{\theta}_{prbs}$','$\hat{\theta}_{act}$',...
    'Location','southeast')

ax2 = nexttile;
set(ax2,'xticklabel',[])
plot(1:MC,theta_sim(3)*ones(MC,1), 1:MC,squeeze(theta_rand(end,5,:)), 1:MC,squeeze(theta_rcdhz(end,5,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{2} = ' num2str(theta_sim(3)) ' \quad \hat{\theta}_{2,prbs} = ' num2str(round(mean(theta_rand(end,5,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(end,5,:)),3)) ') ' ...
    '\quad \hat{\theta}_{2,act} = ' num2str(round(mean(theta_rcdhz(end,5,:)),3)) ' \; (\sigma = ' num2str(round(std(theta_rcdhz(end,5,:)),3)) ')$'])
ylabel('$\theta_{2}$')

ax3 = nexttile;
plot(1:MC,theta_sim(4)*ones(MC,1), 1:MC,squeeze(theta_rand(end,6,:)), 1:MC,squeeze(theta_rcdhz(end,6,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{3} = ' num2str(theta_sim(4)) ' \quad \hat{\theta}_{3,prbs} = ' num2str(round(mean(theta_rand(end,6,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(end,6,:)),3)) ') ' ...
    '\quad \hat{\theta}_{3,act} = ' num2str(round(mean(theta_rcdhz(end,6,:)),3)) ' \; (\sigma = ' num2str(round(std(theta_rcdhz(end,6,:)),3)) ')$'])
xlabel('Simulations')
ylabel('$\theta_{3}$')

%% Plot 2
f2 = figure;
f2.Position = [200 100 1000 500];

t2 = tiledlayout(2,3);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';

ax1 = nexttile([1 3]);
plot(50:interval:T,avg_rand(:,3), 'Color', colors(1,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rcdhz(:,3), '--', 'Color', colors(2,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rand(:,5), 'Color', colors(1,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rcdhz(:,5), '--', 'Color', colors(2,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rand(:,6), 'Color', colors(1,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rcdhz(:,6), '--', 'Color', colors(2,:), 'LineWidth', 1.5)

xlim([50 T])
title(['Averages (' num2str(MC) '$\ $Monte-Carlo Simulations)'])
xlabel('Observations T')
ylabel('$\mathrm{avg}({\theta})$')
legend('$u_{prbs}$','$u_{act}$', 'Location','best')


ax2 = nexttile;
set(ax2,'xticklabel',[])
plot(50:interval:T,std(squeeze(theta_rand(:,3,:)),[],2), 50:interval:T,std(squeeze(theta_rcdhz(:,3,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{1})$', 'Interpreter','latex')
legend('$u_{prbs}$','$u_{act}$', 'Location','northeast', 'Interpreter','latex')

ax3 = nexttile;
set(ax3,'xticklabel',[])
plot(50:interval:T,std(squeeze(theta_rand(:,5,:)),[],2), 50:interval:T,std(squeeze(theta_rcdhz(:,5,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{2})$', 'Interpreter','latex')
legend('$u_{prbs}$','$u_{act}$', 'Location','northeast', 'Interpreter','latex')

ax4 = nexttile;
plot(50:interval:T,std(squeeze(theta_rand(:,6,:)),[],2), 50:interval:T,std(squeeze(theta_rcdhz(:,6,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{3})$', 'Interpreter','latex')
legend('$u_{prbs}$','$u_{act}$',...
    'Location','northeast', 'Interpreter','latex')

%% Plot 3
f3 = figure;
f3.Position = [200 100 1000 500];

t3 = tiledlayout(2,2);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';

ax1 = nexttile;
surf(abs((ones(36,MC)*sys_sim.B(3) - squeeze(theta_rcdhz(:,3,:)))./(ones(36,MC)*sys_sim.B(3))*100) - ...
    abs((ones(36,MC)*sys_sim.B(3) - squeeze(theta_rand(:,3,:)))./(ones(36,MC)*sys_sim.B(3))*100))
view(2)
map = [colors(3,:); colors(1,:); colors(4:5,:)];
colormap(map)
colorbar
clim([-60 20])

ax2 = nexttile;
surf(abs((ones(36,MC)*sys_sim.A(2) - squeeze(theta_rcdhz(:,5,:)))./(ones(36,MC)*sys_sim.A(2))*100) - ...
    abs((ones(36,MC)*sys_sim.A(2) - squeeze(theta_rand(:,5,:)))./(ones(36,MC)*sys_sim.A(2))*100))
view(2)
map = [colors(3,:); colors(1,:); colors(4:5,:)];
colormap(map)
colorbar
clim([-60 20])

ax3 = nexttile;
surf(abs((ones(36,MC)*sys_sim.A(3) - squeeze(theta_rcdhz(:,6,:)))./(ones(36,MC)*sys_sim.A(3))*100) - ...
    abs((ones(36,MC)*sys_sim.A(3) - squeeze(theta_rand(:,6,:)))./(ones(36,MC)*sys_sim.A(3))*100))
view(2)
map = [colors(3,:); colors(1,:); colors(4:5,:)];
colormap(map)
colorbar
clim([-60 20])

%% Plot 4
f4 = figure;
f4.Position = [200 100 1000 500];
t4 = tiledlayout(3,1);
t4.TileSpacing = 'compact';
t4.Padding = 'compact';
title(t4,'Input \& Output Signals for Second-Order System', 'Interpreter','latex')

% nexttile
% plot(1:T,data_id_rand(1).u1, 1:T,data_id_rand(1).y1, 'LineWidth',1.5)
% title('White Noise Input Signal', 'Interpreter','latex')
% xlabel('Observations T', 'Interpreter','latex')
% ylabel('Magnitude', 'Interpreter','latex')
% legend('input','output','Location','northeast', 'Interpreter','latex')
% ylim([-5 5])

nexttile
plot(1:T,data_id(1).rand.InputData, 1:T,data_id(1).rand.OutputData, 'LineWidth',1.5)
title('PRBS Input Signal', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('Magnitude', 'Interpreter','latex')
legend('input','output','Location','northeast', 'Interpreter','latex')
ylim([-5 5])


nexttile
plot(1:T,data_id(1).rcdhz.InputData, 1:T,data_id(1).rcdhz.OutputData, 'LineWidth',1.5)
title('D-Optimal designed Input Signal ($H=10$)', 'Interpreter','latex')
xlabel('Observations T', 'Interpreter','latex')
ylabel('Magnitude', 'Interpreter','latex')
legend('input','output','Location','northeast', 'Interpreter','latex')
ylim([-5 5])


%% Plot 5