tic
%% Initializaiton
theta_sim = [0.5, -1.5, 0.7];
var_noise = 0.05;
var_input = 0.25;
MC = 100; % For correlation matrix set MC=1
T = 500;
interval = 50;
t_eval = interval:interval:T;
e = sqrt(var_noise)*randn(T,MC); % white measurement noise (used for simulation)
y = zeros(T,MC);
y_aux = zeros(T,MC);
u = zeros(T,MC);
% Given systen structure
init_sys = idtf([0 NaN 0],[1 NaN(1,2)],1);
init_sys.Structure.Numerator.Free = [0 1 0];
I = zeros(length(theta_sim),length(theta_sim),MC);
cost = zeros(MC,1);

data_id = cell(MC,1);
sys_id = zeros(floor(T/interval-50/interval)+1,6,MC);

%% Input Signal Generation
% White Noise Input Signal
% u = sqrt(var_input)*randn(T,MC);

% Designed Feedback Input Signal
psi = sqrt(1/theta_sim(1)^2)*randn(T,MC); % random part
u(1:2,:) = zeros(2,MC);

%% System Simulation & Parameter Identificatio
for i = 1:MC
    idx = 1;
    for j = 3:T-1
        y(j,i) = theta_sim*[u(j-1,i); -y(j-1:-1:j-2,i)] + e(j,i);
        % y_aux(j,i) = theta_sim*[u(j-1,i); -y_aux(j-1:-1:j-2,i)];
        u(j,i) = 1/theta_sim(1)*(theta_sim(2)*y(j,i) + theta_sim(3)*y(j-1,i)) + psi(j,i);
        if any(t_eval(:) == j) && j >= 50
            data_id{i,1} = iddata(y(1:j,i,1),u(1:j,i,1),1);
            id_struct = oe(data_id{i,1},init_sys);
            sys_id(idx,:,i) = [id_struct.B id_struct.F];
            idx = idx + 1;
        end
    end
    data_id{i,1} = iddata(y(1:j,i,1),u(1:j,i,1),1);
    id_struct = oe(data_id{i,1},init_sys);
    sys_id(idx,:,i) = [id_struct.B id_struct.F];
    [I(:,:,i),cost(i)] = fisher_information(u(:,i),y(:,i));
end

%% Plots
layout_options
layout = layout_options;
colors = layout.colors;
theta_rand = sys_id;

% Figure 1
f1 = figure;
f1.Units = 'centimeters';
f1.Position = [8 4 15 10];
t1 = tiledlayout(3,1);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

ax1 = nexttile;
set(ax1,'xticklabel',[])
plot(1:MC,theta_sim(1)*ones(MC,1), 1:MC,squeeze(theta_rand(end,2,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{1} = ' num2str(theta_sim(1)) ' \quad \hat{\theta}_{1} = ' num2str(round(mean(theta_rand(end,2,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(1,2,:)),3)) ')$'])
ylabel('$\theta_{1}$')
legend('$\theta$','$\hat{\theta}_{opt}$', 'Location','southeast', 'FontSize', 8.5)

ax2 = nexttile;
set(ax2,'xticklabel',[])
plot(1:MC,theta_sim(2)*ones(MC,1), 1:MC,squeeze(theta_rand(end,5,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{2} = ' num2str(theta_sim(2)) ' \quad \hat{\theta}_{2} = ' num2str(round(mean(theta_rand(end,5,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(1,5,:)),3)) ')$'])
ylabel('$\theta_{2}$')

ax3 = nexttile;
plot(1:MC,theta_sim(3)*ones(MC,1), 1:MC,squeeze(theta_rand(end,6,:)), 'LineWidth', 1.5)
xlim([1 MC])
title(['$ \theta_{3} = ' num2str(theta_sim(3)) ' \quad \hat{\theta}_{3} = ' num2str(round(mean(theta_rand(end,6,:)),3)) '\; (\sigma = ' num2str(round(std(theta_rand(1,6,:)),3)) ')$'])
xlabel('Simulations M')
ylabel('$\theta_{3}$')

% Figure 2
f2 = figure;
f2.Units = 'centimeters';
f2.Position = [8 4 15 10];
t2 = tiledlayout(2,3);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';

avg_rand = mean(theta_rand,3);

ax1 = nexttile([1 3]);
plot(50:interval:T,avg_rand(:,2), 'Color', colors(1,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rand(:,5), 'Color', colors(1,:), 'LineWidth', 1.5)
hold all
plot(50:interval:T,avg_rand(:,6), 'Color', colors(1,:), 'LineWidth', 1.5)

xlim([50 T])
ylim([-1.6 1.2])
title(['Averages (' num2str(MC) '$\ $Monte-Carlo Simulations)'])
xlabel('Observations N')
ylabel('$\mathrm{avg}({\theta})$')
legend('$u_{prbs}$','$u_{act}$', 'Location','best', 'FontSize', 8.5)


ax2 = nexttile;
set(ax2,'xticklabel',[])
plot(50:interval:T,std(squeeze(theta_rand(:,2,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations N', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{1})$', 'Interpreter','latex')

ax3 = nexttile;
set(ax3,'xticklabel',[])
plot(50:interval:T,std(squeeze(theta_rand(:,5,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations N', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{2})$', 'Interpreter','latex')

ax4 = nexttile;
plot(50:interval:T,std(squeeze(theta_rand(:,6,:)),[],2), 'LineWidth', 1.5)
xlim([50 T])
title('Standard Deviation', 'Interpreter','latex')
xlabel('Observations N', 'Interpreter','latex')
ylabel('$\sigma({\theta}_{3})$', 'Interpreter','latex')
legend('$u_{prbs}$','$u_{act}$', 'Location','northeast', 'Interpreter','latex', 'FontSize', 8.5)


%% Correlation Matrix
% Extrac relevant data for exemplary system
y_current = y(3:end,5); % y(t)
u_current = u(3:end,5); % u(t)
y_prev = y(2:end-1,5); % y(t-1)
% u_prev = u(2:end-1, :); % u(t-1)

% Data Matrix
data_matrix = [y_current(:), u_current(:), y_prev(:),];

% Correlation
corr_matrix = corrcoef(data_matrix);

% Correlation Matrix as a Heatmap
figure;
imagesc(corr_matrix);
colorbar;
axis equal tight;
title('Correlation Matrix');
set(gca, 'XTick', 1:4, 'XTickLabel', {'y(t)', 'u(t)', 'y(t-1)', 'u(t-1)'});
set(gca, 'YTick', 1:4, 'YTickLabel', {'y(t)', 'u(t)', 'y(t-1)', 'u(t-1)'});

for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
            'HorizontalAlignment', 'center', 'Color', 'black', 'FontSize', 12);
    end
end

%% Fisher information
function [I, cost] = fisher_information(u,y)
N = length(y);
I = zeros(3,3);
for t = 3:N
    dL_dtheta1 = u(t-1);
    dL_dtheta2 = -y(t-1);
    dL_dtheta3 = -y(t-2);
    grad_L = [dL_dtheta1; dL_dtheta2; dL_dtheta3];
    I = I + (grad_L * grad_L');
    cost = log(det(I));
end
end
toc