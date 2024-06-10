%%================================%%
%%%%   Thesis Plots             %%%%
%%================================%%

export = 1; % 0 or 1 (for exporting plots to PDF)

%% Figure 3.1 - Probability Density Functions
layout_options;

f1 = figure;
f1.Units = 'centimeters';
f1.Position = [8 4 14 8];

mu_range = -2:0.05:10;
mu = 4;
sigma1 = 1;
sigma2 = 4;
obsv = 20;

y1 = normpdf(mu_range,mu,sigma1);
y2 = normpdf(mu_range,mu,sigma2);
plot(mu_range,y1,mu_range,y2, 'LineWidth',1.5)
xlabel('X')
ylabel('Probability Density p')
ylim([0 0.42])
legend(sprintf('$\\sigma_{1}^2 = %i$', sigma1),sprintf('$\\sigma_{2}^2 = %i$', sigma2),...
    'Location','northeast')


%% Figure 3.2 - Likelihood Functions & Score Functions
data1 = normrnd(mu,sigma1,obsv,1);
data2 = normrnd(mu,sigma2,obsv,1);

log_likelihood1 = zeros(length(mu_range),obsv);
log_likelihood2 = zeros(length(mu_range),obsv);
score1 = zeros(length(mu_range)-1,obsv);
score2 = zeros(length(mu_range)-1,obsv);

for j = 1:obsv
    for i = 1:length(mu_range)
        mu = mu_range(i);
        log_likelihood1(i,j) = -0.5*length(data1)*log(2*pi) - length(data1)*log(sigma1)...
            - 0.5*sum(((data1(j)-mu)/sigma1).^2);
        log_likelihood2(i,j) = -0.5*length(data2)*log(2*pi) - length(data2)*log(sigma2)...
            - 0.5*sum(((data2(j)-mu)/sigma2).^2);
        score1(:,j) = diff(log_likelihood1(:,j))./diff(mu_range)';
        score2(:,j) = diff(log_likelihood2(:,j))./diff(mu_range)';
    end

end
mu_d = (mu_range(2:end)+mu_range(1:(end-1)))/2;

pd1 = fitdist(score1(121,:)','Normal');
pd2 = fitdist(score2(121,:)','Normal');

mu_score_range = -3:0.01:3;
mu_score = 0; % corrected! (true values not exactly zero)
sigma_score1 = pd1.sigma;
sigma_score2 = pd2.sigma;

y_score1 = normpdf(mu_score_range,mu_score,sigma_score1);
y_score2 = normpdf(mu_score_range,mu_score,sigma_score2);

f2 = figure;
f2.Units = 'centimeters';
f2.Position = [8 4 14 12];
t1 = tiledlayout(3,2);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

nexttile
hold on;
for j = 1:obsv
    plot(mu_range,log_likelihood1(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlim([mu_range(1) mu_range(end)])
ylim([-55 -15])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood $\mathcal{L}$');
yticks([-55 -35 -15])
% legend(sprintf('$\\sigma_{1} = %i$', sigma1))

nexttile
hold on
for j = 1:obsv
    plot(mu_range,log_likelihood2(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off;
xlim([mu_range(1) mu_range(end)])
ylim([-55 -45])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood $\mathcal{L}$');
yticks([-55 -50 -45])
% legend(sprintf('$\\sigma_{2} = %i$', sigma2))

nexttile
hold on
for j = 1:obsv
    plot(mu_d,score1(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off;
xlabel('Parameter $\mu$');
ylabel('Score $\partial \mathcal{L} / \partial \theta$')
xlim([0 8])
ylim([-4 4])
% yticks([-10 0 10])

nexttile
hold on
for j = 1:obsv
    plot(mu_d,score2(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off;
xlabel('Parameter $\mu$');
ylabel('Score $\partial \mathcal{L} / \partial \theta$')
xlim([0 8])
ylim([-4 4])
% yticks([-1 0 1])

nexttile
plot(mu_score_range,y_score1, 'LineWidth',1.5)
xlabel('Data')
ylabel('pdf')
xlim([-2.5 2.5])
ylim([0 2])

nexttile
plot(mu_score_range,y_score2, 'LineWidth',1.5)
xlabel('Data')
ylabel('pdf')
xlim([-2.5 2.5])
ylim([0 2])


%% Figure 3.3 - Confidence Ellipses
% Ellipse parameters
a = 5; % axis 1
b = 3; % axis 2
h = 2; % x-coordinate (center)
k = 1; % y-coordinate (center)
theta = pi/4; % rotation [rad]

% Ellipse
t = linspace(0,2*pi,100);
x = h + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
y = k + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);

% Bounding Box (A-optimality)
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

rect_x = [x_min x_max x_max x_min x_min];
rect_y = [y_min y_min y_max y_max y_min];

% (Scaled) Axis (D-/E-optimality)
u1 = [a*cos(theta); a*sin(theta)];
u2 = [-b*sin(theta); b*cos(theta)];

major_x = [h, h+u1(1)]; % endpoints
major_y = [k, k+u1(2)];
minor_length = sqrt((x_max-x_min)^2 + (y_max-y_min)^2)/2;
u2_norm = u2/norm(u2)*minor_length;
minor_x = [h, h+u2_norm(1)];
minor_y = [k, k+u2_norm(2)];

% Plot
f3 = figure;
f3.Units = 'centimeters';
f3.Position = [8 4 14 8];
plot(x,y,'LineWidth',1.5,'Color',[0.2431 0.2667 0.2980]);

hold on

fill(x,y, [1.0000 0.8353 0.0000],'FaceAlpha', 0.3)

plot(rect_x,rect_y, '--','LineWidth',1,'Color',[0.0000 0.7451 1.0000])

plot(h,k, 'o','MarkerSize',4,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);
quiver(h,k,u1(1),u1(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.3176 0.6196],'LineWidth',1);
quiver(h,k,u2_norm(1),u2_norm(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.7451 1.0000],'LineWidth',1);

% Labels
text(4.5,-2, 'A-optimality','FontSize',11, 'HorizontalAlignment','center', 'Color',[1.0000 0.8353 0.0000]);
text(-0.4,4.6, 'D-optimality','FontSize',11, 'HorizontalAlignment','center', 'Color',[0.0000 0.7451 1.0000]);
text(5.1,2.8, 'E-optimality','FontSize',11 ,'HorizontalAlignment','center', 'Color',[0.0000 0.3176 0.6196]);

% Format
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xlim([-2.5 6.5])
ylim([-3.7 5.7])
xlabel('$\theta_{1}$')
ylabel('$\theta_{2}$')
box off

axp = get(gca,'Position');
xs=axp(1);
xe=axp(1)+axp(3)+0.04;
ys=axp(2);
ye=axp(2)+axp(4)+0.05;
annotation('arrow', [xs xe],[ys ys], 'LineWidth',1,'HeadStyle','vback3',...
    'HeadWidth',8,'HeadLength',8);
annotation('arrow', [xs xs],[ys ye], 'LineWidth',1,'HeadStyle','vback3',...
    'HeadWidth',8,'HeadLength',8);
hold off


%% Export Graphics
if export == 1
    f1name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_1.pdf');
    exportgraphics(f1,f1name,'ContentType','vector');
    f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_2.pdf');
    exportgraphics(f2,f2name,'ContentType','vector');
    f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_3.pdf');
    exportgraphics(f3,f3name,'ContentType','vector');
end