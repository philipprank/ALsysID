%%================================%%
%%%%   Thesis Plots             %%%%
%%================================%%
clearvars
export = 1; % 0 or 1 (for exporting plots to PDF)

obsv = 1000;
eval = floor(obsv/20);

%% Figure 3.1 - Probability Density Functions
layout_options;

f1 = figure;
f1.Units = 'centimeters';
f1.Position = [8 4 14 8];

mu_range = -2:0.05:10;
mu = 4;
sigma1 = 1;
sigma2 = 4;

y1 = normpdf(mu_range,mu,sigma1);
y2 = normpdf(mu_range,mu,sigma2);
plot(mu_range,y1,mu_range,y2, 'LineWidth',1.5)
xlabel('X')
ylabel('Probability Density p')
ylim([0 0.42])
legend(sprintf('$\\sigma_{1}^2 = %i$', sigma1),sprintf('$\\sigma_{2}^2 = %i$', sigma2),...
    'Location','northeast')


%% Fi - Likelihood Functions & Score Functions
data1 = normrnd(mu,sigma1,obsv,1);
data2 = normrnd(mu,sigma2,obsv,1);

log_likelihood1 = zeros(length(mu_range),obsv);
log_likelihood2 = zeros(length(mu_range),obsv);
score1 = zeros(length(mu_range)-1,obsv);
score2 = zeros(length(mu_range)-1,obsv);

for j = 1:obsv
    for i = 1:length(mu_range)
        mu = mu_range(i);
        log_likelihood1(i,j) = -0.5*log(2*pi) - 0.5*log(sigma1) - 0.5*sum(((data1(j)-mu)/sigma1).^2);
        log_likelihood2(i,j) = -0.5*log(2*pi) - 0.5*log(sigma2) - 0.5*sum(((data2(j)-mu)/sigma2).^2);
        score1(:,j) = diff(log_likelihood1(:,j))./diff(mu_range)';
        score2(:,j) = diff(log_likelihood2(:,j))./diff(mu_range)';
    end
end
mu_d = (mu_range(2:end) + mu_range(1:(end-1)))/2;

pd1 = fitdist(score1(121,:)','Normal');
pd2 = fitdist(score2(121,:)','Normal');

mu_score_range = -3:0.01:3;
mu_score = 0; % corrected! (true values may not be exactly zero)
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
hold on
for j = 1:eval:obsv
    plot(mu_range,log_likelihood1(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
% plot(mu_range,sum_log_likelihood1, 'LineWidth',1.5, 'Color',[1.0000 0.8353 0.0000]);
hold off
xlim([mu_range(1) mu_range(end)])
ylim([-15 3])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood $f$');
yticks([-15 -10 -5 0])
% legend(sprintf('$\\sigma_{1} = %i$', sigma1))

nexttile
hold on
for j = 1:eval:obsv
    plot(mu_range,log_likelihood2(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlim([mu_range(1) mu_range(end)])
ylim([-15 3])
xlabel('Parameter $\mu$');
ylabel('log-Likelihood f');
yticks([-15 -10 -5 0])
% legend(sprintf('$\\sigma_{2} = %i$', sigma2))

nexttile
hold on
for j = 1:eval:obsv
    plot(mu_d,score1(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlabel('Parameter $\mu$');
ylabel('Score $\partial f / \partial \theta$')
xlim([0 8])
ylim([-4 4])
% yticks([-10 0 10])

nexttile
hold on
for j = 1:eval:obsv
    plot(mu_d,score2(:,j), 'Color',[0.0000 0.3176 0.6196]);
end
hold off;
xlabel('Parameter $\mu$');
ylabel('Score $\partial f / \partial \theta$')
xlim([0 8])
ylim([-4 4])
% yticks([-1 0 1])

nexttile
plot(mu_score_range,y_score1, 'LineWidth',1.5)
xlabel('$\partial f / \partial \theta$ $(\mu = 4)$')
ylabel('pdf')
xlim([-2.5 2.5])
ylim([0 2])

nexttile
plot(mu_score_range,y_score2, 'LineWidth',1.5)
xlabel('$\partial f / \partial \theta$ $(\mu = 4)$')
ylabel('pdf')
xlim([-2.5 2.5])
ylim([0 2])

%% Figure 3.3 - Probability Density Functions (multivariate)
f3 = figure;
f3.Units = 'centimeters';
f3.Position = [8 4 14 8];
t1 = tiledlayout(2,2);
t1.TileSpacing = 'compact';
t1.Padding = 'compact';

mu = [4 4];
rho1 = -0.2;
rho2 = -0.8;
% sigma1_multi = [1 -0.5; -0.5 1];
sigma1_multi = [0.6 rho1; rho1 0.6];
% sigma2_multi = [3 -1.5; -1.5 3];
sigma2_multi = [1 rho2; rho2 1];


x = 0:0.15:8;
y = 0:0.15:8;
[X,Y] = meshgrid(x,y);

F1 = mvnpdf([X(:) Y(:)],mu,sigma1_multi);
F1 = reshape(F1,length(x),length(x));
F2 = mvnpdf([X(:) Y(:)],mu,sigma2_multi);
F2 = reshape(F2,length(x),length(x));

nexttile
surf(X,Y,F1);
xlabel('$X_{1}$');
ylabel('$X_{2}$');
zlabel('Joint pdf p');
xticks([0 2 4 6])
yticks([0 2 4 6])
zlim([0 0.3])

nexttile
surf(X,Y,F2);
xlabel('$X_{1}$');
ylabel('$X_{2}$');
zlabel('Joint pdf p');
xticks([0 2 4 6])
yticks([0 2 4 6])
zlim([0 0.3])

nexttile
contourf(X,Y,F1)
xlabel('$X_{1}$');
ylabel('$X_{2}$');

nexttile
contourf(X,Y,F2)
xlabel('$X_{1}$');
ylabel('$X_{2}$');

%% Figure 3.4 - Likelihood Functions & Score Functions (multivariate)
mu1_multi = [4 4];
mu_multi_range = -2:0.05:10;
% sigma1_multi = [1 -0.5; -0.5 1];
% sigma2_multi = [1 -0.85; -0.85 1];
data1_multi = mvnrnd(mu1_multi,sigma1_multi,obsv);
data2_multi = mvnrnd(mu1_multi,sigma1_multi,obsv);

log_likelihood1_multi = zeros(length(mu_multi_range),length(mu_multi_range),obsv);
log_likelihood2_multi = zeros(length(mu_multi_range),length(mu_multi_range),obsv);
score1_multi = zeros(length(mu_multi_range),length(mu_multi_range),obsv);
score2_multi = zeros(length(mu_multi_range),length(mu_multi_range),obsv);

for k = 1:obsv
    for i = 1:length(mu_multi_range)
        for j = 1:length(mu_multi_range)
            mu_multi = [mu_multi_range(i) mu_multi_range(j)];
            log_likelihood1_multi(i,j,k) = -log(2*pi) - 0.5*log(det(sigma1_multi))...
                - 0.5*sum((data1_multi(k,:) - mu_multi)*inv(sigma1_multi)*(data1_multi(k,:) - mu_multi)');
            log_likelihood2_multi(i,j,k) = -log(2*pi) - 0.5*log(det(sigma2_multi))...
                - 0.5*sum((data2_multi(k,:) - mu_multi)*inv(sigma2_multi)*(data2_multi(k,:) - mu_multi)');
        end
    end
end
[Gx1,Gy1] = gradient(log_likelihood1_multi,0.05);
Gx1 = squeeze(Gx1(121,121,:));
Gy1 = squeeze(Gy1(121,121,:));
[Gy2,Gx2] = gradient(log_likelihood2_multi,0.05);
Gx2 = squeeze(Gx2(121,121,:));
Gy2 = squeeze(Gy2(121,121,:));

[counts1,bins1] = hist3([Gx1 Gy1], 'Edges',...
    {linspace(1.2*min(Gx1),1.2*max(Gx1),18) linspace(1.2*min(Gy1),1.2*max(Gy1),18)});
[counts2,bins2] = hist3([Gx2 Gy2], 'Edges',...
    {linspace(1.2*min(Gx2),1.2*max(Gx2),18) linspace(1.2*min(Gy2),1.2*max(Gy2),18)});

% mu,sigma needs more observations for exact results (e.g. obsv>1000)
mu_score1_multi = mean([Gx1 Gy1]);
% sigma_score1_multi = cov([Gx1 Gy1]);
mu_score2_multi = mean([Gx2 Gy2]);
% sigma_score2_multi = cov([Gx2 Gy2]);
% use corrected mu & sigma for plots!
mu_score_multi = [0 0];
sigma_score1_multi = inv(sigma1_multi);
sigma_score2_multi = inv(sigma2_multi);

% Covariance Ellipse
[R1, eigval1] = eig(sigma_score1_multi);
[R2, eigval2] = eig(sigma_score2_multi);
theta = linspace(0,2*pi,100);
ellipse_x_r1 = sqrt(eigval1(1,1))*cos(theta);
ellipse_y_r1 = sqrt(eigval1(2,2))*sin(theta);
ellipse_x_r2 = sqrt(eigval2(1,1))*cos(theta);
ellipse_y_r2 = sqrt(eigval2(2,2))*sin(theta);
ellipse1 = [ellipse_x_r1; ellipse_y_r1]'*R1';
ellipse2 = [ellipse_x_r2; ellipse_y_r2]'*R2';

f4 = figure;
f4.Units = 'centimeters';
f4.Position = [8 4 14 12];
t2 = tiledlayout(3,2);
t2.TileSpacing = 'compact';
t2.Padding = 'compact';

nexttile
hold on
for j = 1:eval:obsv
    contour(mu_multi_range,mu_multi_range,log_likelihood1_multi(:,:,j),[-4 -4],...
        'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlabel('Parameter $\mu_{1}$');
ylabel('Parameter $\mu_{2}$');

nexttile
hold on
for j = 1:eval:obsv
    contour(mu_multi_range,mu_multi_range,log_likelihood2_multi(:,:,j),[-4 -4],...
        'Color',[0.0000 0.3176 0.6196]);
end
hold off
xlabel('Parameter $\mu_{1}$');
ylabel('Parameter $\mu_{2}$');

nexttile
h1 = imagesc(linspace(-2.3,2.3,18),linspace(-2.3,2.3,16),counts1);
xlabel('$\partial f / \partial \theta_{1}$ $(\mu_{1} = 4)$')
ylabel('$\partial f / \partial \theta_{2}$ $(\mu_{2} = 4)$')
set(gca,'YDir','normal')

nexttile
h2 = imagesc(linspace(-2.3,2.3,18),linspace(-2.3,2.3,16),counts2);
xlabel('$\partial f / \partial \theta_{1}$ $(\mu_{1} = 4)$')
ylabel('$\partial f / \partial \theta_{2}$ $(\mu_{2} = 4)$')
set(gca,'YDir','normal')

nexttile
plot(ellipse1(:,1), ellipse1(:,2), 'Color',[0.0000 0.3176 0.6196], 'LineWidth', 1.5);
xlabel('$\partial f / \partial \theta_{1}$ $(\mu_{1} = 4)$')
ylabel('$\partial f / \partial \theta_{2}$ $(\mu_{2} = 4)$')
xlim([-2.3 2.3])
ylim([-2.3 2.3])

nexttile
plot(ellipse2(:,1), ellipse2(:,2), 'Color',[0.0000 0.3176 0.6196], 'LineWidth', 1.5);
xlabel('$\partial f / \partial \theta_{1}$ $(\mu_{1} = 4)$')
ylabel('$\partial f / \partial \theta_{2}$ $(\mu_{2} = 4)$')
xlim([-2.3 2.3])
ylim([-2.3 2.3])

%% Figure 3.5 - Confidence Ellipses (Optimality Criteria)
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
f5 = figure;
f5.Units = 'centimeters';
f5.Position = [8 4 14 8];
plot(x,y,'LineWidth',1.5,'Color',[0.2431 0.2667 0.2980]);

hold on

fill(x,y, [1.0000 0.8353 0.0000],'FaceAlpha', 0.3)

plot(rect_x,rect_y, '--','LineWidth',1,'Color',[0.0000 0.7451 1.0000])

quiver(h,k,u1(1),u1(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.3176 0.6196],'LineWidth',1);
quiver(h,k,u2_norm(1),u2_norm(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.7451 1.0000],'LineWidth',1);
plot(h,k, 'o','MarkerSize',4,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);

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

%% Figure 4.1 - Excitation Signals (PRBS & PRMS)
f6 = figure;
f6.Units = 'centimeters';
f6.Position = [8 4 14 8];
t3 = tiledlayout(2,1);
t3.TileSpacing = 'compact';
t3.Padding = 'compact';
obsv_input = 100;

nexttile
u_prbs = idinput(obsv_input);
plot(u_prbs, 'LineWidth',1.5)
xlabel('t')
ylabel('u')
ylim([-1.3 1.3])
yticks([-1 1])
xticklabels({})
yticklabels({'$u_{min}$','$u_{max}$'})

nexttile
amp_levels = linspace(-1,1,6);
num_levels = length(amp_levels);
prms_indices = randi(num_levels, [1, obsv_input]);
prms_signal = amp_levels(prms_indices);
plot(prms_signal, 'Linewidth',1.5);
xlabel('t')
ylabel('u')
ylim([-1.2 1.2])
yticks([-1 1])
xticklabels({})
yticklabels({'$u_{min}$','$u_{max}$'})

%% Figure 4.2 - Excitation Signals (Multisine & Chirp)
% t = 0:1:obsv_input;
% f0 = 0; % starting frequency [Hz]
% f1 = 1/2; % end frequency [Hz]
% u_chrip = chirp(t,f0,obsv_input,f1)';
% 
% f7 = figure;
% f7.Units = 'centimeters';
% f7.Position = [8 4 14 8];
% t4 = tiledlayout(2,1);
% t4.TileSpacing = 'compact';
% t4.Padding = 'compact';
% 
% nexttile
% fs = 1;            % Abtastrate
% T = 200;           % Gesamtzeit
% t = 0:1/fs:T-1;    % Zeitvektor
% 
% freq_sin = [0.1, 1, 5];
% amp_sin = [1, 1, 1];
% phase_sin = [0, pi/4, pi/2];
% 
% signal = zeros(T,1);
% for i = 1:length(freq_sin)
%     u_multisine = u_multisine + amp_sin(i) * sin(2*pi*freq_sin(i)*t + phase_sin(i));
% end
% plot(u_multisine, 'LineWidth',1.5)
% xlabel('t')
% ylabel('u')
% xlim([0 200])
% ylim([-1.3 1.3])
% yticks([-1 1])
% xticklabels({})
% yticklabels({'$u_{min}$','$u_{max}$'})
% 
% nexttile
% plot(u_chrip, 'LineWidth',1.5)
% xlabel('t')
% ylabel('u')
% xlim([0 200])
% ylim([-1.3 1.3])
% yticks([-1 1])
% xticklabels({})
% yticklabels({'$u_{min}$','$u_{max}$'})

%% Figure 5.1 -- Active Designed Signal vs. PRBS
% Second-Order OE-System
N = 50;
theta = [-1.3 0.6 0.15];

%% Figure 5.2 -- Uncertainty Ellipses


%% Export Graphics
if export == 1
    f1name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_1.pdf');
    exportgraphics(f1,f1name,'ContentType','vector');
    f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_2.pdf');
    exportgraphics(f2,f2name,'ContentType','vector');
    f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_3.pdf');
    exportgraphics(f3,f3name,'ContentType','vector');
    f4name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_4.pdf');
    exportgraphics(f4,f4name,'ContentType','vector');
    f5name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_5.pdf');
    exportgraphics(f5,f5name,'ContentType','vector');

end