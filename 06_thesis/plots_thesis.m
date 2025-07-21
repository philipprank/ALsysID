%%================================%%
%%%%   Thesis Plots             %%%%
%%================================%%
clearvars
layout_options;
export = 0; % 0 or 1 (for exporting plots to PDF)

obsv = 1000;
eval = floor(obsv/20);

%% Figure 3.1 - Probability Density Functions
f3_2 = figure;
f3_2.Units = 'centimeters';
f3_2.Position = [8 4 14 8];

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

f3_2 = figure;
f3_2.Units = 'centimeters';
f3_2.Position = [8 4 14 12];
t3_2 = tiledlayout(3,2);
t3_2.TileSpacing = 'compact';
t3_2.Padding = 'compact';

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
f3_3 = figure;
f3_3.Units = 'centimeters';
f3_3.Position = [8 4 14 8];
t3_3 = tiledlayout(2,2);
t3_3.TileSpacing = 'compact';
t3_3.Padding = 'compact';

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

f3_4 = figure;
f3_4.Units = 'centimeters';
f3_4.Position = [8 4 14 12];
t3_4 = tiledlayout(3,2);
t3_4.TileSpacing = 'compact';
t3_4.Padding = 'compact';

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
f3_5 = figure;
f3_5.Units = 'centimeters';
f3_5.Position = [8 4 11 11/1.78];
plot(x,y,'LineWidth',1.5,'Color',[0.2431 0.2667 0.2980]);

hold on

fill(x,y, [1.0000 0.8353 0.0000],'FaceAlpha', 0.3)

plot(rect_x,rect_y, '--','LineWidth',1,'Color',[0.0000 0.7451 1.0000])

quiver(h,k,u1(1),u1(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.3176 0.6196],'LineWidth',1);
quiver(h,k,u2_norm(1),u2_norm(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.7451 1.0000],'LineWidth',1);
plot(h,k, 'o','MarkerSize',4,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);

% Labels
text(4.5,-2, 'A-optimality','FontSize',8, 'HorizontalAlignment','center', 'Color',[1.0000 0.8353 0.0000]);
text(-0.4,4.6, 'D-optimality','FontSize',8, 'HorizontalAlignment','center', 'Color',[0.0000 0.7451 1.0000]);
text(5.1,2.8, 'E-optimality','FontSize',8 ,'HorizontalAlignment','center', 'Color',[0.0000 0.3176 0.6196]);

% Format
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xlim([-2.5 6.5])
ylim([-3.7 5.7])
xlabel('$\theta_{1}$', 'FontSize',8)
ylabel('$\theta_{2}$', 'FontSize',8)
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

%% Figure 3.6 - Convex Function
f = @(x) x.^2 - 4*x + 1;
x = linspace(-1,5,500);

x1 = 0.5;
x2 = 4;
y1 = f(x1);
y2 = f(x2);

f3_6 = figure;
f3_6.Units = 'centimeters';
f3_6.Position = [8 4 14 8];
plot(x,f(x), 'LineWidth',1.5);
hold on;

plot(x1, y1, 'o','MarkerSize',4,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);
plot(x2, y2, 'o','MarkerSize',4,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);
plot([x1 x2], [y1 y2], '--', 'LineWidth',1);
text(x1-0.15,y1-0.15,{'$(x_{1}|f(x_{1}))$'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(x2,y2+0.1,{'$(x_{2}|f(x_{2}))$'}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

xlim([-1 5]);
ylim([-1 10]);
hold off;

% Format
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
xlim([-2.5 6.5])
ylim([-3.7 5.7])
xlabel('x');
ylabel('y');
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
f4_1 = figure;
f4_1.Units = 'centimeters';
f4_1.Position = [8 4 14 8];
t4_1 = tiledlayout(2,1);
t4_1.TileSpacing = 'compact';
t4_1.Padding = 'compact';
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
% f8 = figure;
% f8.Units = 'centimeters';
% f8.Position = [8 4 14 8];
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


%% Active Designed Signal vs. PRBS
tic
MC = 10;
N = 200;
a = [1,-0.7];
b = [0,0.1];
init_sys = idtf([0 NaN],[1 NaN],1);
init_sys.Structure.Numerator.Free = [0 1];
sys = idpoly(a,b,1,1,1,[],1);
data_id_prbs = cell(MC,1);
data_id_opt = cell(MC,1);
na = length(sys.A) - 1;
nb = length(sys.B) - 1;
sys_id_prbs = zeros(1,na+nb+2,MC);
sys_id_opt = zeros(1,na+nb+2,MC);
y_prbs  = zeros(N,1);
y_opt = zeros(N,1);
e = 0.1*randn(N,MC);
c_amp = 1;

for i = 1:MC
    u_prbs = idinput(N);
    N_vec = 0:1:N-1;
    % u_prbs = sqrt(2)*sin(2*pi*0.1*N_vec)';
    % Period = 10;
    % NumPeriod = 20;
    % u_prbs = idinput([Period,NumChannel,NumPeriod]);
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxIterations', 5000,...
        'MaxFunctionEvaluations',30000);
    [u_opt, eval(i)] = fmincon(@cost,idinput(N),[],[],[],[],-c_amp*ones(N,1),c_amp*ones(N,1),[],options);

    % Parameter Estimation
    for j = max(na,nb)+1:N
        y_prbs(j) = sim_system(sys,u_prbs(j-1:-1:j-nb,1),y_prbs(j-1:-1:j-na,1),e(j,i));
        y_opt(j) = sim_system(sys,u_opt(j-1:-1:j-nb,1),y_opt(j-1:-1:j-na,1),e(j,i));
    end
    data_id_prbs{i,1} = iddata(y_prbs(1:N,1),u_prbs(1:N,1),1);
    data_id_opt{i,1} = iddata(y_opt(1:N,1),u_opt(1:N,1),1);
    id_struct_prbs = oe(data_id_prbs{i,1},init_sys);
    id_struct_opt = oe(data_id_opt{i,1},init_sys);
    sys_id_prbs(1,:,i) = [id_struct_prbs.B id_struct_prbs.F];
    sys_id_opt(1,:,i) = [id_struct_opt.B id_struct_opt.F];
end

f3_2 = figure;
f3_2.Units = 'centimeters';
f3_2.Position = [8 4 11 11/1.78];
t3_2 = tiledlayout(2,1);
t3_2.TileSpacing = 'compact';
t3_2.Padding = 'compact';

nexttile
plot(1:N,u_prbs, 'LineWidth',1.5)
ax = gca;
ax.FontSize = 8; 
ylim([-1.2 1.2])
% ylabel ('u')
yticks([-1 1])
% yticklabels({'$-c_{amp}$','$c_{amp}$'})

nexttile
plot(1:N,u_opt, 'LineWidth',1.5)
ylim([-1.2 1.2])
xlabel('Observations T', 'FontSize',8)
% ylabel ('u')
yticks([-1 1])
% yticklabels({'$-c_{amp}$','$c_{amp}$'})

function obj = cost(u)
% parameters
N = 200;
a1 = -0.7;
b1 = 0.1;
num_params = 2;

y = zeros(1,N);
dy_da1 = zeros(1,N);
dy_db1 = zeros(1,N);
dy_theta = zeros(num_params,N);

% Initial conditions
y(1) = 0;
y(2) = 0;
dy_da1(1) = 0;
dy_da1(2) = 0;
dy_db1(1) = 0;
dy_db1(2) = 0;

% Recursive expressions
for t = 3:N
    y(t) = -a1*y(t-1) + b1*u(t-1);

    dy_da1(t) = -y(t-1) - a1*dy_da1(t-1);
    dy_db1(t) = u(t-1) - a1*dy_db1(t-1);
    dy_theta(:,t) = [dy_da1(t); dy_db1(t)];
end

% Fisher information matrix & cost function
I = zeros(num_params,num_params);
for t = 1:N
    I = I + dy_theta(:,t)*dy_theta(:,t)';
end
% Optimality Criteria
obj = -log(det(I));
% obj = trace(I^(-1));
% obj = -min(eig(I));
end
toc

%% Scatter Plot & Uncertainty Ellipses (Active Designed Signal vs. PRBS)
% Uncertainty Ellipse
function plot_ellipse(data)
data_est = squeeze([data(:,2,:), data(:,4,:)])';
cova = cov(data_est);
[eigenvec, eigenval] = eig(cova);
[largest_eigenvec_ind_c,~] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
largest_eigenval = max(max(eigenval));

if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
else
    smallest_eigenval = max(eigenval(:,1));
end

angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
if(angle < 0)
    angle = angle + 2*pi;
end
avg = mean(data_est);

chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0 = avg(1);
Y0 = avg(2);
a = chisquare_val*sqrt(largest_eigenval);
b = chisquare_val*sqrt(smallest_eigenval);

ellipse_x_r  = a*cos(theta_grid);
ellipse_y_r  = b*sin(theta_grid);
R = [ cos(phi) sin(phi); -sin(phi) cos(phi)];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0, 'LineWidth',1.5, 'Color',[0.0000 0.3176 0.6196])
end

% Plots
f3_3 = figure;
f3_3.Units = 'centimeters';
f3_3.Position = [8 4 11 11/1.78];
t3_3 = tiledlayout(1,2);
t3_3.TileSpacing = 'compact';
t3_3.Padding = 'compact';

nexttile
for i = 1:MC
    plot(sys_id_prbs(:,2,i), sys_id_prbs(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_prbs);
xlabel('$\hat{\theta}_{1}$')
ylabel('$\hat{\theta}_{2}$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

nexttile
for i = 1:MC
    plot(sys_id_opt(:,2,i), sys_id_opt(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_opt);
xlabel('$\hat{\theta}_{1}$')
ylabel('$\hat{\theta}_{2}$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

f5name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_5.pdf');
exportgraphics(f3_5,f5name,'ContentType','vector');

%% Export Figures
f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_3.pdf');
exportgraphics(f3_2,f2name,'ContentType','vector');
% f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_4.pdf');
% exportgraphics(f3_3,f3name,'ContentType','vector');
%% Figure 5.3 -- Convex Relaxed Signal
% CVX needs to be installed and on the path
tic
MC = 300;
N = 200;
a = [1,-0.7];
b = [0,0.1];
a1 = -0.7;
b1 = 0.1;
init_sys = idtf([0 NaN],[1 NaN],1);
init_sys.Structure.Numerator.Free = [0 1];
sys = idpoly(a,b,1,1,1,[],1);
data_id_prbs = cell(MC,1);
data_id_opt = cell(MC,1);
na = length(sys.A) - 1;
nb = length(sys.B) - 1;
sys_id_prbs = zeros(1,na+nb+2,MC);
sys_id_opt = zeros(1,na+nb+2,MC);
y_prbs  = zeros(N,1);
y_opt = zeros(N,1);
e = 0.1*randn(N,MC);
c_amp = 1;

% Convex Optimization Problem
elem_a = zeros(1,N);
elem_b = zeros(1,N);
for j = 1:N-1
    elem_a(j+1) = (-1)^(j+1)*(j-1)*a1^(j-2)*b1;
    elem_b(j+1) = (-1)^(j+1)*a1^(j-1);
end
F1 = toeplitz(elem_a,zeros(1,N));
F2 = toeplitz(elem_b,zeros(1,N));
Q_11 = F1'*F1;
Q_12 = F1'*F2;
Q_21 = F2'*F1;
Q_22 = F2'*F2;

cvx_begin
variable U(N,N) semidefinite
minimize(trace(inv([trace(U*Q_11) trace(U*Q_12); trace(U*Q_21) trace(U*Q_22)])))
subject to
for k = 1:N
    U(k,k) <= c_amp;
end
cvx_end

% [V, D] = eig(U);
% eigenvalues = diag(D);
% [~, maxIndex] = max(eigenvalues);
% u_opt = V(:, maxIndex);

for i = 1:MC
    u_prbs = idinput(N);
    % Optimal Input Signal Generation
    D = chol(U,'upper');
    Xi  = randn(N,1);
    u_opt = diag(1)*sign(D'*Xi);

    % Parameter Estimation
    for j = max(na,nb)+1:N
        y_prbs(j) = sim_system(sys,u_prbs(j-1:-1:j-nb,1),y_prbs(j-1:-1:j-na,1),e(j,i));
        y_opt(j) = sim_system(sys,u_opt(j-1:-1:j-nb,1),y_opt(j-1:-1:j-na,1),e(j,i));
    end
    data_id_prbs{i,1} = iddata(y_prbs(1:N,1),u_prbs(1:N,1),1);
    data_id_opt{i,1} = iddata(y_opt(1:N,1),u_opt(1:N,1),1);
    id_struct_prbs = oe(data_id_prbs{i,1},init_sys);
    id_struct_opt = oe(data_id_opt{i,1},init_sys);
    sys_id_prbs(1,:,i) = [id_struct_prbs.B id_struct_prbs.F];
    sys_id_opt(1,:,i) = [id_struct_opt.B id_struct_opt.F];
end

f5_3 = figure;
f5_3.Units = 'centimeters';
f5_3.Position = [8 4 14 8];
t5_3 = tiledlayout(2,1);
t5_3.TileSpacing = 'compact';
t5_3.Padding = 'compact';

nexttile
plot(1:N,u_prbs, 'LineWidth',1.5)
ylim([-1.2 1.2])
% ylabel ('u')
yticks([-1 1])
yticklabels({'$-c_{amp}$','$c_{amp}$'})

nexttile
plot(1:N,u_opt, 'LineWidth',1.5)
ylim([-1.2 1.2])
xlabel('Observations N')
% ylabel ('u')
yticks([-1 1])
yticklabels({'$-c_{amp}$','$c_{amp}$'})
toc

%% Figure 5.4 -- Scatter Plot & Uncertainty Ellipses (Active Designed Signal vs. PRBS)
% Uncertainty Ellipse
f5_4 = figure;
f5_4.Units = 'centimeters';
f5_4.Position = [8 4 14 8];
t5_4 = tiledlayout(1,2);
t5_4.TileSpacing = 'compact';
t5_4.Padding = 'compact';

nexttile
for i = 1:MC
    plot(sys_id_prbs(:,2,i), sys_id_prbs(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_prbs);
xlabel('$a$')
ylabel('$b$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

nexttile
for i = 1:MC
    plot(sys_id_opt(:,2,i), sys_id_opt(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_opt);
xlabel('$a$')
ylabel('$b$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

%% Export Figure 5.3 & Figure 5.4
f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_3.pdf');
exportgraphics(f5_3,f2name,'ContentType','vector');
f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_4.pdf');
exportgraphics(f5_4,f3name,'ContentType','vector');


%% Figure 5.5 -- Receding Horizon Input Design
tic
% load u_opt.mat u_opt
% load e.mat e
MC = 300;
N = 200;
H = 8; % horizon length
a = [1,-0.7];
b = [0,0.1];
init_sys = idtf([0 NaN],[1 NaN],1);
init_sys.Structure.Numerator.Free = [0 1];
sys = idpoly(a,b,1,1,1,[],1);
data_id_prbs = cell(MC,1);
data_id_opt = cell(MC,1);
na = length(sys.A) - 1;
nb = length(sys.B) - 1;
sys_id_prbs = zeros(1,na+nb+2,MC);
sys_id_opt = zeros(1,na+nb+2,MC);
y_prbs  = zeros(N,1);
u_opt = zeros(N,MC);
u_opt(1,:) = idinput(MC)';
% u_opt = [ones(1,20) -ones(1,10) ones(1,20) -ones(1,30) ones(1,20) ones(1,10) -ones(1,30) ones(1,20) -ones(1,40)]';
y_opt = zeros(N,1);
e = 0.1*randn(N,MC);
c_amp = 1;
J = zeros(MC,1);

M = 0:2^H-1;
u_all = dec2bin(M,H) - '0';
u_all(u_all == 0) = -1;
u_all = flip(u_all);


%%
for i = 1:MC
    u_prbs = idinput(N);
    N_vec = 0:1:N-1;
    I_p = zeros(2,2);

    % Parameter Estimation
    for j = max(na,nb)+1:N
        u_hat = receding_horizon_cost(I_p,M,u_all,u_opt(j-1,i),y_opt(j-1),H,a,b);
        u_test(:,j) = u_hat;
        u_opt(j,i) = u_hat(2);
        y_prbs(j) = sim_system(sys,u_prbs(j-1:-1:j-nb,1),y_prbs(j-1:-1:j-na,1),e(j,i));
        y_opt(j) = sim_system(sys,u_opt(j-1:-1:j-nb,i),y_opt(j-1:-1:j-na,1),e(j,i));
        I_p = I_p + [-y_opt(j); u_opt(j,i)]*[-y_opt(j); u_opt(j,i)]';
    end
    J(i) = log(det(I_p));
    data_id_prbs{i,1} = iddata(y_prbs(1:N,1),u_prbs(1:N,1),1);
    data_id_opt{i,1} = iddata(y_opt(1:N,1),u_opt(1:N,i),1);
    id_struct_prbs = oe(data_id_prbs{i,1},init_sys);
    id_struct_opt = oe(data_id_opt{i,1},init_sys);
    sys_id_prbs(1,:,i) = [id_struct_prbs.B id_struct_prbs.F];
    sys_id_opt(1,:,i) = [id_struct_opt.B id_struct_opt.F];
end

f5_5 = figure;
f5_5.Units = 'centimeters';
f5_5.Position = [8 4 14 8];
t5_5 = tiledlayout(2,1);
t5_5.TileSpacing = 'compact';
t5_5.Padding = 'compact';

nexttile
plot(1:N, u_prbs, 'LineWidth',1.5)
ylim([-1.2 1.2])
yticks([-1 1])
yticklabels({'$-c_{amp}$','$c_{amp}$'})

nexttile
plot(1:N, u_opt(:,1), 'LineWidth', 1.5)
ylim([-1.2 1.2])
xlabel('Observations N')
yticks([-1 1])
yticklabels({'$-c_{amp}$','$c_{amp}$'})

function u_opt_H = receding_horizon_cost(I_p,M,u_all,u_p,y_p,H,a,b)
a1 = a(2);
b1 = b(2);
y = zeros(1,H-1);
dy_da1 = zeros(1,H);
dy_db1 = zeros(1,H);
dy_theta = zeros(2,H);

% Initial conditions
y(1) = y_p;
dy_da1(1) = 0;
dy_db1(1) = 0;
for k = 1:length(M)
    I_f = zeros(2,2);
    u = [u_p u_all(k,:)];
    for t = 2:H
        y(t) = -a1*y(t-1) + b1*u(t-1);

        dy_da1(t) = -y(t-1) - a1*dy_da1(t-1);
        dy_db1(t) = u(t-1) - a1*dy_db1(t-1);
        dy_theta(:,t) = [dy_da1(t); dy_db1(t)];
        I_f = I_f + dy_theta(:,t)*dy_theta(:,t)';
    end
    I = I_p + I_f;
    if k == 1
        u_opt_H = u;
        comp = I; % comp: variable for internal comparison
    else
        if log(det(I)) > log(det(comp))
            u_opt_H = u;
            comp = I;
        end
    end
end
end

toc

%% Figure 5.6 -- Scatter Plot & Uncertainty Ellipses (Receding Horizon vs. PRBS)
% Uncertainty Ellipse
f5_6 = figure;
f5_6.Units = 'centimeters';
f5_6.Position = [8 4 14 8];
t5_6 = tiledlayout(1,2);
t5_6.TileSpacing = 'compact';
t5_6.Padding = 'compact';

nexttile
for i = 1:MC
    plot(sys_id_prbs(:,2,i), sys_id_prbs(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_prbs);
xlabel('$a$')
ylabel('$b$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

nexttile
for i = 1:MC
    plot(sys_id_opt(:,2,i), sys_id_opt(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_opt);
xlabel('$a$')
ylabel('$b$')
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

%%
% Gegebenes Signal
Fs = 1; % Abtastfrequenz
signal = randn(1,200); % Beispielsignal (hier ein zufälliges Signal)

% Berechnung der spektralen Leistungsdichte mit pwelch
[pxx, f] = pwelch(signal,[],[],[],Fs);

% Plotten der spektralen Leistungsdichte
figure;
plot(f, 10*log10(pxx)); % PSD in dB/Hz
title('Spektrale Leistungsdichte des Signals');
xlabel('Frequenz (Hz)');
ylabel('Leistungsdichte (dB/Hz)');
grid on;

%% Figure 5.6 -- Frequency Domain Input Design
tic
N = 20;
a1 = -0.7;
b1 = 0.1;
c_pow = 1;
cvx_begin
variables r0 r1 r2
minimize(-log_det([b1^(2)*r0 -a1*b1*r0-b1*r1; -a1*b1*r0-b1*r1 r0+2*a1*r1+a1^2*r2]))
subject to
[r0 r1 r2; r1 r0 r1; r2 r1 r0] == semidefinite(3);
r0 <= c_pow;
cvx_end
toc


%% Export Graphics
if export == 1
    f1name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_1.pdf');
    exportgraphics(f3_3,f1name,'ContentType','vector');
    f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_2.pdf');
    exportgraphics(f3_3,f2name,'ContentType','vector');
    f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_3.pdf');
    exportgraphics(f3_3,f3name,'ContentType','vector');
    f4name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_4.pdf');
    exportgraphics(f3_4,f4name,'ContentType','vector');
    f5name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_5.pdf');
    exportgraphics(f3_5,f5name,'ContentType','vector');
    f6name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_3_6.pdf');
    exportgraphics(f3_6,f6name,'ContentType','vector');
    f7name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_1.pdf');
    exportgraphics(f3_3,f7name,'ContentType','vector');
    f8name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_2.pdf');
    exportgraphics(f3_3,f8name,'ContentType','vector');
    f2name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_3.pdf');
    exportgraphics(f5_3,f2name,'ContentType','vector');
    f3name = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_5_4.pdf');
    exportgraphics(f5_4,f3name,'ContentType','vector');
end