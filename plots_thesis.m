layout_options
clear; clc; close all;
rng(8);
%% Figure 2.2 - Confidence Ellipses (Optimality Criteria)
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
u1_wn = [a*cos(theta); a*sin(theta)];
u2_wn = [-b*sin(theta); b*cos(theta)];

major_x = [h, h+u1_wn(1)]; % endpoints
major_y = [k, k+u1_wn(2)];
minor_length = sqrt((x_max-x_min)^2 + (y_max-y_min)^2)/2;
u2_norm = u2_wn/norm(u2_wn)*minor_length;
minor_x = [h, h+u2_norm(1)];
minor_y = [k, k+u2_norm(2)];

% Plot
f2_2 = figure;
f2_2.Units = 'centimeters';
f2_2.Position = [8 4 11 11/2];
plot(x,y,'LineWidth',1,'Color',[0.2431 0.2667 0.2980]);

hold on
fill(x,y, [1.0000 0.8353 0.0000],'FaceAlpha', 0.5)
plot(rect_x,rect_y, '--','LineWidth',1,'Color',[0.0000 0.7451 1.0000])

quiver(h,k,u1_wn(1),u1_wn(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.3176 0.6196],'LineWidth',1);
quiver(h,k,u2_norm(1),u2_norm(2),0, 'MaxHeadSize',0.25,'Color',[0.0000 0.7451 1.0000],'LineWidth',1);
plot(h,k, 'o','MarkerSize',3,'MarkerFaceColor',[0.2431 0.2667 0.2980],'MarkerEdgeColor',[0.2431 0.2667 0.2980]);

% Labels
text(4.5,-2, 'A-optimality','FontSize',8, 'HorizontalAlignment','center', 'Color',[1.0000 0.8353 0.0000]);
text(-0.4,4.6, 'D-optimality','FontSize',8, 'HorizontalAlignment','center', 'Color',[0.0000 0.7451 1.0000]);
text(5.1,2.8, 'E-optimality','FontSize',8, 'HorizontalAlignment','center', 'Color',[0.0000 0.3176 0.6196]);

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
annotation('arrow', [xs xe],[ys ys], 'LineWidth',0.8,'HeadStyle','vback3',...
    'HeadWidth',6,'HeadLength',6);
annotation('arrow', [xs xs],[ys ye], 'LineWidth',0.8,'HeadStyle','vback3',...
    'HeadWidth',6,'HeadLength',6);
hold off

%% Export Figure 2.2 - PDF & TikZ
pdfFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_2.pdf');
exportgraphics(f2_2, pdfFile, 'ContentType', 'vector');
% Export mit matlab2tikz:
tikzFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_2.tex');
matlab2tikz(tikzFile, 'figurehandle', f2_2, 'standalone', true, 'extraPreamble', '\usepackage{lmodern}');

%% Active Designed Signal vs. PRBS
tic
MC = 300;
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

%% Figure 2.3 - Scatter Plot & Uncertainty Ellipses
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

plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0, 'LineWidth',1, 'Color',[0.0000 0.3176 0.6196])
end

% Plots
f2_3 = figure;
f2_3.Units = 'centimeters';
f2_3.Position = [8 4 11 11/1.78];
t2_3 = tiledlayout(1,2);
t2_3.TileSpacing = 'compact';
t2_3.Padding = 'compact';

nexttile
for i = 1:MC
    plot(sys_id_prbs(:,2,i), sys_id_prbs(:,4,i),'.','Color',[0.2431 0.2667 0.2980],...
        'MarkerSize',5)
    hold all
end
plot_ellipse(sys_id_prbs);
ax = gca;
ax.FontSize = 8;
xlabel('$\hat{\theta}_{1}$', 'FontSize',8)
ylabel('$\hat{\theta}_{2}$', 'FontSize',8)
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
ax = gca;
ax.FontSize = 8;
xlabel('$\hat{\theta}_{1}$', 'FontSize',8)
ylabel('$\hat{\theta}_{2}$', 'FontSize',8)
xlim([0.05 0.15])
yticks([-0.9 -0.8 -0.7 -0.6 -0.5])
ylim([-0.9 -0.5])

%% Export Figure 2.3 - PDF & TikZ
pdfFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_3.pdf');
exportgraphics(f2_3, pdfFile, 'ContentType', 'vector');
% Export mit matlab2tikz:
tikzFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_3.tex');
matlab2tikz(tikzFile, 'figurehandle', f2_3, 'standalone', true, 'extraPreamble', '\usepackage{lmodern}');

%% Figure 2.4 - Comparison PRBS vs. Designed Input Signal
f2_4 = figure;
f2_4.Units = 'centimeters';
f2_4.Position = [8 4 11 11/1.78];
t2_4 = tiledlayout(2,1);
t2_4.TileSpacing = 'compact';
t2_4.Padding = 'compact';

nexttile
plot(1:N,u_prbs, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
ylim([-1.2 1.2])
% ylabel ('u')
yticks([-1 1])
% yticklabels({'$-c_{amp}$','$c_{amp}$'})

nexttile
plot(1:N,u_opt, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
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

%% Export Figure 2.4 - PDF & TikZ
pdfFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_4.pdf');
exportgraphics(f2_4, pdfFile, 'ContentType', 'vector');
% Export mit matlab2tikz:
tikzFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_2_4.tex');
matlab2tikz(tikzFile, 'figurehandle', f2_4, 'standalone', true, 'extraPreamble', '\usepackage{lmodern}');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beispiel: FIR-Filter mit 10 Koeffizienten passend zu r0, r1, r2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;
rng(8);  % reproducability

%% Parameter
a = -0.7;                  % OE-Modellparameter
b = 0.1;                   % OE-Modellparameter
sigma_epsilon = 0.1;       % Störvarianz
T = 1000;                 % Länge des Zeitsignals
K = 3;                    % Anzahl der FIR-Koeffizienten (h0...h9)

%% CVX-Optimierung: Bestimme optimale Autokorrelationswerte r0, r1, r2
cvx_begin quiet
variables r0 r1 r2

% Fisher-Informationsmatrix für OE-Modell
M = [ b^2*r0,           -a*b*r0 - b*r1;
    -a*b*r0 - b*r1,  (1 + a^2)*r0 + 2*a*r1 ];

% Zielfunktion: D-Optimalität (maximiere log_det(M))
minimize(-log_det(M))

% Energie-Beschränkung
2*a^2*r2 + 2*a*(1+a^2)*r1 + (1 + a^4 + 4*a^2)*r0 <= 1;

% PSD-Bedingung: Toeplitzmatrix der Autokorrelation positiv semidefinit
T_r = [ r0, r1, r2;
    r1, r0, r1;
    r2, r1, r0 ];
T_r == semidefinite(3);
cvx_end

fprintf('Optimierte Autokorrelationswerte:\n');
fprintf('r0 = %.4f, r1 = %.4f, r2 = %.4f\n\n', r0, r1, r2);

%% Ziel-ACF
r0 = 1.0;
r1 = -0.5;
r2 = -0.8;

%% Yule-Walker-Gleichung lösen (für AR(2))
R = [r0, r1;
    r1, r0];
r_vec = [r1; r2];

a = R \ r_vec;
a1 = a(1);
a2 = a(2);
fprintf('Analytisch berechnete AR(2)-Koeffizienten:\n');
fprintf('a1 = %.4f, a2 = %.4f\n', a1, a2);

%% Rauschvarianz berechnen (damit ACF bei lag 0 = r0 passt)
sigma_eta_sq = r0 - a1*r1 - a2*r2;
fprintf('Benötigte Varianz des weißen Rauschens: %.4f\n', sigma_eta_sq);

%% AR(2)-Signal erzeugen
T = 10000;
eta = sqrt(sigma_eta_sq) * randn(T,1);
u = filter(1, [1, -a1, -a2], eta);

%% Zeitsignal & ACF
subplot(2,1,1);
plot(u); title('Generiertes AR(2)-Signal'); grid on;

[acf_est, lags] = xcorr(u, 20, 'unbiased');
acf_pos = acf_est(lags >= 0);
lags_pos = lags(lags >= 0);

subplot(2,1,2);
stem(lags_pos, acf_pos, 'filled');
xlabel('Lag'); ylabel('ACF');
title('Empirisch geschätzte Autokorrelationsfunktion'); grid on;

fprintf('\nEmpirische ACF-Werte:\n');
for k = 0:2
    fprintf('r_%d ≈ %.4f\n', k, acf_pos(k+1));
end

%%
clear; clc; close all;

%% Systemdefinition aus Bild
% A-Polynome
A11 = [1, -1.9659, 0.9663];
A12 = [0, 0.0471, -0.0455];
A21 = [0, 0.0058, -0.0068];
A22 = [1, -1.9228, 0.9250];

% B-Polynome
B11 = [0, 42.1906, -40.0763];
B12 = [0, 89.4773, -88.2050];
B21 = [0, 73.5355, -71.5468];
B22 = [0, 79.0570, -78.9397];

T = 1000;                        % Anzahl Zeitschritte
u = idinput(T, 'prbs', [0 20], [-1 1]);                 % Weißes Rauschen: u(:,1)=u1, u(:,2)=u2
y = zeros(T, 2);                 % Initialisierung für y1 und y2

for t = 3:T
    % Output 1 berechnen
    y(t,1) = -A11(2)*y(t-1,1) - A11(3)*y(t-2,1) ...
        -A12(2)*y(t-1,2) - A12(3)*y(t-2,2) ...
        +B11(2)*u(t-1,1) + B11(3)*u(t-2,1) ...
        +B12(2)*u(t-1,2) + B12(3)*u(t-2,2);

    % Output 2 berechnen
    y(t,2) = -A22(2)*y(t-1,2) - A22(3)*y(t-2,2) ...
        -A21(2)*y(t-1,1) - A21(3)*y(t-2,1) ...
        +B21(2)*u(t-1,1) + B21(3)*u(t-2,1) ...
        +B22(2)*u(t-1,2) + B22(3)*u(t-2,2);
end

figure;
subplot(2,1,1);
plot(y(:,1), 'b');
title('Output y_1');

subplot(2,1,2);
plot(y(:,2), 'r');
title('Output y_2');

%%
T = 4000;
Ts = 1;

% Weißes Rauschen
u_white = randn(T,2)*0.0005;

% Gefiltertes langsames Rauschen
[b,a] = butter(2, 0.01);  % Tiefpass
u_slow = filter(b,a,u_white);

% Plot
figure;
subplot(2,1,1); plot(u_white); title('Standard Weißes Rauschen');
subplot(2,1,2); plot(u_slow); title('Gefiltertes (langsames) Rauschen');

%%
clear; clc; close all;
layout_options

%% Systemdefinition – ARX-Polynome
% y1-Output
A11 = [1 -1.9659 0.9663];
A12 = [0 0.0471 -0.0455];
B11 = [0 42.1906 -40.0763];
B12 = [0 89.4773 -88.2050];

% y2-Output
A21 = [0 0.0058 -0.0068];
A22 = [1 -1.9228 0.9250];
B21 = [0 73.5355 -71.5468];
B22 = [0 79.0570 -78.9397];

%% Simulationseinstellungen
T = 800;
u1_wn = zeros(T,1);
u2_wn = zeros(T,1);
y1 = zeros(T,1);
y2 = zeros(T,1);

%% Step Response: u1 = 1, u2 = 0.0005
u1_wn(:) = 0.0005;
u2_wn(:) = 0;
y1_u1 = zeros(T,1);
y2_u1 = zeros(T,1);
y1 = zeros(T,1); y2 = zeros(T,1);

for t = 3:T
    y1(t) = -A11(2)*y1(t-1) - A11(3)*y1(t-2) ...
        -A12(2)*y2(t-1) - A12(3)*y2(t-2) ...
        +B11(2)*u1_wn(t-1) + B11(3)*u1_wn(t-2) ...
        +B12(2)*u2_wn(t-1) + B12(3)*u2_wn(t-2);
    y2(t) = -A22(2)*y2(t-1) - A22(3)*y2(t-2) ...
        -A21(2)*y1(t-1) - A21(3)*y1(t-2) ...
        +B21(2)*u1_wn(t-1) + B21(3)*u1_wn(t-2) ...
        +B22(2)*u2_wn(t-1) + B22(3)*u2_wn(t-2);
    y1_u1(t) = y1(t);
    y2_u1(t) = y2(t);
end

%% Step Response: u1 = 0, u2 = 0.0005
u1_wn(:) = 0;
u2_wn(:) = 0.0005;
y1_u2 = zeros(T,1);
y2_u2 = zeros(T,1);
y1 = zeros(T,1); y2 = zeros(T,1);

for t = 3:T
    y1(t) = -A11(2)*y1(t-1) - A11(3)*y1(t-2) ...
        -A12(2)*y2(t-1) - A12(3)*y2(t-2) ...
        +B11(2)*u1_wn(t-1) + B11(3)*u1_wn(t-2) ...
        +B12(2)*u2_wn(t-1) + B12(3)*u2_wn(t-2);
    y2(t) = -A22(2)*y2(t-1) - A22(3)*y2(t-2) ...
        -A21(2)*y1(t-1) - A21(3)*y1(t-2) ...
        +B21(2)*u1_wn(t-1) + B21(3)*u1_wn(t-2) ...
        +B22(2)*u2_wn(t-1) + B22(3)*u2_wn(t-2);
    y1_u2(t) = y1(t);
    y2_u2(t) = y2(t);
end

%% Figure - Step Response
f4_2 = figure;
f4_2.Units = 'centimeters';
f4_2.Position = [8 4 11 11/1.78];

subplot(2,2,1)
plot(0:T-1, y1_u1, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
yticks([0 0.4 0.8])
xlim([0, 800])
ylim([0, 0.8])
title('$u_{1} \to y_{1}$', 'FontSize',8)
xlabel('Time $T$ [s]', 'FontSize',8);
ylabel('$y_{1}$ [m]', 'FontSize',8);
grid on;

subplot(2,2,2)
plot(0:T-1, y2_u1, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
yticks([0 0.4 0.8])
xlim([0, 800])
ylim([0, 0.8])
title('$u_{1} \to y_{2}$', 'FontSize',8)
xlabel('Time $T$ [s]', 'FontSize',8);
ylabel('$y_{2}$ [m]', 'FontSize',8);
grid on;

subplot(2,2,3)
plot(0:T-1, y1_u2, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
xlim([0, 800])
ylim([0, 1])
title('$u_{2} \to y_{1}$', 'FontSize',8);
xlabel('Time $T$ [s]', 'FontSize',8);
ylabel('$y_{1}$ [m]', 'FontSize',8);
grid on;

subplot(2,2,4)
plot(0:T-1, y2_u2, 'LineWidth',1)
ax = gca;
ax.FontSize = 8;
xlim([0, 800])
ylim([0, 0.4])
title('$u_{2} \to y_{2}$', 'FontSize',8)
xlabel('Time $T$ [s]', 'FontSize',8);
ylabel('$y_{2}$ [m]', 'FontSize',8);
grid on;

%% Export Figure
f4 = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_4_2.pdf');
exportgraphics(f4_2,f4,'ContentType','vector');

%% =========================================================================
%%  Quadruple-Tank Process
%% =========================================================================
%% Initialization
clear; clc; close all;
rng(54);

%% Parameters
A11 = [1, -1.9659, 0.9663];
A12 = [0, 0.0471, -0.0455];
B11 = [0, 42.1906, -40.0763];
B12 = [0, 89.4773, -88.2050];

A21 = [0, 0.0058, -0.0068];
A22 = [1, -1.9228,  0.9250];
B21 = [0, 73.5355, -71.5468];
B22 = [0, 79.0570, -78.9397];

theta1_true = [A11(2:3), A12(2:3), B11(2:3), B12(2:3)]';
theta2_true = [A22(2:3), A21(2:3), B21(2:3), B22(2:3)]';
theta_true = [theta1_true; theta2_true];

%% Settings
% Simulation Parameters
T = 12000; % Observation length
Ts = 1; % Sample time
M = 200; % Number of simulations
interval = 100:200:T;
numSteps = numel(interval);
noiseStd1 = 0.000005;
noiseStd2 = 0.000005;
e1 = noiseStd1*randn(T,1,M);
e2 = noiseStd2*randn(T,1,M);

% Input Signal Parameters
du = 2; % Input dimension
Range = [-0.0005, 0.0005]; % Amplitude constraints

y1_sim_wn = zeros(T,M);
y2_sim_wn = zeros(T,M);
y_sim_wn = zeros(T,2,M);
u1_wn = zeros(T,M);
u2_wn = zeros(T,M);
u_wn = zeros(T,2,M);

y1_sim_prbs = zeros(T,M);
y2_sim_prbs = zeros(T,M);
y_sim_prbs = zeros(T,2,M);

y1_sim_periodic = zeros(T,M);
y2_sim_periodic = zeros(T,M);
y_sim_periodic = zeros(T,2,M);

y1_sim_oracle = zeros(T,1);
y2_sim_oracle = zeros(T,1);
y_sim_oracle = zeros(T,2);

errorNorm_WN = zeros(M,numSteps);
errorNorm_PRBS = zeros(M,numSteps);
errorNorm_periodic = zeros(M,numSteps);
errorNorm_oracle = zeros(M,numSteps);

% Oracle Input Calculation
A_arx = [-A11(2:end) -A12(2:end); -A21(2:end) -A22(2:end)];
B_arx = [B11(2:end) B12(2:end); B21(2:end) B22(2:end)];
A = [A_arx B_arx; eye(2) zeros(2,6); zeros(2,8); zeros(2,4) eye(2) zeros(2)];
B = [zeros(4,2); eye(2); zeros(2)];
sigma_u = 2.5e-07;
k = 100;
seed = 10;
[u_opt,u_df] = design_input(A,B,k,sigma_u,seed);
u_oracle = repmat(u_opt, ceil(T/k), 1);
u_oracle = u_oracle(1:T,:);
u1_oracle = u_oracle(:,1);
u2_oracle = u_oracle(:,2);

%% Simulation
for runIdx = 1:M
    % White Noise Parameters & Signal Generation
    Band = [0.002, 1];
    u_wn(:,:,runIdx) = idinput([T,du],'rgs',Band,Range); % rbs: 'random binary signal'
    u1_wn = u_wn(:,1,runIdx);
    u2_wn = u_wn(:,2,runIdx);
    for t = 3:T
        y1_sim_wn(t,runIdx) = -A11(2)*y1_sim_wn(t-1,runIdx) - A11(3)*y1_sim_wn(t-2,runIdx) ...
            - A12(2)*y2_sim_wn(t-1,runIdx) - A12(3)*y2_sim_wn(t-2,runIdx) ...
            + B11(2)*u1_wn(t-1) + B11(3)*u1_wn(t-2) + B12(2)*u2_wn(t-1) + B12(3)*u2_wn(t-2);
        y2_sim_wn(t,runIdx) = -A22(2)*y2_sim_wn(t-1,runIdx) - A22(3)*y2_sim_wn(t-2,runIdx) ...
            - A21(2)*y1_sim_wn(t-1,runIdx) - A21(3)*y1_sim_wn(t-2,runIdx) ...
            + B21(2)*u1_wn(t-1) + B21(3)*u1_wn(t-2) + B22(2)*u2_wn(t-1) + B22(3)*u2_wn(t-2);
        y_sim_wn(t,:,runIdx) = [y1_sim_wn(t,runIdx), y2_sim_wn(t,runIdx)];
    end
    
    % Simulation - PRBS
    % PRBS Parameters & Signal Generation
    Band = [0.002, 1];
    u_prbs = idinput([T,du],'prbs',Band,Range); % prbs: 'pseudorandom binary signal'
    u1_prbs = u_prbs(:,1);
    u2_prbs = u_prbs(:,2);
    for t = 3:T
        y1_sim_prbs(t,1,runIdx) = -A11(2)*y1_sim_prbs(t-1,1,runIdx) - A11(3)*y1_sim_prbs(t-2,1,runIdx) ...
            - A12(2)*y2_sim_prbs(t-1,1,runIdx) - A12(3)*y2_sim_prbs(t-2,1,runIdx) ...
            + B11(2)*u1_prbs(t-1) + B11(3)*u1_prbs(t-2) + B12(2)*u2_prbs(t-1) + B12(3)*u2_prbs(t-2);
        y2_sim_prbs(t,1,runIdx) = -A22(2)*y2_sim_prbs(t-1,1,runIdx) - A22(3)*y2_sim_prbs(t-2,1,runIdx) ...
            - A21(2)*y1_sim_prbs(t-1,1,runIdx) - A21(3)*y1_sim_prbs(t-2,1,runIdx) ...
            + B21(2)*u1_prbs(t-1) + B21(3)*u1_prbs(t-2) + B22(2)*u2_prbs(t-1) + B22(3)*u2_prbs(t-2);
        y_sim_prbs(t,:,runIdx) = [y1_sim_prbs(t,1,runIdx), y2_sim_prbs(t,1,runIdx)];
    end
    
    % Simulation - Periodic Signal
    % Multisine Parameters & Signal Generation
    numComp = 10; % Number of frequency components for periodic signal
    freq_limit = [0.001 1]; % Frequency limits for periodic signal
    [u_periodic, info] = multisine_generation(freq_limit,1,T,numComp,du);
    freq = info(1).frequencies;
    currVar = var(u_periodic);
    scale = sqrt(2.5*10^(-7)./currVar);
    u1_periodic = u_periodic(:,1);
    u2_periodic = u_periodic(:,2);
    for t = 3:T
        y1_sim_periodic(t,1,runIdx) = -A11(2)*y1_sim_periodic(t-1,1,runIdx) - A11(3)*y1_sim_periodic(t-2,1,runIdx) ...
            - A12(2)*y2_sim_periodic(t-1,1,runIdx) - A12(3)*y2_sim_periodic(t-2,1,runIdx) ...
            + B11(2)*u1_periodic(t-1) + B11(3)*u1_periodic(t-2) + B12(2)*u2_periodic(t-1) + B12(3)*u2_periodic(t-2);
        y2_sim_periodic(t,1,runIdx) = -A22(2)*y2_sim_periodic(t-1,1,runIdx) - A22(3)*y2_sim_periodic(t-2,1,runIdx) ...
            - A21(2)*y1_sim_periodic(t-1,1,runIdx) - A21(3)*y1_sim_periodic(t-2,1,runIdx) ...
            + B21(2)*u1_periodic(t-1) + B21(3)*u1_periodic(t-2) + B22(2)*u2_periodic(t-1) + B22(3)*u2_periodic(t-2);
        y_sim_periodic(t,:,runIdx) = [y1_sim_periodic(t,1,runIdx), y2_sim_periodic(t,1,runIdx)];
    end

    %% Parameter Estimation using White Noise Excitation
    y1_meas_wn(:,1,runIdx) = y1_sim_wn(:,1,runIdx) + e1(:,1,runIdx);
    y2_meas_wn(:,1,runIdx) = y2_sim_wn(:,1,runIdx) + e2(:,1,runIdx);

    Phi1 = zeros(T-2, 8);
    Y1_wn = zeros(T-2,1);
    Phi2 = zeros(T-2, 8);
    Y2_wn = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_wn(t-1,1,idx), -y1_meas_wn(t-2,1,idx), -y2_meas_wn(t-1,1,idx), -y2_meas_wn(t-2,1,idx), u1_wn(t-1), u1_wn(t-2), u2_wn(t-1), u2_wn(t-2)];
        Y1_wn(idx) = y1_meas_wn(t);
        Phi2(idx,:) = [-y2_meas_wn(t-1), -y2_meas_wn(t-2), -y1_meas_wn(t-1), -y1_meas_wn(t-2), u1_wn(t-1), u1_wn(t-2), u2_wn(t-1), u2_wn(t-2)];
        Y2_wn(idx) = y2_meas_wn(t);
    end
 
    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_wn = Phi1(rows,:)\Y1_wn(rows);
        theta2_est_wn = Phi2(rows,:)\Y2_wn(rows,:);
        errorNorm_WN(runIdx, sIdx) = norm([theta1_est_wn; theta2_est_wn] - theta_true, 2);
    end
    y_meas_wn = [y1_meas_wn y2_meas_wn];

    %% Parameter Estimation using PRBS Excitation
    y1_meas_prbs(:,1,runIdx) = y1_sim_prbs(:,1,runIdx) + e1(:,1,runIdx);
    y2_meas_prbs(:,1,runIdx) = y2_sim_prbs(:,1,runIdx) + e2(:,1,runIdx);

    Phi1 = zeros(T-2, 8);
    Y1_prbs = zeros(T-2,1);
    Phi2 = zeros(T-2, 8);
    Y2_prbs = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_prbs(t-1), -y1_meas_prbs(t-2), -y2_meas_prbs(t-1), -y2_meas_prbs(t-2), u1_prbs(t-1), u1_prbs(t-2), u2_prbs(t-1), u2_prbs(t-2)];
        Y1_prbs(idx) = y1_meas_prbs(t);
        Phi2(idx,:) = [-y2_meas_prbs(t-1), -y2_meas_prbs(t-2), -y1_meas_prbs(t-1), -y1_meas_prbs(t-2), u1_prbs(t-1), u1_prbs(t-2), u2_prbs(t-1), u2_prbs(t-2)];
        Y2_prbs(idx) = y2_meas_prbs(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_prbs = Phi1(rows,:)\Y1_prbs(rows);
        theta2_est_prbs = Phi2(rows,:)\Y2_prbs(rows,:);
        errorNorm_PRBS(runIdx, sIdx) = norm([theta1_est_prbs; theta2_est_prbs] - theta_true, 2);
    end
    y_meas_prbs = [y1_meas_prbs y2_meas_prbs];

    %% Parameter Estimation using Periodic Excitation
    y1_meas_periodic(:,1,runIdx) = y1_sim_periodic(:,1,runIdx) + e1(:,1,runIdx);
    y2_meas_periodic(:,1,runIdx) = y2_sim_periodic(:,1,runIdx) + e2(:,1,runIdx);

    Phi1 = zeros(T-2,8);
    Y1_periodic = zeros(T-2,1);
    Phi2 = zeros(T-2,8);
    Y2_periodic = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_periodic(t-1), -y1_meas_periodic(t-2), -y2_meas_periodic(t-1), -y2_meas_periodic(t-2), u1_periodic(t-1), u1_periodic(t-2), u2_periodic(t-1), u2_periodic(t-2)];
        Y1_periodic(idx) = y1_meas_periodic(t);
        Phi2(idx,:) = [-y2_meas_periodic(t-1), -y2_meas_periodic(t-2), -y1_meas_periodic(t-1), -y1_meas_periodic(t-2), u1_periodic(t-1), u1_periodic(t-2), u2_periodic(t-1), u2_periodic(t-2)];
        Y2_periodic(idx) = y2_meas_periodic(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_periodic = Phi1(rows,:)\Y1_periodic(rows);
        theta2_est_periodic = Phi2(rows,:)\Y2_periodic(rows,:);
        errorNorm_periodic(runIdx, sIdx) = norm([theta1_est_periodic; theta2_est_periodic] - theta_true, 2);
    end
    y_meas_periodic = [y1_meas_periodic y2_meas_periodic];

    %% Parameter Estimation using Oracle Excitation
    y1_meas_oracle = y1_sim_oracle + e1;
    y2_meas_oracle = y2_sim_oracle + e2;

    Phi1 = zeros(T-2,8);
    Y1_oracle = zeros(T-2,1);
    Phi2 = zeros(T-2,8);
    Y2_oracle = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_oracle(t-1), -y1_meas_oracle(t-2), -y2_meas_oracle(t-1), -y2_meas_oracle(t-2), u1_oracle(t-1), u1_oracle(t-2), u2_oracle(t-1), u2_oracle(t-2)];
        Y1_oracle(idx) = y1_meas_oracle(t);
        Phi2(idx,:) = [-y2_meas_oracle(t-1), -y2_meas_oracle(t-2), -y1_meas_oracle(t-1), -y1_meas_oracle(t-2), u1_oracle(t-1), u1_oracle(t-2), u2_oracle(t-1), u2_oracle(t-2)];
        Y2_oracle(idx) = y2_meas_oracle(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_oracle = Phi1(rows,:)\Y1_oracle(rows);
        theta2_est_oracle = Phi2(rows,:)\Y2_oracle(rows,:);
        errorNorm_oracle(runIdx, sIdx) = norm([theta1_est_oracle; theta2_est_oracle] - theta_true, 2);
    end
    y_meas_oracle = [y1_meas_oracle y2_meas_oracle];
    
    %% Simulation - Active Designed Signal
    for t = 3:T
        y1_sim_oracle(t) = -A11(2)*y1_sim_oracle(t-1) - A11(3)*y1_sim_oracle(t-2) ...
            - A12(2)*y2_sim_oracle(t-1) - A12(3)*y2_sim_oracle(t-2) ...
            + B11(2)*u1_oracle(t-1) + B11(3)*u1_oracle(t-2) ...
            + B12(2)*u2_oracle(t-1) + B12(3)*u2_oracle(t-2);
        y2_sim_oracle(t) = -A22(2)*y2_sim_oracle(t-1) - A22(3)*y2_sim_oracle(t-2) ...
            - A21(2)*y1_sim_oracle(t-1) - A21(3)*y1_sim_oracle(t-2) ...
            + B21(2)*u1_oracle(t-1) + B21(3)*u1_oracle(t-2) ...
            + B22(2)*u2_oracle(t-1) + B22(3)*u2_oracle(t-2);
        y_sim_oracle(t,:) = [y1_sim_oracle(t); y2_sim_oracle(t)];
    end
end

%% Optimal Input Calculation (Time-Domain)
% % Berechnung des optimalen MIMO-Inputs mittels fmincon
% c_amp = 2*0.0005;  % Amplitudenbegrenzung (analog zu Range)
% u0 = idinput([T,du], 'rgs', Band, Range);
% u0_flat = u0(:);  % Umwandlung in einen langen Vektor, da fmincon einen Vektor erwartet
%
% options = optimoptions('fmincon',...
%     'Display','iter',...
%     'Algorithm','sqp',...
%     'MaxIterations',5000,...
%     'MaxFunctionEvaluations',30000);
%
% [u_opt_flat, fval] = fmincon(@(u_flat) cost_mimo(u_flat, T, du, c_amp), ...
%     u0_flat, [], [], [], [], -c_amp*ones(T*du,1), c_amp*ones(T*du,1), [], options);
%
% % Rückumformen in die Matrixform (T x du)
% u_opt = reshape(u_opt_flat, T, du);
%
% % Plot des optimalen Inputs für beide Kanäle
% figure;
% subplot(2,1,1);
% plot(u_opt(:,1));
% title('Optimaler Input – Kanal 1');
% xlabel('Zeit (Samples)');
% ylabel('Amplitude');
%
% subplot(2,1,2);
% plot(u_opt(:,2));
% title('Optimaler Input – Kanal 2');
% xlabel('Zeit (Samples)');
% ylabel('Amplitude');

%% Optimal Input Calculation (Frequency-Domain)


%% Simulation - Optimal Input
% y1_sim_opt = zeros(T,1);
% y2_sim_opt = zeros(T,1);
%
% for t = 3:T
%     y1_sim_opt(t) = -A11(2)*y1_sim_opt(t-1) - A11(3)*y1_sim_opt(t-2) ...
%                     - A12(2)*y2_sim_opt(t-1) - A12(3)*y2_sim_opt(t-2) ...
%                     + B11(2)*u_opt(t-1,1) + B11(3)*u_opt(t-2,1) ...
%                     + B12(2)*u_opt(t-1,2) + B12(3)*u_opt(t-2,2);
%     y2_sim_opt(t) = -A22(2)*y2_sim_opt(t-1) - A22(3)*y2_sim_opt(t-2) ...
%                     - A21(2)*y1_sim_opt(t-1) - A21(3)*y1_sim_opt(t-2) ...
%                     + B21(2)*u_opt(t-1,1) + B21(3)*u_opt(t-2,1) ...
%                     + B22(2)*u_opt(t-1,2) + B22(3)*u_opt(t-2,2);
% end
%
% figure;
% subplot(2,1,1);
% plot(y1_sim_opt);
% title('Systemausgang y1 mit optimalem Input');
% xlabel('Zeit (Samples)');
% ylabel('Amplitude');
%
% subplot(2,1,2);
% plot(y2_sim_opt);
% title('Systemausgang y2 mit optimalem Input');
% xlabel('Zeit (Samples)');
% ylabel('Amplitude');

%% Monte Carlo Simulations
M = 200; % 1000
interval = 100:200:T;
numSteps = numel(interval);
noiseStd1 = 0.000005;
noiseStd2 = 0.000005;

errorNorm_WN   = zeros(M,numSteps);
errorNorm_PRBS = zeros(M,numSteps);
errorNorm_periodic = zeros(M,numSteps);
errorNorm_oracle = zeros(M,numSteps);

for runIdx = 1:M
    % Measurement Noise
    e1 = noiseStd1 * randn(T,1);
    e2 = noiseStd2 * randn(T,1);

    %% Parameter Estimation using White Noise Excitation
    y1_meas_wn = y1_sim_wn + e1;
    y2_meas_wn = y2_sim_wn + e2;

    Phi1 = zeros(T-2, 8);
    Y1_wn = zeros(T-2,1);
    Phi2 = zeros(T-2, 8);
    Y2_wn = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_wn(t-1), -y1_meas_wn(t-2), -y2_meas_wn(t-1), -y2_meas_wn(t-2), u1_wn(t-1), u1_wn(t-2), u2_wn(t-1), u2_wn(t-2)];
        Y1_wn(idx) = y1_meas_wn(t);
        Phi2(idx,:) = [-y2_meas_wn(t-1), -y2_meas_wn(t-2), -y1_meas_wn(t-1), -y1_meas_wn(t-2), u1_wn(t-1), u1_wn(t-2), u2_wn(t-1), u2_wn(t-2)];
        Y2_wn(idx) = y2_meas_wn(t);
    end
 
    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_wn = Phi1(rows,:)\Y1_wn(rows);
        theta2_est_wn = Phi2(rows,:)\Y2_wn(rows,:);
        errorNorm_WN(runIdx, sIdx) = norm([theta1_est_wn; theta2_est_wn] - theta_true, 2);
    end
    y_meas_wn = [y1_meas_wn y2_meas_wn];

    %% Parameter Estimation using PRBS Excitation
    y1_meas_prbs = y1_sim_prbs + e1;
    y2_meas_prbs = y2_sim_prbs + e2;

    Phi1 = zeros(T-2, 8);
    Y1_prbs = zeros(T-2,1);
    Phi2 = zeros(T-2, 8);
    Y2_prbs = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_prbs(t-1), -y1_meas_prbs(t-2), -y2_meas_prbs(t-1), -y2_meas_prbs(t-2), u1_prbs(t-1), u1_prbs(t-2), u2_prbs(t-1), u2_prbs(t-2)];
        Y1_prbs(idx) = y1_meas_prbs(t);
        Phi2(idx,:) = [-y2_meas_prbs(t-1), -y2_meas_prbs(t-2), -y1_meas_prbs(t-1), -y1_meas_prbs(t-2), u1_prbs(t-1), u1_prbs(t-2), u2_prbs(t-1), u2_prbs(t-2)];
        Y2_prbs(idx) = y2_meas_prbs(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_prbs = Phi1(rows,:)\Y1_prbs(rows);
        theta2_est_prbs = Phi2(rows,:)\Y2_prbs(rows,:);
        errorNorm_PRBS(runIdx, sIdx) = norm([theta1_est_prbs; theta2_est_prbs] - theta_true, 2);
    end
    y_meas_prbs = [y1_meas_prbs y2_meas_prbs];

    %% Parameter Estimation using Periodic Excitation
    y1_meas_periodic = y1_sim_periodic + e1;
    y2_meas_periodic = y2_sim_periodic + e2;

    Phi1 = zeros(T-2,8);
    Y1_periodic = zeros(T-2,1);
    Phi2 = zeros(T-2,8);
    Y2_periodic = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_periodic(t-1), -y1_meas_periodic(t-2), -y2_meas_periodic(t-1), -y2_meas_periodic(t-2), u1_periodic(t-1), u1_periodic(t-2), u2_periodic(t-1), u2_periodic(t-2)];
        Y1_periodic(idx) = y1_meas_periodic(t);
        Phi2(idx,:) = [-y2_meas_periodic(t-1), -y2_meas_periodic(t-2), -y1_meas_periodic(t-1), -y1_meas_periodic(t-2), u1_periodic(t-1), u1_periodic(t-2), u2_periodic(t-1), u2_periodic(t-2)];
        Y2_periodic(idx) = y2_meas_periodic(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_periodic = Phi1(rows,:)\Y1_periodic(rows);
        theta2_est_periodic = Phi2(rows,:)\Y2_periodic(rows,:);
        errorNorm_periodic(runIdx, sIdx) = norm([theta1_est_periodic; theta2_est_periodic] - theta_true, 2);
    end
    y_meas_periodic = [y1_meas_periodic y2_meas_periodic];

    %% Parameter Estimation using Oracle Excitation
    y1_meas_oracle = y1_sim_oracle + e1;
    y2_meas_oracle = y2_sim_oracle + e2;

    Phi1 = zeros(T-2,8);
    Y1_oracle = zeros(T-2,1);
    Phi2 = zeros(T-2,8);
    Y2_oracle = zeros(T-2,1);

    for t = 3:T
        idx = t-2;
        Phi1(idx,:) = [-y1_meas_oracle(t-1), -y1_meas_oracle(t-2), -y2_meas_oracle(t-1), -y2_meas_oracle(t-2), u1_oracle(t-1), u1_oracle(t-2), u2_oracle(t-1), u2_oracle(t-2)];
        Y1_oracle(idx) = y1_meas_oracle(t);
        Phi2(idx,:) = [-y2_meas_oracle(t-1), -y2_meas_oracle(t-2), -y1_meas_oracle(t-1), -y1_meas_oracle(t-2), u1_oracle(t-1), u1_oracle(t-2), u2_oracle(t-1), u2_oracle(t-2)];
        Y2_oracle(idx) = y2_meas_oracle(t);
    end

    for sIdx = 1:numSteps
        i = interval(sIdx);
        rows = 1:(i-2);
        theta1_est_oracle = Phi1(rows,:)\Y1_oracle(rows);
        theta2_est_oracle = Phi2(rows,:)\Y2_oracle(rows,:);
        errorNorm_oracle(runIdx, sIdx) = norm([theta1_est_oracle; theta2_est_oracle] - theta_true, 2);
    end
    y_meas_oracle = [y1_meas_oracle y2_meas_oracle];
end

%% (Empirical) Covariance Matrix
dy = size(y_meas_wn,2);
du = size(u_wn,2);
dx = 2*(dy+du);
cov_matrix_wn = zeros(dx,dx,T);
cov_matrix_prbs = zeros(dx,dx,T);
cov_matrix_periodic = zeros(dx,dx,T);
cov_matrix_oracle = zeros(dx,dx,T);
x = [y_meas_wn(1,:), zeros(1,dy), u_wn(1,:), zeros(1,du)]';
cov_matrix_wn(:,:,1) = x*x';
x = [y_meas_prbs(1,:), zeros(1,dy), u_prbs(1,:), zeros(1,du)]';
cov_matrix_prbs(:,:,1) = x*x';
x = [y_meas_periodic(1,:), zeros(1,dy), u_periodic(1,:), zeros(1,du)]';
cov_matrix_periodic(:,:,1) = x*x';
x = [y_meas_oracle(1,:), zeros(1,dy), u_oracle(1,:), zeros(1,du)]';
cov_matrix_oracle(:,:,1) = x*x';

for i = 2:T
    yprev = y_meas_wn(i-1,:);
    uprev = u_wn(i-1,:);
    x = [y_meas_wn(i,:), yprev, u_wn(i,:), uprev]';
    cov_matrix_wn(:,:,i) = cov_matrix_wn(:,:,i-1)+x*x';
    yprev = y_meas_prbs(i-1,:);
    uprev = u_prbs(i-1,:);
    x = [y_meas_prbs(i,:), yprev, u_prbs(i,:), uprev]';
    cov_matrix_prbs(:,:,i) = cov_matrix_prbs(:,:,i-1) + x*x';
    yprev = y_meas_periodic(i-1,:);
    uprev = u_periodic(i-1,:);
    x = [y_meas_periodic(i,:), yprev, u_periodic(i,:), uprev]';
    cov_matrix_periodic(:,:,i) = cov_matrix_periodic(:,:,i-1) + x*x';
    yprev = y_meas_oracle(i-1,:);
    uprev = u_oracle(i-1,:);
    x = [y_meas_oracle(i,:), yprev, u_oracle(i,:), uprev]';
    cov_matrix_oracle(:,:,i) = cov_matrix_oracle(:,:,i-1) + x*x';
end

lambda_min_wn = zeros(T,1);
lambda_min_prbs = zeros(T,1);
lambda_min_periodic = zeros(T,1);
lambda_min_oracle = zeros(T,1);
cost_wn = zeros(T,1);
cost_prbs = zeros(T,1);
cost_periodic = zeros(T,1);
cost_oracle = zeros(T,1);

for i = 1:T
    lambda_min_wn(i) = min(eig(cov_matrix_wn(:,:,i)));
    lambda_min_prbs(i) = min(eig(cov_matrix_prbs(:,:,i)));
    lambda_min_periodic(i) = min(eig(cov_matrix_periodic(:,:,i)));
    lambda_min_oracle(i) = min(eig(cov_matrix_oracle(:,:,i)));
    % cost_wn(i) = max(eig(cov_matrix_wn(:,:,i)));
    % cost_prbs(i) = max(eig(cov_matrix_prbs(:,:,i)));
    % cost_periodic(i) = max(eig(cov_matrix_periodic(:,:,i)));
    % cost_oracle(i) = max(eig(cov_matrix_oracle(:,:,i)));
end

f4_x = figure;
f4_x.Units = 'centimeters';
f4_x.Position = [8 4 11 11/1.78];
plot(1:T,lambda_min_wn, 'LineWidth',1);
hold on
plot(1:T,lambda_min_prbs, 'LineWidth',1);
plot(1:T,lambda_min_periodic, 'LineWidth',1);
plot(1:T,lambda_min_oracle, 'LineWidth',1);
xlabel('i');
ylabel('Minimum Eigenvalue');
legend('WN','PRBS','Periodic','Oracle', 'location','northwest');

% f4_xx = figure;
% f4_xx.Units = 'centimeters';
% f4_xx.Position = [8 4 11 11/1.78];
% plot(1:T,cost_wn);
% hold on;
% plot(1:T,cost_prbs);
% plot(1:T,cost_periodic);
% xlabel('i');
% ylabel('Cost');
% legend('WN','PRBS','Periodic');

%% Error across all simulations
meanError_WN = mean(errorNorm_WN, 1);
meanError_PRBS = mean(errorNorm_PRBS, 1);
meanError_periodic = mean(errorNorm_periodic, 1);
meanError_oracle = mean(errorNorm_oracle, 1);

%% Plot – Comparison WN vs. PRBS vs. Multisine vs. Oracle
f4_3 = figure;
f4_3.Units = 'centimeters';
f4_3.Position = [8 4 11 11/1.78];
semilogy(interval,meanError_WN, 'o-', 'LineWidth', 1);
hold on;
semilogy(interval,meanError_PRBS, 's-', 'LineWidth', 1);
hold on
semilogy(interval,meanError_periodic, 'x-', 'LineWidth', 1);
hold on
semilogy(interval,meanError_oracle, '+-', 'LineWidth', 1);
ax = gca;
ax.FontSize = 8;
xlim([300 12000])
xlabel('Time $T$ [s]', 'FontSize',8);
ylabel('$\Vert\hat\theta-\theta\Vert_{2}$', 'FontSize',8);
legend('White noise', 'PRBS', 'Multisine', 'Oracle', 'Location', 'northeast', 'FontSize',8.5);
grid on;

%% Export Figure
f5 = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis','fig_4_3.pdf');
exportgraphics(f4_3,f5,'ContentType','vector');

% Power Spectral Densities (WN, PRBS, Oracle)
f4_4 = figure;
f4_4.Units = 'centimeters';
f4_4.Position = [8 4 11 11/1.78];
% PSD of Input Signals
nfft = 1024;
[pxx1,f1] = pwelch(u_wn(:,1), hann(nfft), nfft/2, nfft, 1);
[pxx2,f2] = pwelch(u_prbs(:,1), hann(nfft), nfft/2, nfft, 1);
N = 20;
pxx1_smooth = movmean(pxx1,N);
pxx2_smooth = movmean(pxx2,N);
U = fft(u_periodic);
pxx3 = abs(U).^2 / T;
f3 = (0:T-1)/T;
[pxx4,f4] = pwelch(u_oracle(:,1), hann(nfft), nfft/2, nfft, 1);
pxx4_smooth = movmean(pxx4,N);

semilogx(f1,10*log10(pxx1_smooth), 'LineWidth',1);
hold on;
semilogx(f2,10*log10(pxx2_smooth), 'LineWidth',1);
hold on
semilogx(f3(1:T/2),10*log10(pxx3(1:T/2,1)), 'LineWidth',1);
hold on
semilogx(f4,10*log10(pxx4_smooth), 'LineWidth',1);ax = gca;
ax.FontSize = 8;
xlim([0.001, 0.50])
ylim([-85 -40])
xlabel('Frequency f [Hz]', 'FontSize',8);
ylabel('PSD [dB/Hz]', 'FontSize',8);
legend('White noise', 'PRBS', 'Multisine', 'Location', 'southwest', 'FontSize',8.5);
grid on

%% Export Figure
f6 = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis','fig_4_4.pdf');
exportgraphics(f4_4,f6,'ContentType','vector');

%% Cost Function
function J = cost_mimo(u_flat, T, du, c_amp)
u = reshape(u_flat, T, du);
% Zielvarianz für ein gleichverteiltes Signal im Intervall [-c_amp, c_amp]:
target_variance = (c_amp^2)/3;
J = 0;
for k = 1:du
    var_diff = var(u(:,k)) - target_variance;
    J = J + var_diff^2;
end
end

% %% Parameter
% fs             = 1;                 % Abtastrate   [Hz]
% T              = 1000;              % Periodenlänge / #Samples
% frequencyLimits = [0.01 0.5];       % zulässiges Band [Hz]
% numComponents  = 10;                % Zahl der Einzeltöne
% nSignals       = 1;                 % ***NEU***: d   (z.B. 4 Signale)
%
% %% Multisinus erzeugen
% [y,info] = multisine_generation(frequencyLimits,fs,T,numComponents,nSignals);
%
% %% Plot
% t = (0:T-1)/fs;
% figure
% plot(t, y.')                       % '.': Zeilen als Kurven
% xlabel('Zeit [s]'), ylabel('Amplitude')
% title(sprintf('Generiertes Multisinus‑Signal (%d × T)',nSignals))
% grid on, legend("Sig "+string(1:nSignals))
%

%% Mulitsine Signal Generation
function [y,info] = multisine_generation(freqBand,fs,T,nComp,nSig)
f0 = fs/T;
kMin = max(1,ceil(freqBand(1)/f0)); % Smallest Frequency
kMax = floor(freqBand(2)/f0); % Largest Frequency
if kMax < kMin
    error('Frequency band is outside the allowable range.')
end
if nComp > (kMax-kMin+1)
    error('Too many components for the selected frequency band.');
end
info = struct('frequencies',[],'magnitude',[],'phase',[]);
info = repmat(info,nSig,1);
magTemplate   = ones(1,nComp)/nComp;
phaseTemplate = pi*(0:nComp-1).*(1:nComp)/nComp;

t = (0:T-1)/fs; % Time vector
Y = zeros(nSig, T);
for s = 1:nSig
    binList = kMin:kMax;
    p = 1./binList;       % log-uniform distribution for frequency components
    p = p/sum(p);         % Normalize
    k = datasample(binList, nComp, 'Replace', false, 'Weights', p);
    k = sort(k);
    freqs = k*f0; % Frequency [Hz]

    info(s).frequencies = freqs;
    info(s).magnitude = magTemplate;
    info(s).phase = phaseTemplate;

    % Sum of Sinusoids
    for j = 1:nComp
        Y(s,:) = Y(s,:) + magTemplate(j)*cos(2*pi*freqs(j)*t + phaseTemplate(j));
    end
    Y(s,:) = Y(s,:)/max(abs(Y(s,:)));
end
y = Y';
end

%% Optimal Input Design
function [u_opt,u_df] = design_input(A,B,k,sigma_u,seed)
if nargin<5
    seed=0;
end
rng(22+seed);
d_u = size(B,2);
omega = (2*pi/k)*(0:k-1);
u_df0 = randn(d_u,k);
u_df0(:,1) = 0;
objective = @(U) -min_eig(U,A,B,omega,k,d_u);
nonlcon = @(U)energy_con(U,sigma_u,k);
Aeq = zeros(d_u,d_u*k);
for j = 1:d_u
    Aeq(j,(j-1)*k+1) = 1;
end
beq = zeros(d_u,1);
opts = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'MaxIterations',40000,'MaxFunctionEvaluations',40000);
% warning('on','MATLAB:rankDeficientMatrix')
% dbstop if warning MATLAB:rankDeficientMatrix
u_df = fmincon(objective,u_df0,[],[],Aeq,beq,[],[],nonlcon,opts);
u_opt = zeros(k,d_u);
for j = 1:d_u
    u_t = real(ifft(u_df(j,:),k))*sqrt(k);
    u_opt(:,j) = u_t/std(u_t,1)*sqrt(sigma_u);
end
end

function lambda_min = min_eig(U,A,B,omega,k,d_u)
n = size(A,1);
U = reshape(U,d_u,k);
Gamma = zeros(n);
for j = 1:k
    z = exp(1i*omega(j));
    H = (z*eye(n)-A)\B;
    Gamma = Gamma+real(H*U(:,j)*U(:,j)'*H');
end
lambda_min = min(eig(Gamma));
end

function [c,ceq] = energy_con(U,sigma_u,k)
c = sum(U(:).^2)-k*sigma_u;
ceq = [];
end

