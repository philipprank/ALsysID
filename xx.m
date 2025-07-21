%% ----------------------------------------------------------------------
%%  Examples – Parameter Estimation for ARX Models
%% ----------------------------------------------------------------------
clear; clc; close all;
rng(20)
layout = layout_options;

%% Model Definition
% ----- SISO 1st Order ---------------------------------------------------
% x_t = [y(t-1) u(t-1)]
A1 = [0.95 -0.05; 0 0];
B1 = [0; 1];
theta_1 = [0.95; -0.05];

% ----- SISO 2nd Order ---------------------------------------------------
% x_t = [y(t-1) y(t-2) u(t-1) u(t-2)]
A2 = [1.5 -0.9  0.1 -0.2;
    1    0    0    0 ;
    0    0    0    0 ;
    0    0    1    0];
B2 = [0; 0; 1; 0];
theta_2 = [1.5 -0.9  0.1 -0.2];

% ----- MIMO 2nd Order ---------------------------------------------------
% x_t = [y1(t-1) y2(t-1) y1(t-2) y2(t-2) u1(t-1) u2(t-1) u1(t-2) u2(t-2)]
A4 = [0.7007  -0.7007   0.1363   0.4396   0.1100  -0.1200   0.1300  -0.1400;
    0.7007   0.7007   0.2413  -0.3996   0.1500  -0.1600   0.1700  -0.1800;
    1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
    0.0000   1.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
    0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000;
    0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000   0.0000;
    0.0000   0.0000   0.0000   0.0000   0.0000   1.0000   0.0000   0.0000];
B4 = [0 0; 0 0; 0 0; 0 0;
    1 0; 0 1; 0 0; 0 0];
theta_4 = [0.7007  -0.7007   0.1363   0.4396   0.1100  -0.1200   0.1300  -0.1400;
    0.7007   0.7007   0.2413  -0.3996   0.1500  -0.1600   0.1700  -0.1800];

%% Model Selection
A = A2;
B = B2;
theta_true = theta_2;
p = 2;
q = 2;

%% Dimensions
if isvector(theta_true)
    dy = 1;
    theta_true = theta_true(:).'; % SISO
else
    dy = size(theta_true,1); % MIMO
end
du = size(B,2);
dx = p*dy + q*du;
theta_true = reshape(theta_true(:).',dy,dx);
disp('-------------------------------------------------------------')
disp(['Output Dimension  dy = ' num2str(dy) ...
    ', Input Dimension du = ' num2str(du)])
disp(['Order (p,q) = (' num2str(p) ',' num2str(q) ')'])
disp(['θ-Size = ' mat2str(size(theta_true))])
disp('-------------------------------------------------------------')

%% Parameters
T = 1000; % Number of observations (1400)
M = 3; % Monte Carlo trials
k = 100; % Period length
sigma_u = 1; % Variance (Input)
sigma_e = sqrt(0.001); % Variance (Noise)
update_interval = 10;
ival = update_interval:update_interval:T;
numSt = numel(ival);

% Input (White Noise) Parameters
Range = [-1 1];
Band = [0 1];

% Preallocation
u_wn = zeros(T,du,M);
u_orac = zeros(T,du,M);
u_idalg = zeros(T,du,M);

y_meas_wn = zeros(T,dy,M);
y_meas_orac = zeros(T,dy,M);
y_meas_idalg = zeros(T,dy,M);

theta_hat_wn = zeros(dy*dx,numSt,M);
theta_hat_orac = zeros(dy*dx,numSt,M);
theta_hat_idalg = zeros(dy*dx,numSt,M);

errorNorm_wn = zeros(M,numSt);
errorNorm_orac = zeros(M,numSt);
errorNorm_idalg = zeros(M,numSt);

e = sigma_e*randn(T,dy,M);
firstIdx = max(p,q)+1;

%% Monte Carlo Simulations
updatePts = [100 200 400 800];
for sim = 1:M
    % ===== White Noise Excitation/Simulation =============================
    u_wn(:,:,sim) = sqrt(1/du)*idinput([T du],'rgs',Band,Range);
    y_meas_wn(:,:,sim) = simulate_ARX(u_wn(:,:,sim),theta_true, ...
        e(:,:,sim),p,q,firstIdx);
    [theta_hat_wn(:,:,sim),errorNorm_wn(sim,:)] = ...
        estimate_ARX(y_meas_wn(:,:,sim),u_wn(:,:,sim), ...
        theta_true,p,q,ival,firstIdx);

    % ===== Oracle Excitation/Simulation ==================================
    [u_period,U_full_loc,~] = oracle_u(A,B,k,sigma_u);
    u_orac(:,:,sim) = repmat(u_period,ceil(T/k),1);
    y_meas_orac(:,:,sim) = simulate_ARX(u_orac(:,:,sim),theta_true, ...
        e(:,:,sim),p,q,firstIdx);
    [theta_hat_orac(:,:,sim),errorNorm_orac(sim,:)] = ...
        estimate_ARX(y_meas_orac(:,:,sim),u_orac(:,:,sim), ...
        theta_true,p,q,ival,firstIdx);

    % ===== IDAlg Excitation/Simulation ===================================
    [u_idalg(:,:,sim),y_meas_idalg(:,:,sim),A_est] = ...
        idalg_u(theta_true,p,q, ...
        sigma_u,e(:,:,sim), ...
        updatePts,k,du,dy,T,Band,Range);
    [theta_hat_idalg(:,:,sim),errorNorm_idalg(sim,:)] = ...
        estimate_ARX(y_meas_idalg(:,:,sim),u_idalg(:,:,sim), ...
        theta_true,p,q,ival,firstIdx);
end

%% Theoretical Covariance Matrix
Gamma_theo = build_Gamma(U_full_loc,A,B,k);

col_wn = [0 0 0];
col_orac = [1 0 0];
col_idalg = [0 0 1];

% 10- und 90-Perzentile
p10_wn = prctile(errorNorm_wn ,10,1);
p90_wn = prctile(errorNorm_wn ,90,1);
p10_orac = prctile(errorNorm_orac,10,1);
p90_orac = prctile(errorNorm_orac,90,1);
p10_idalg = prctile(errorNorm_idalg,10,1);
p90_idalg = prctile(errorNorm_idalg,90,1);

%% Figure 4.1. - Estimation Errors & Error Bounds
f4_1 = figure;
f4_1.Units = 'centimeters';
f4_1.Position = [8 4 11 11/1.78];
hold on
x = ival;
hF1 = fill([x fliplr(x)],[p10_wn fliplr(p90_wn)],layout.colors(1,:),'EdgeColor','none','FaceAlpha',0.25);
hF2 = fill([x fliplr(x)],[p10_orac fliplr(p90_orac)],layout.colors(2,:),'EdgeColor','none','FaceAlpha',0.25);
hF3 = fill([x fliplr(x)],[p10_idalg fliplr(p90_idalg)],layout.colors(3,:),'EdgeColor','none','FaceAlpha',0.25);

hL1 = plot(x,mean(errorNorm_wn,1),'LineWidth',1.2);
hL2 = plot(x,mean(errorNorm_orac,1),'LineWidth',1.2);
hL3 = semilogy(x,mean(errorNorm_idalg,1),'LineWidth',1.2);
ylim([0 0.05])
hold off;
ax.FontSize = 9;
xlabel('$T$','FontSize',9);
ylabel('$\Vert\hat\theta - \theta\Vert_{2}$','FontSize',9);
legend([hL1 hL2 hL3], ...
    'White Noise', ...
    'Oracle', ...
    'IDAlg', ...
    'Location','best');
grid on

%% Export Figure 4.1 - PDF & TikZ
pdfFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_4_01.pdf');
exportgraphics(f4_1, pdfFile, 'ContentType','vector');
% Export with matlab2tikz:
tikzFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_4_01.tex');
matlab2tikz(tikzFile, 'figurehandle',f4_1,'standalone', true, 'extraPreamble', '\usepackage{lmodern}');

%% Empirical Covariance Matrix
dx_cov = p*dy + q*du;
lambda_min = struct('wn',[],'orac',[],'idalg',[]);
lambda_max = struct('wn',[],'orac',[],'idalg',[]);

lambda_min.wn = zeros(T,M);
lambda_max.wn = zeros(T,M);
lambda_min.orac = zeros(T,M);
lambda_max.orac = zeros(T,M);
lambda_min.idalg = zeros(T,M);
lambda_max.idalg = zeros(T,M);

Gamma_wn_all = zeros(dx_cov,dx_cov,M);
Gamma_orac_all = zeros(dx_cov,dx_cov,M);
Gamma_idalg_all = zeros(dx_cov,dx_cov,M);

for m = 1:M
    G_wn = zeros(dx_cov);
    G_orac = zeros(dx_cov);
    G_idalg = zeros(dx_cov);
    for t = firstIdx:T
        % White Noise
        x = [reshape(flipud(y_meas_wn(t-p:t-1,:,m)).',[],1) ;
            reshape(flipud(u_wn(t-q:t-1,:,m)).',[],1)];
        G_wn = G_wn + x*x';
        lambda_min.wn(t,m) = min(eig(G_wn));
        lambda_max.wn(t,m) = max(eig(G_wn));

        % Oracle
        x = [reshape(flipud(y_meas_orac(t-p:t-1,:,m)).',[],1) ;
            reshape(flipud(u_orac(t-q:t-1,:,m)).',[],1)];
        G_orac = G_orac + x*x';
        lambda_min.orac(t,m) = min(eig(G_orac));
        lambda_max.orac(t,m) = max(eig(G_orac));

        % IDAlg
        x = [reshape(flipud(y_meas_idalg(t-p:t-1,:,m)).',[],1) ;
            reshape(flipud(u_idalg(t-q:t-1,:,m)).',[],1)];
        G_idalg = G_idalg + x*x';
        lambda_min.idalg(t,m) = min(eig(G_idalg));
        lambda_max.adapt(t,m) = max(eig(G_idalg));
    end
    Gamma_wn_all (:,:,m) = G_wn;
    Gamma_orac_all(:,:,m) = G_orac;
    Gamma_idalg_all(:,:,m) = G_idalg;
end

f4_2 = figure;
f4_2.Units = 'centimeters';
f4_2.Position = [8 4 11 11/1.78];
plot_eig(lambda_min,lambda_max,T,firstIdx)

%% Export Figure 4.1 - PDF & TikZ
pdfFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_4_02.pdf');
exportgraphics(f4_2, pdfFile, 'ContentType','vector');
% Export with matlab2tikz:
tikzFile = fullfile('C:\Users\phili\OneDrive - bwedu\ALsysID\06_thesis\plots','fig_4_02.tex');
matlab2tikz(tikzFile, 'figurehandle',f4_2,'standalone', true, 'extraPreamble', '\usepackage{lmodern}');

%% Asymptotic Distribution
errorNorm_vec_wn = zeros(M,numel(theta_true));
errorNorm_vec_orac = zeros(M,numel(theta_true));
errorNorm_vec_idalg = zeros(M,numel(theta_true));

for sim = 1:M
    diff_wn = theta_hat_wn(:,end,sim) - theta_true(:);
    diff_orac = theta_hat_orac(:,end,sim) - theta_true(:);
    diff_idalg = theta_hat_idalg(:,end,sim) - theta_true(:);
    errorNorm_vec_wn(sim,:) = sqrt(T)*diff_wn;
    errorNorm_vec_orac(sim,:) = sqrt(T)*diff_orac;
    errorNorm_vec_idalg(sim,:) = sqrt(T)*diff_idalg;
end

figure; hold on
histogram(errorNorm_vec_wn(:,1),40,'Normalization','pdf','FaceColor',layout.colors(1,:));
histogram(errorNorm_vec_orac(:,1),40,'Normalization','pdf','FaceColor',layout.colors(2,:));
histogram(errorNorm_vec_idalg(:,1),40,'Normalization','pdf','FaceColor',layout.colors(3,:));
legend('White-noise','Oracle','IDAlg')

%% Disrepancy between Γ_emp vs Γ_theo
Gamma_emp_opt = mean(Gamma_orac_all,3)/(T-firstIdx+1);  % Oracle
Gamma_emp_idalg  = mean(Gamma_idalg_all,3)/(T-firstIdx+1);  % IDAlg

err_F_opt = norm(Gamma_emp_opt-Gamma_theo,'fro');
rel_F_opt = err_F_opt/norm(Gamma_theo,'fro');
err_F_idalg  = norm(Gamma_emp_idalg -Gamma_theo,'fro');
rel_F_ad  = err_F_idalg/norm(Gamma_theo,'fro');
fprintf('\nOracle  : Frobenius-Fehler %.3e  (%.1f %%)\n', err_F_opt,100*rel_F_opt);
fprintf('Adaptive: Frobenius-Fehler %.3e  (%.1f %%)\n', err_F_idalg ,100*rel_F_ad );


%% Auxiliary Functions
function y = simulate_ARX(u,theta,e,p,q,firstIdx)
[T,du] = size(u);
dy = size(e,2);
y = zeros(T,dy);
for t = firstIdx:T
    y_lags = reshape(flipud(y(t-p:t-1,:)).',[],1);
    u_lags = reshape(flipud(u(t-q:t-1,:)).',[],1);
    phi = [y_lags; u_lags];
    y(t,:) = (theta*phi).' + e(t,:);
end
end

function [theta_hat,err] = estimate_ARX(y,u,theta_true,p,q,ival,firstIdx)
[T,dy] = size(y);
du = size(u,2);
dx = p*dy + q*du;
Ntot = T-firstIdx+1;
Phi = zeros(Ntot,dx);
Y = zeros(Ntot,dy);
for t = firstIdx:T
    j = t-firstIdx+1;
    y_lags = reshape(flipud(y(t-p:t-1,:)).',[],1);
    u_lags = reshape(flipud(u(t-q:t-1,:)).',[],1);
    Phi(j,:) = [y_lags; u_lags].';
    Y(j,:) = y(t,:);
end
numSt = numel(ival);
theta_hat = zeros(dy*dx,numSt);
err = zeros(1,numSt);
for k = 1:numSt
    N = ival(k) - firstIdx+1;
    theta_est = (Phi(1:N,:)\Y(1:N,:)).';
    theta_hat(:,k) = theta_est(:);
    err(k) = norm(theta_est(:) - theta_true(:));
end
end

function plot_eig(lambda_min,lambda_max,T,firstIdx)
lm = @(S) mean(S,2);
plot(firstIdx:T,lm(lambda_min.wn(firstIdx:T,:)),'LineWidth',1.2);
hold on
plot(firstIdx:T,lm(lambda_min.orac(firstIdx:T,:)),'LineWidth',1.2)
plot(firstIdx:T,lm(lambda_min.idalg(firstIdx:T,:)),'LineWidth',1.2)
xlabel('$T$');
ylabel('$\lambda_{\min}\left(\sum_{t=1}^{T}x_{t} x_{t}^{\mathsf{T}}\right)$');
grid on
legend('White-noise','Oracle','IDAlg','Location','northwest')
end

function [u_opt,U_full,U_pos] = oracle_u(A,B,k,sigma_u)
d_u = size(B,2);
m = k/2-1;
x0 = randn(2*d_u*m,1)/sqrt(2);
objective = @(x) -min_eig(assemble_U(x,d_u,k),A,B,k);
nonlcon = @(x) energy_con(assemble_U(x,d_u,k),sigma_u,k);
opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter', ...
    'MaxIterations',4000,'MaxFunctionEvaluations',4000);
U_opt = fmincon(objective,x0,[],[],[],[],[],[],nonlcon,opts);
U_full = assemble_U(U_opt,d_u,k);
U_pos = U_full(:,2:m+1);
u_opt = ifft(U_full,[],2,'symmetric').';
end

function lam = min_eig(U,A,B,k)
n = size(A,1); omega = (2*pi/k)*(0:k-1);
Gamma = zeros(n);
for j = 1:k
    z = exp(1i*omega(j)); H = (z*eye(n)-A)\B;
    Gamma = Gamma + real(H*U(:,j)*U(:,j)'*H');
end
lam = min(eig(Gamma));
end

function U = assemble_U(x,d_u,k)
m = k/2-1;
x = reshape(x,2*d_u,m);
Upos = x(1:d_u,:) + 1i*x(d_u+1:end,:);
U = zeros(d_u,k);
U(:,2:m+1) = Upos;
U(:,k/2+2:end) = conj(fliplr(Upos));
end

function [c,ceq] = energy_con(U,sigma_u,k)
c = sum(abs(U(:)).^2) - sigma_u*k^2;
ceq = [];
end

function Gamma = build_Gamma(U_df,A,B,k)
d_u = size(B,2);
omega = (2*pi/k)*(0:k-1);
U = reshape(U_df,d_u,k);
n = size(A,1);
Gamma = zeros(n);
for j = 1:k
    z = exp(1i*omega(j));
    H = (z*eye(n)-A)\B;
    Gamma = Gamma + real(H*U(:,j)*U(:,j)'*H');
end
Gamma = Gamma/(k^2);
end

function [u_idalg,y,A_est] = idalg_u(theta_true,p,q, ...
    sigma_u,e,updatePts,k,du,dy,T,Band,Range)
firstIdx = max(p,q)+1;
u_idalg = zeros(T,du);
y = zeros(T,dy);
% White Noise as
u_idalg(1:updatePts(1),:) = sqrt(1/du)* ...
    idinput([updatePts(1) du],'rgs',Band,Range);

currentPeriod = [];
idxPeriod = 1;
nextUpdIdx = 1;

for t = 1:T
    if nextUpdIdx <= numel(updatePts) && t == updatePts(nextUpdIdx)+1
        % Parameter Estimation
        theta_est = estimate_theta_online(y(1:t-1,:),u_idalg(1:t-1,:),p,q,firstIdx,dy,du);
        [A_est,B_est] = build_state_matrices(theta_est,p,q,dy,du);
        [u_period,~,~] = oracle_u(A_est,B_est,k,sigma_u);
        currentPeriod = u_period;
        idxPeriod = 1;
        nextUpdIdx = nextUpdIdx + 1;
    end

    if t > updatePts(1)
        u_idalg(t,:) = currentPeriod(idxPeriod,:);
        idxPeriod = idxPeriod + 1;
        if idxPeriod > k, idxPeriod = 1; end
    end

    if t >= firstIdx
        y_lags = reshape(flipud(y(t-p:t-1,:)).',[],1);
        u_lags = reshape(flipud(u_idalg(t-q:t-1,:)).',[],1);
        phi = [y_lags; u_lags];
        y(t,:) = (theta_true*phi).' + e(t,:);
    end
end
end

function theta_est = estimate_theta_online(y,u,p,q,firstIdx,dy,du)
Tcurr = size(y,1);
dx = p*dy + q*du;
N = Tcurr-firstIdx+1;
Phi = zeros(N,dx);
Y = zeros(N,dy);
for t = firstIdx:Tcurr
    j = t-firstIdx+1;
    y_lags = reshape(flipud(y(t-p:t-1,:)).',[],1);
    u_lags = reshape(flipud(u(t-q:t-1,:)).',[],1);
    Phi(j,:) = [y_lags; u_lags].';
    Y(j,:) = y(t,:);
end
theta_est = (Phi\Y).';
end

function [A_est,B_est] = build_state_matrices(theta_est,p,q,dy,du)
dx = p*dy + q*du;
A_est = zeros(dx);
A_est(1:dy,:) = reshape(theta_est,dy,dx);
for i = 1:dy*(p-1)
    A_est(dy+i,i) = 1;
end
row_off = p*dy;
for j = 1:du*(q-1)
    A_est(row_off+du+j,row_off+j) = 1;
end
B_est = zeros(dx,du);
B_est(row_off+1:row_off+du,:) = eye(du);
end
