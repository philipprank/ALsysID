%%  Quadruple-Tank Process – Monte-Carlo-Analyse
clear; clc; close all;
rng(32);

%% True ARX Model Parameters
A11 = [1, -1.9659,  0.9663];
A12 = [0, 0.0471, -0.0455];
B11 = [0, 42.1906, -40.0763];
B12 = [0, 89.4773, -88.2050];

A21 = [0, 0.0058, -0.0068];
A22 = [1, -1.9228, 0.9250];
B21 = [0, 73.5355, -71.5468];
B22 = [0, 79.0570, -78.9397];

theta1_true = [A11(2:3), A12(2:3), B11(2:3), B12(2:3)].';
theta2_true = [A22(2:3), A21(2:3), B21(2:3), B22(2:3)].';
theta_true = [theta1_true; theta2_true];
p = numel(theta_true);

%% Simulation Parameters
T = 10000; % Number of observations
M = 200; % Number of simulations
steps = 200;
ival = 100:steps:T;
numSteps = numel(ival);

noiseStd1 = 5e-6;                % Variance u_{1}
noiseStd2 = 5e-6;                % Variance u_{1}
e1 = noiseStd1*randn(T,M);       % Measurement noise y_{1}
e2 = noiseStd2*randn(T,M);       % Measurement noise y_{2}

% Parameters White Noise & PRBS
du = 2; % Input dimension
Range = [-5e-4, 5e-4];
Band = [0.002, 1];
% Parameters Multisine
numComp = 10; % Number of frequency components for periodic signal
freq_limit = [0.001 1]; % Frequency limits for periodic signal
freq = zeros(M,numComp);
% Parameters Oracle
A_arx = [-A11(2:end) -A12(2:end); -A21(2:end) -A22(2:end)];
B_arx = [B11(2:end) B12(2:end); B21(2:end) B22(2:end)];
A = [A_arx B_arx; eye(2) zeros(2,6); zeros(2,8); zeros(2,4) eye(2) zeros(2)];
B = [zeros(4,2); eye(2); zeros(2)];
sigma_u = 2.5e-07;
k = 1000;
[u_opt,u_df] = design_input(A,B,k,sigma_u);
u_oracle = repmat(u_opt,ceil(T/k),1);
u_oracle = u_oracle(1:T,:);
u1_oracle = u_oracle(:,1);
u2_oracle = u_oracle(:,2);

%% Pre-Allocation (White Noise)
u_wn = zeros(T,du,M);
y1_sim_wn = zeros(T,M);
y2_sim_wn = zeros(T,M);
y1_meas_wn = zeros(T,M);
y2_meas_wn = zeros(T,M);
theta1_hat_WN = zeros(8,numSteps,M);
theta2_hat_WN = zeros(8,numSteps,M);
errorNorm_WN = zeros(M,numSteps);

%% Pre-Allocation (PRBS)
u_prbs = zeros(T,du,M);
y1_sim_prbs = zeros(T,M);
y2_sim_prbs = zeros(T,M);
y1_meas_prbs = zeros(T,M);
y2_meas_prbs = zeros(T,M);
theta1_hat_PRBS = zeros(8,numSteps,M);
theta2_hat_PRBS = zeros(8,numSteps,M);
errorNorm_PRBS = zeros(M,numSteps);

%% Pre-Allocation (Mulitsine)
u_sine = zeros(T,du,M);
y1_sim_sine = zeros(T,M);
y2_sim_sine = zeros(T,M);
y1_meas_sine = zeros(T,M);
y2_meas_sine = zeros(T,M);
theta1_hat_sine = zeros(8,numSteps,M);
theta2_hat_sine = zeros(8,numSteps,M);
errorNorm_sine = zeros(M,numSteps);

%% Pre-Allocation (Oralce)
y1_sim_oracle = zeros(T,M);
y2_sim_oracle = zeros(T,M);
y1_meas_oracle = zeros(T,M);
y2_meas_oracle = zeros(T,M);
theta1_hat_oracle = zeros(8,numSteps,M);
theta2_hat_oracle = zeros(8,numSteps,M);
errorNorm_oracle = zeros(M,numSteps);

%% Monte-Carlo-Simulations
for sim = 1:M
    %% White Noise Excitation
    u_wn(:,:,sim) = idinput([T,du],'rgs',Band,Range);
    u1 = u_wn(:,1,sim);
    u2 = u_wn(:,2,sim);
    for t = 3:T
        y1_sim_wn(t,sim) = -A11(2)*y1_sim_wn(t-1,sim) ...
            - A11(3)*y1_sim_wn(t-2,sim) ...
            - A12(2)*y2_sim_wn(t-1,sim) ...
            - A12(3)*y2_sim_wn(t-2,sim) ...
            + B11(2)*u1(t-1) + B11(3)*u1(t-2) ...
            + B12(2)*u2(t-1) + B12(3)*u2(t-2) ...
            + e1(t,sim);

        y2_sim_wn(t,sim) = -A22(2)*y2_sim_wn(t-1,sim) ...
            - A22(3)*y2_sim_wn(t-2,sim) ...
            - A21(2)*y1_sim_wn(t-1,sim) ...
            - A21(3)*y1_sim_wn(t-2,sim) ...
            + B21(2)*u1(t-1) + B21(3)*u1(t-2) ...
            + B22(2)*u2(t-1) + B22(3)*u2(t-2) ...
            + e2(t,sim);
    end

    % Add Measurenment Noise
    y1_meas_wn(:,sim) = y1_sim_wn(:,sim);
    y2_meas_wn(:,sim) = y2_sim_wn(:,sim);

    Phi1 = zeros(T-2,8);
    Phi2 = zeros(T-2,8);
    Y1   = zeros(T-2,1);
    Y2   = zeros(T-2,1);

    for t = 3:T
        k = t-2;
        Phi1(k,:) = [-y1_meas_wn(t-1,sim), -y1_meas_wn(t-2,sim), ...
            -y2_meas_wn(t-1,sim), -y2_meas_wn(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y1(k) = y1_meas_wn(t,sim);

        Phi2(k,:) = [-y2_meas_wn(t-1,sim), -y2_meas_wn(t-2,sim), ...
            -y1_meas_wn(t-1,sim), -y1_meas_wn(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y2(k) = y2_meas_wn(t,sim);
    end

    % Least Squares Estimation
    for sIdx = 1:numSteps
        N = ival(sIdx) - 2;
        theta1_est = Phi1(1:N,:)\Y1(1:N);
        theta2_est = Phi2(1:N,:)\Y2(1:N);

        theta1_hat_WN(:,sIdx,sim) = theta1_est;
        theta2_hat_WN(:,sIdx,sim) = theta2_est;

        errorNorm_WN(sim,sIdx) = norm([theta1_est;theta2_est] - theta_true,2);
    end

    %% PRBS Excitation
    u_prbs(:,:,sim) = idinput([T,du],'prbs',Band,Range);
    u1 = u_prbs(:,1,sim);
    u2 = u_prbs(:,2,sim);

    for t = 3:T
        y1_sim_prbs(t,sim) = -A11(2)*y1_sim_prbs(t-1,sim) ...
            - A11(3)*y1_sim_prbs(t-2,sim) ...
            - A12(2)*y2_sim_prbs(t-1,sim) ...
            - A12(3)*y2_sim_prbs(t-2,sim) ...
            + B11(2)*u1(t-1) + B11(3)*u1(t-2) ...
            + B12(2)*u2(t-1) + B12(3)*u2(t-2) ...
            + e1(t,sim);

        y2_sim_prbs(t,sim) = -A22(2)*y2_sim_prbs(t-1,sim) ...
            - A22(3)*y2_sim_prbs(t-2,sim) ...
            - A21(2)*y1_sim_prbs(t-1,sim) ...
            - A21(3)*y1_sim_prbs(t-2,sim) ...
            + B21(2)*u1(t-1) + B21(3)*u1(t-2) ...
            + B22(2)*u2(t-1) + B22(3)*u2(t-2) ...
            + e2(t,sim);
    end

    y1_meas_prbs(:,sim) = y1_sim_prbs(:,sim);
    y2_meas_prbs(:,sim) = y2_sim_prbs(:,sim);

    Phi1 = zeros(T-2,8);
    Phi2 = zeros(T-2,8);
    Y1   = zeros(T-2,1);
    Y2   = zeros(T-2,1);
    for t = 3:T
        k = t-2;
        Phi1(k,:) = [-y1_meas_prbs(t-1,sim), -y1_meas_prbs(t-2,sim), ...
            -y2_meas_prbs(t-1,sim), -y2_meas_prbs(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y1(k) = y1_meas_prbs(t,sim);

        Phi2(k,:) = [-y2_meas_prbs(t-1,sim), -y2_meas_prbs(t-2,sim), ...
            -y1_meas_prbs(t-1,sim), -y1_meas_prbs(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y2(k) = y2_meas_prbs(t,sim);
    end

    for sIdx = 1:numSteps
        N = ival(sIdx) - 2;
        theta1_est = Phi1(1:N,:)\Y1(1:N);
        theta2_est = Phi2(1:N,:)\Y2(1:N);

        theta1_hat_PRBS(:,sIdx,sim) = theta1_est;
        theta2_hat_PRBS(:,sIdx,sim) = theta2_est;

        errorNorm_PRBS(sim,sIdx) = norm([theta1_est;theta2_est]-theta_true,2);
    end

    %% Multisine Excitation
    [u, info] = multisine_generation(freq_limit,1,T,numComp,du);
    freq(sim,:) = info(1).frequencies;
    currVar = var(u);
    scale = sqrt(2.5*10^(-7)./currVar);
    u_sine(:,:,sim) = scale.*u;
    u1 = scale(1)*u(:,1);
    u2 = scale(2)*u(:,2);

    for t = 3:T
        y1_sim_sine(t,sim) = -A11(2)*y1_sim_sine(t-1,sim) ...
            - A11(3)*y1_sim_sine(t-2,sim) ...
            - A12(2)*y2_sim_sine(t-1,sim) ...
            - A12(3)*y2_sim_sine(t-2,sim) ...
            + B11(2)*u1(t-1) + B11(3)*u1(t-2) ...
            + B12(2)*u2(t-1) + B12(3)*u2(t-2) ...
            + e1(t,sim);

        y2_sim_sine(t,sim) = -A22(2)*y2_sim_sine(t-1,sim) ...
            - A22(3)*y2_sim_sine(t-2,sim) ...
            - A21(2)*y1_sim_sine(t-1,sim) ...
            - A21(3)*y1_sim_sine(t-2,sim) ...
            + B21(2)*u1(t-1) + B21(3)*u1(t-2) ...
            + B22(2)*u2(t-1) + B22(3)*u2(t-2) ...
            + e2(t,sim);
    end

    y1_meas_sine(:,sim) = y1_sim_sine(:,sim);
    y2_meas_sine(:,sim) = y2_sim_sine(:,sim);

    Phi1 = zeros(T-2,8);
    Phi2 = zeros(T-2,8);
    Y1   = zeros(T-2,1);
    Y2   = zeros(T-2,1);
    for t = 3:T
        k = t-2;
        Phi1(k,:) = [-y1_meas_sine(t-1,sim), -y1_meas_sine(t-2,sim), ...
            -y2_meas_sine(t-1,sim), -y2_meas_sine(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y1(k) = y1_meas_sine(t,sim);

        Phi2(k,:) = [-y2_meas_sine(t-1,sim), -y2_meas_sine(t-2,sim), ...
            -y1_meas_sine(t-1,sim), -y1_meas_sine(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y2(k) = y2_meas_sine(t,sim);
    end

    for sIdx = 1:numSteps
        N = ival(sIdx) - 2;
        theta1_est = Phi1(1:N,:)\Y1(1:N);
        theta2_est = Phi2(1:N,:)\Y2(1:N);

        theta1_hat_sine(:,sIdx,sim) = theta1_est;
        theta2_hat_sine(:,sIdx,sim) = theta2_est;

        errorNorm_sine(sim,sIdx) = norm([theta1_est;theta2_est]-theta_true,2);
    end

    %% Input Design Exciatation
    u1 = u_oracle(:,1);
    u2 = u_oracle(:,2);
    for t = 3:T
        y1_sim_oracle(t,sim) = -A11(2)*y1_sim_oracle(t-1,sim) ...
            - A11(3)*y1_sim_oracle(t-2,sim) ...
            - A12(2)*y2_sim_oracle(t-1,sim) ...
            - A12(3)*y2_sim_oracle(t-2,sim) ...
            + B11(2)*u1(t-1) + B11(3)*u1(t-2) ...
            + B12(2)*u2(t-1) + B12(3)*u2(t-2) ...
            + e1(t,sim);

        y2_sim_oracle(t,sim) = -A22(2)*y2_sim_oracle(t-1,sim) ...
            - A22(3)*y2_sim_oracle(t-2,sim) ...
            - A21(2)*y1_sim_oracle(t-1,sim) ...
            - A21(3)*y1_sim_oracle(t-2,sim) ...
            + B21(2)*u1(t-1) + B21(3)*u1(t-2) ...
            + B22(2)*u2(t-1) + B22(3)*u2(t-2) ...
            + e2(t,sim);
    end

    y1_meas_oracle(:,sim) = y1_sim_oracle(:,sim);
    y2_meas_oracle(:,sim) = y2_sim_oracle(:,sim);

    Phi1 = zeros(T-2,8);
    Phi2 = zeros(T-2,8);
    Y1   = zeros(T-2,1);
    Y2   = zeros(T-2,1);
    for t = 3:T
        k = t-2;
        Phi1(k,:) = [-y1_meas_oracle(t-1,sim), -y1_meas_oracle(t-2,sim), ...
            -y2_meas_oracle(t-1,sim), -y2_meas_oracle(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y1(k) = y1_meas_oracle(t,sim);

        Phi2(k,:) = [-y2_meas_oracle(t-1,sim), -y2_meas_oracle(t-2,sim), ...
            -y1_meas_oracle(t-1,sim), -y1_meas_oracle(t-2,sim), ...
            u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
        Y2(k) = y2_meas_oracle(t,sim);
    end

    for sIdx = 1:numSteps
        N = ival(sIdx) - 2;
        theta1_est = Phi1(1:N,:)\Y1(1:N);
        theta2_est = Phi2(1:N,:)\Y2(1:N);

        theta1_hat_oracle(:,sIdx,sim) = theta1_est;
        theta2_hat_oracle(:,sIdx,sim) = theta2_est;

        errorNorm_oracle(sim,sIdx) = norm([theta1_est;theta2_est]-theta_true,2);
    end
end

%% Plot - Estimation Error
meanError_wn   = mean(errorNorm_WN ,1);
meanError_prbs = mean(errorNorm_PRBS,1);
meanError_sine = mean(errorNorm_sine,1);
meanError_oracle = mean(errorNorm_oracle,1);

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [8 4 11 11/1.78];
semilogy(ival,meanError_wn ,'o-','LineWidth',1);
hold on;
semilogy(ival,meanError_prbs,'s-','LineWidth',1);
semilogy(ival,meanError_sine,'x-','LineWidth',1);
semilogy(ival,meanError_oracle,'o-','LineWidth',1);
hold off;
ax.FontSize = 8;
xlabel('Time T','FontSize',8);
ylabel('$\Vert\hat\theta - \theta\Vert_{2}$','FontSize',8);
legend('White Noise','PRBS','Multisine','Oralce', 'Location','northEast');
grid on;

%% Plot - Power Spectral Density
fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position = [8 4 11 11/1.78];
% PSD of Input Signals
nfft = 1024;
[pxx1,f1] = pwelch(u_wn(:,1), hann(nfft), nfft/2, nfft, 1);
[pxx2,f2] = pwelch(u_prbs(:,1), hann(nfft), nfft/2, nfft, 1);
[pxx3,f3] = pwelch(u_oracle(:,1), hann(nfft), nfft/2, nfft, 1);
N = 20;
pxx1_smooth = movmean(pxx1,N);
pxx2_smooth = movmean(pxx2,N);
pxx3_smooth = movmean(pxx3,N);
U = fft(u_sine);
pxx4 = abs(U).^2/T;
f4 = (0:T-1)/T;

semilogx(f1,10*log10(pxx1_smooth), 'LineWidth',1);
hold on;
semilogx(f2,10*log10(pxx2_smooth), 'LineWidth',1);
semilogx(f4(1:T/2),10*log10(pxx4(1:T/2,1)), 'LineWidth',1);
semilogx(f3,10*log10(pxx3_smooth), 'LineWidth',1);
hold off;
ax.FontSize = 8;
xlim([0.001, 0.50])
ylim([-85 -40])
xlabel('Frequency f [Hz]', 'FontSize',8);
ylabel('PSD [dB/Hz]', 'FontSize',8);
legend('White Noise', 'PRBS','Multisine','Oracle', 'Location','southwest', 'FontSize',8.5);

%% Plot – Empirical Covariance Matrix
dy = 2;  du = 2;
dx = 2*(dy+du);

lambda_min_wn   = zeros(T,M);
lambda_min_prbs = zeros(T,M);
lambda_min_sine = zeros(T,M);
lambda_min_oracle = zeros(T,M);

for m = 1:M
    Gamma_wn   = zeros(dx,dx);
    Gamma_prbs = zeros(dx,dx);
    Gamma_sine = zeros(dx,dx);
    Gamma_oracle = zeros(dx,dx);

    for t = 1:T
        % White-Noise Signal
        y_now = [y1_meas_wn(t,m), y2_meas_wn(t,m)];
        u_now = [u_wn(t,1,m), u_wn(t,2,m)];
        if t==1
            y_prev = [0 0];
            u_prev = [0 0];
        else
            y_prev = [y1_meas_wn(t-1,m), y2_meas_wn(t-1,m)];
            u_prev = [u_wn(t-1,1,m), u_wn(t-1,2,m)];
        end

        phi = [y_now, y_prev, u_now, u_prev].';
        Gamma_wn = Gamma_wn + phi*phi.';
        lambda_min_wn(t,m) = min(eig(Gamma_wn));

        % PRBS Signal
        y_now = [y1_meas_prbs(t,m), y2_meas_prbs(t,m)];
        u_now = [u_prbs(t,1,m), u_prbs(t,2,m)];
        if t==1
            y_prev = [0 0];
            u_prev = [0 0];
        else
            y_prev = [y1_meas_prbs(t-1,m), y2_meas_prbs(t-1,m)];
            u_prev = [u_prbs(t-1,1,m), u_prbs(t-1,2,m)];
        end

        phi = [y_now, y_prev, u_now, u_prev].';
        Gamma_prbs = Gamma_prbs + phi*phi.';
        lambda_min_prbs(t,m) = min(eig(Gamma_prbs));

        % Multisine Signal
        y_now = [y1_meas_sine(t,m) , y2_meas_sine(t,m)];
        u_now = [u_sine(t,1,m) , u_sine(t,2,m)];
        if t==1
            y_prev = [0 0];
            u_prev = [0 0];
        else
            y_prev = [y1_meas_sine(t-1,m) , y2_meas_sine(t-1,m)];
            u_prev = [u_sine(t-1,1,m) , u_sine(t-1,2,m)];
        end

        phi = [y_now, y_prev, u_now, u_prev].';
        Gamma_sine = Gamma_sine + phi*phi.';
        lambda_min_sine(t,m) = min(eig(Gamma_sine));

        % Oralce Signal
        y_now = [y1_meas_oracle(t,m) , y2_meas_oracle(t,m)];
        u_now = [u_oracle(t,1) , u_oracle(t,2)];
        if t==1
            y_prev = [0 0];
            u_prev = [0 0];
        else
            y_prev = [y1_meas_oracle(t-1,m) , y2_meas_oracle(t-1,m)];
            u_prev = [u_oracle(t-1,1) , u_oracle(t-1,2)];
        end
        
        phi = [y_now, y_prev, u_now, u_prev].';
        Gamma_oracle = Gamma_oracle + phi*phi.';
        lambda_min_oracle(t,m) = min(eig(Gamma_oracle));
    end
end

% Mean over All Monte Carlo Runs
lambda_min_wn_mean     = mean(lambda_min_wn, 2);
lambda_min_prbs_mean   = mean(lambda_min_prbs,2);
lambda_min_sine_mean   = mean(lambda_min_sine,2);
lambda_min_oracle_mean = mean(lambda_min_oracle,2);

f3 = figure;
f3.Units = 'centimeters';
f3.Position = [8 4 11 11/1.78];

plot(1:T,lambda_min_wn_mean, 'LineWidth',1);
hold on;
plot(1:T,lambda_min_prbs_mean, 'LineWidth',1);
plot(1:T,lambda_min_sine_mean, 'LineWidth',1);
plot(1:T,lambda_min_oracle_mean, 'LineWidth',1);
hold off;

xlabel('T');
ylabel('Minimum Eigenvalue');
legend('WN','PRBS','Periodic','Oracle', 'Location','northwest');
grid on;

%% (Optional) Δθ_T = √T(θ̂_T – θ)
scaledErr_WN   = zeros(M,p);
scaledErr_PRBS = zeros(M,p);
for sim = 1:M
    theta_hat_T_WN   = [theta1_hat_WN(:,numSteps,sim);  ...
        theta2_hat_WN(:,numSteps,sim)];
    theta_hat_T_PRBS = [theta1_hat_PRBS(:,numSteps,sim);...
        theta2_hat_PRBS(:,numSteps,sim)];

    scaledErr_WN(sim,:)   = sqrt(T)*(theta_hat_T_WN  - theta_true).';
    scaledErr_PRBS(sim,:) = sqrt(T)*(theta_hat_T_PRBS- theta_true).';
end

%% 8.  Statistiken:  Mean Values & Covariance Matrix
mu_WN   = mean(scaledErr_WN,1);
mu_PRBS = mean(scaledErr_PRBS,1);

Sigma_WN   = cov(scaledErr_WN);
Sigma_PRBS = cov(scaledErr_PRBS);

fprintf('\nMax. |Mittelwert|  (WN)  : %.3g\n',max(abs(mu_WN)));
fprintf('Max. |Mittelwert|  (PRBS): %.3g\n',max(abs(mu_PRBS)));

%% (optional) Theoretische Kovarianz ≈ σ²·Γ_∞^{-1}
% Hier exemplarisch Gamma_hat aus dem ersten White-Noise-Versuch:
% sim = 1;
% u1 = u_wn(:,1,sim); u2 = u_wn(:,2,sim);
% Phi1 = zeros(T-2,8);  Phi2 = zeros(T-2,8);
% for t = 3:T
%     k = t-2;
%     Phi1(k,:) = [-y1_meas_wn(t-1,sim), -y1_meas_wn(t-2,sim), ...
%         -y2_meas_wn(t-1,sim), -y2_meas_wn(t-2,sim), ...
%         u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
%     Phi2(k,:) = [-y2_meas_wn(t-1,sim), -y2_meas_wn(t-2,sim), ...
%         -y1_meas_wn(t-1,sim), -y1_meas_wn(t-2,sim), ...
%         u1(t-1), u1(t-2), u2(t-1), u2(t-2)];
% end
% Gamma_hat = ( [Phi1, zeros(T-2,8) ; zeros(T-2,8), Phi2].' * ...
%     [Phi1, zeros(T-2,8) ; zeros(T-2,8), Phi2] ) / T;   % 16×16
% Sigma_theo = noiseStd1^2 * inv(Gamma_hat);
%
% fprintf('Rel. Frobenius-Fehler Sigma_WN vs. Theorie: %.3g\n', ...
%     norm(Sigma_WN-Sigma_theo,'fro')/norm(Sigma_theo,'fro'));

% %% Additional Plot
% paramIdx = 5;
% figure;
% histogram(scaledErr_WN(:,paramIdx),40,'Normalization','pdf');
% hold on;
% x = linspace(min(scaledErr_WN(:,paramIdx)), ...
%     max(scaledErr_WN(:,paramIdx)),200);
% pdf_est = normpdf(x, mu_WN(paramIdx), sqrt(Sigma_WN(paramIdx,paramIdx)));
% plot(x,pdf_est,'LineWidth',1.2); hold off;
% title(sprintf('Parameter %d – White noise',paramIdx));
% xlabel('$\Delta\theta_i = \sqrt{T}(\hat\theta_{T,i}-\theta_i)$');
% grid on;
%
% figure;
% qqplot(scaledErr_WN(:,paramIdx));
% title(sprintf('QQ-Plot  Parameter %d – White noise',paramIdx));

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

function [u_opt,u_df] = design_input(A,B,k,sigma_u,seed,nStarts)
% DESIGN_INPUT_MSRAND  –  globales Input-Design per MultiStart (SQP)
%   [u_opt,u_df] = design_input_MSrand(A,B,k,sigma_u,seed,nStarts)
%
%   A,B        – Systemmatrizen
%   k          – Periodenlänge (Anzahl äquidistanter Frequenzen)
%   sigma_u    – gewünschte Eingangsvarianz
%   seed       – RNG-Seed  (Default = 0)                → rng(22+seed)
%   nStarts    – Anzahl MultiStart-Punkte (Default = 10)
%
%   u_opt      – optimales Zeitbereichssignal  (k×d_u)
%   u_df       – optimaler Punkt im Frequenzbereich    (d_u·k × 1)

% --------------------- Default-Werte & RNG ---------------------------
if nargin < 5
    seed = 0;
end
if nargin < 6
    nStarts = 1;
end               % #Zufallspunkte
rng(22 + seed);                                % reproduzierbar

% --------------------- Problem-Größen -------------------------------
d_u   = size(B,2);                 % #Eingänge
d     = d_u * k;                   % Optimierungsdimension
omega = (2*pi/k) * (0:k-1);        % diskrete Kreisfrequenzen

% --------------------- Ziel & Nebenbedingungen ----------------------
objective = @(U) -min_eig(U,A,B,omega,k,d_u);        % maximiert λ_min
nonlcon   = @(U)  energy_con(U,sigma_u,k);            % Energie ≤ k·σᵤ

Aeq = zeros(d_u,d_u*k);                              % u(:,1) = 0
for j = 1:d_u
    Aeq(j,(j-1)*k+1) = 1;
end
beq = zeros(d_u,1);

opts = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
        'MaxIterations',20000,'MaxFunctionEvaluations',20000);

% --------------------- MultiStart-Setup ------------------------------
problem = createOptimProblem('fmincon', ...
          'objective', objective,         ...
          'x0',        zeros(d,1),        ...  % Dummy – wird ersetzt
          'Aeq',       Aeq, 'beq', beq,    ...
          'nonlcon',   nonlcon,           ...
          'options',   opts);

% Zufallsstartpunkte IN der Box ±sqrt(k*sigma_u)
rs = RandomStartPointSet('NumStartPoints', nStarts, ...
                         'ArtificialBound', sqrt(k*sigma_u));

ms = MultiStart('UseParallel', true, 'Display', 'iter');

[u_df, ~] = run(ms, problem, rs);     % global bestes Ergebnis

% --------------------- Zeitbereich-Signal rekonstruieren -------------
u_opt = zeros(k,d_u);
Umat  = reshape(u_df, d_u, k);
for j = 1:d_u
    u_t = real(ifft(Umat(j,:), k)) * sqrt(k);      % Parseval-Skalierung
    u_opt(:,j) = u_t / std(u_t,1) * sqrt(sigma_u); % Varianz-Normierung
end
end
% ---------------------------------------------------------------------
function lambda_min = min_eig(U,A,B,omega,k,d_u)
n   = size(A,1);
U   = reshape(U, d_u, k);
Gam = zeros(n);
for j = 1:k
    z = exp(1i*omega(j));
    H = (z*eye(n) - A) \ B;
    Gam = Gam + real(H*U(:,j)*U(:,j).'*H');
end
lambda_min = min(eig(Gam));
end
% ---------------------------------------------------------------------
function [c,ceq] = energy_con(U,sigma_u,k)
c   = sum(U(:).^2) - k*sigma_u;   % ≤ 0  … Energiebeschränkung
ceq = [];                         % keine Gleichungen
end
