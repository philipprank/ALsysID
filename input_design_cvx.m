% =====================================================================
%  active_input_design_siso1_matrix_trace.m
%  ---------------------------------------------------------------------
%  Active Input Design (SISO, ARX‑1) – Fisher‑Matrix via Spur‑Form (3.26)
%    · Optimales (Oracle‑) Signal per SDP (Variable U ∈ H₊^{k×k})
%    · White‑Noise‑Signal als Referenz
%    · Least‑Squares‑Schätzung θ̂_T und Fehlerkurve ‖θ̂_T−θ*‖₂
%  ---------------------------------------------------------------------
%  REQUIREMENTS
%    · MATLAB
%    · CVX + SDPT3 (oder SeDuMi), danach einmal cvx_setup ausführen
% =====================================================================

clear; clc; close all;
rng(10);

%% 1) System definition --------------------------------------------------
A = [0.95 -0.05; 0 0];    % 2‑state Realisation von ARX(1,1)
B = [0; 1];               % SISO
p = 1; q = 1;
true_theta = [0.95 -0.05];

%% 2) Experiment parameters ---------------------------------------------
T        = 1400;          % Trajektorienlänge
k        = 100;           % Periodenlänge
sigma_u  = 1;             % Eingangs‑Varianzbudget
sigma_e  = sqrt(0.001);   % Störstandardabweichung
firstIdx = max(p,q)+1;
ival     = firstIdx:10:T; % Checkpoints für Fehlerkurve

%% 3a) Oracle input via trace‑SDP ---------------------------------------
[u_per, U_opt] = design_input_sdp_matrix_trace(A,B,k,sigma_u);
u_oracle   = repmat(u_per, ceil(T/k), 1);
u_oracle   = u_oracle(1:T);

%% 3b) White‑noise input -------------------------------------------------
u_wn = randn(T,1);        % Varianz 1

%% 4) Simulate trajectories ---------------------------------------------
noise_or = sigma_e*randn(T,1);
noise_wn = sigma_e*randn(T,1);

y_or = simulate_ARX(u_oracle, true_theta, noise_or, p, q, firstIdx);
y_wn = simulate_ARX(u_wn,     true_theta, noise_wn, p, q, firstIdx);

%% 5) Parameter estimation ----------------------------------------------
[~, err_or] = estimate_ARX(y_or, u_oracle, true_theta, p, q, ival, firstIdx);
[~, err_wn] = estimate_ARX(y_wn, u_wn,     true_theta, p, q, ival, firstIdx);

%% 6) Plot parameter error curves ---------------------------------------
figure;
semilogy(ival, err_or,'b-','LineWidth',1.5); hold on;
semilogy(ival, err_wn,'r--','LineWidth',1.5);
grid on;
xlabel('T'); ylabel('|\thetâ_T − θ*|_2');
legend('oracle input','white noise','Location','southwest');
title('Parameter‑Schätzfehler über die Zeit');

%% 7) First 300 samples --------------------------------------------------
figure;
subplot(2,2,1); plot(u_oracle(1:300)); grid on; ylabel('u_{oracle}');
subplot(2,2,2); plot(u_wn(1:300));     grid on; ylabel('u_{wn}');
subplot(2,2,3); plot(y_or(1:300));     grid on; ylabel('y_{oracle}');
subplot(2,2,4); plot(y_wn(1:300));     grid on; ylabel('y_{wn}');
sgtitle('Inputs & Outputs (erste 300 Samples)');


% =====================================================================
%  Local functions (alles in derselben Datei erlaubt)
% =====================================================================

function [u_per, U_opt] = design_input_sdp_matrix_trace(A,B,k,sigma_u)
% GENERALISCHES TRACE‑SDP (SISO/MIMO, ARX‑1)
%   A (n×n), B (n×du), k Periodenlänge, sigma_u Eingangs‑Budget

[n,du] = size(B);

% 1) Frequenzgitter und M_ell vorab berechnen
omega = 2*pi*(0:k-1)/k;
Mcell = cell(k,1);
for ell = 1:k
    Mcell{ell} = (exp(1j*omega(ell))*eye(n) - A)\B;  % n×du
end

% 2) Q{i,j} aufbauen: block‑diagonal aus den k Beiträgen
Q = cell(n,n);
for i = 1:n
    for j = 1:n
        Q{i,j} = zeros(k,k);
        for ell = 1:k
            Hij = Mcell{ell}(i,:)*Mcell{ell}(j,:)' ;    % Skalar bei du=1
            Q{i,j}(ell,ell) = (1/k) * Hij;
        end
    end
end

% 3) CVX‑Formulierung des SDP
cvx_solver sdpt3
cvx_begin sdp
cvx_precision high
variable U(k,k) hermitian semidefinite
variable t                      % λ_min
variable Gamma(n,n) hermitian   % Fisher‑Matrix

% Aufbau von Gamma(i,j) = trace(U*Q{i,j})
for i = 1:n
    for j = 1:n
        Gamma(i,j) == real( trace( U * Q{i,j} ) );
    end
end

maximize( t )
subject to
Gamma - t*eye(n) == semidefinite(n);   % E‑Kriterium
trace(U)         <= sigma_u * k;       % Energie‑Budget
U(1,1)           == 0;                 % kein DC‑Anteil
cvx_end

U_opt = U;

% 4) Rank‑1‑Rekonstruktion und periodisches Signal
[vecU,lambda] = eigs(U_opt,1,'la');        % größter Eigenwert
u_hat = sqrt(real(lambda)) * vecU;         % k×1 (komplex)
% reelle, konjugiert symmetrie gewährleisten
u_full = [u_hat; conj(flipud(u_hat(2:end-1)))];
u_per  = real(ifft(u_full,'symmetric'));
u_per  = u_per(1:k);  % Länge k
end

function y = simulate_ARX(u,theta,e,p,q,firstIdx)
T = length(u); y = zeros(T,1);
for t = firstIdx:T
    phi   = [flipud(y(t-p:t-1)); flipud(u(t-q:t-1))];
    y(t)  = theta*phi + e(t);
end
end

function [theta_hat,err] = estimate_ARX(y,u,theta_true,p,q,ival,firstIdx)
T  = length(y); dx = p+q; N = T-firstIdx+1;
Phi = zeros(N,dx); Y = zeros(N,1);
for t = firstIdx:T
    j = t-firstIdx+1;
    Phi(j,:) = [flipud(y(t-p:t-1)); flipud(u(t-q:t-1))]';
    Y(j)     = y(t);
end
theta_hat = zeros(dx,length(ival)); err = zeros(size(ival));
for kidx = 1:length(ival)
    Ni = ival(kidx)-firstIdx+1;
    theta_est        = Phi(1:Ni,:) \ Y(1:Ni);
    theta_hat(:,kidx)= theta_est;
    err(kidx)        = norm(theta_est - theta_true');
end
end
