%% ------------------------------------------------------------------------
%%  Active Input Design – LS-Schätzung nur für A (B bekannt)      main_A.m
%% ------------------------------------------------------------------------
clear; clc; close all
rng default                                   % Reproduzierbarkeit

%% --------------------- 1) MODELL-DEFINITION -----------------------------
% ----- SISO, Ordnung n = 2 ----------------------------------------------
A1 = [ 0.9  -0.4 ;
    1.0   0   ];
B1 = [ 0.2 ;
    0.5 ];

% ----- SISO, Ordnung n = 3 ----------------------------------------------
A2 = [ 0.8  -0.3   0.2 ;
    1.0   0     0   ;
    0     1.0   0   ];
B2 = [ 0.4 ;
    0   ;
    0   ];

% ----- MIMO (2 × 2), Ordnung n = 4 --------------------------------------
A3 = [ 0.70   0.15  -0.05   0.00 ;   % dominante Eigenwerte 0.7 ± 0.1j
    1.00   0.00   0.00   0.00 ;
    0.00   1.00   0.00   0.00 ;
    0.00   0.00   1.00   0.00 ];

B3 = [ 0.25   0.10 ;   % Eingang 1 → stärkerer Einfluss auf x1, x2
    0.10   0.20 ;
    0.00   0.30 ;   % Eingang 2 → wirkt eher auf x3, x4
    0.00   0.40 ];

A4 = [  0.72   0.10  -0.05   0.00   0.00 ;   % dominante EV: 0.72 ± 0.12j
    1.00   0.00   0.00   0.00   0.00 ;   % + zwei Realteile 0.6, 0.5
    0.00   1.00   0.00   0.00   0.00 ;
    0.00   0.00   1.00   0.00   0.00 ;
    0.00   0.00   0.00   1.00   0.00 ];

B4 = [  0.30   0.15   0.00   0.05 ;   % Eingang 1 & 2 wirken stark auf x₁–x₂
    0.12   0.25   0.10   0.00 ;
    0.00   0.20   0.30   0.10 ;   % Eingang 3 dominiert x₃–x₅
    0.00   0.00   0.25   0.20 ;
    0.00   0.00   0.15   0.35 ];

%% --------------------- 2) MODELL WÄHLEN ---------------------------------
% *** Hier (nur) A:=…, B:=… umschalten, z. B. A=A3, B=B3 für das MIMO-Beispiel
A = A3;                 %  A1 | A2 | A3
B = B3;                 %  B1 | B2 | B3
% ------------------------------------------------------------------------

%% --------------------- 3) DIMENSIONS/INFO -------------------------------
[n,~]  = size(A);          % Zustandsdimension
[~,du] = size(B);          % Eingangsdimension
theta_true = A(:);         % Nur A wird geschätzt
theta_len  = numel(theta_true);

fprintf('-------------------------------------------------------------\n')
fprintf('n  (states)   = %d\n',n)
fprintf('du (inputs)   = %d\n',du)
fprintf('|A|_vec       = %d Parameter\n',theta_len)
fprintf('-------------------------------------------------------------\n')

%% --------------------- 4) SIM- & MC-PARAMETER ---------------------------
T  = 400;                    % Datenlänge
M  = 20;                       % Monte-Carlo-Runs
k  = 100;                      % Periode „Oracle“-Signal
sigma_u = 1;                   % Varianz Eingang
sigma_e = sqrt(0.001);          % Messrauschen (auf Zuständen)

update_interval = 2;
ival  = update_interval:update_interval:T;
numSt = numel(ival);

Range      = [-1 1];
Band       = [0 1];            % idinput-Band
freq_band  = [0.01 0.5];       % Multisine-Band
numComp    = 4;

%% --------------------- 5) SPEICHER --------------------------------------
u_wn   = zeros(T,du,M);
u_opt  = zeros(T,du,M);
u_sine = zeros(T,du,M);

x_wn   = zeros(T,n,M);
x_opt  = zeros(T,n,M);
x_sine = zeros(T,n,M);

Ahat_wn   = zeros(theta_len,numSt,M);
Ahat_opt  = zeros(theta_len,numSt,M);
Ahat_sine = zeros(theta_len,numSt,M);

err_wn   = zeros(M,numSt);
err_opt  = zeros(M,numSt);
err_sine = zeros(M,numSt);

freq = zeros(M,numComp);                    % Multisine-Frequenzen
e    = sigma_e*randn(T,n,M);                % Messrauschen
x0   = zeros(n,1);                          % Startzustand

%% --------------------- 6) MONTE-CARLO-SCHLEIFE --------------------------
for sim = 1:M
    % ===== (a) White-Noise-Eingang =======================================
    u_wn(:,:,sim) = sqrt(1/du)*idinput([T du],'rgs',Band,Range);
    x_wn(:,:,sim) = simulate_states(u_wn(:,:,sim),A,B,e(:,:,sim),x0);

    [Ahat_wn(:,:,sim),err_wn(sim,:)] = ...
        estimate_A_ls(x_wn(:,:,sim),u_wn(:,:,sim),B,theta_true,ival);

    % ===== (b) Oracle-Eingang ============================================
    [u_per,U_full_loc] = design_input(A,B,k,sigma_u);
    u_opt(:,:,sim) = repmat(u_per,ceil(T/k),1);
    x_opt(:,:,sim) = simulate_states(u_opt(:,:,sim),A,B,e(:,:,sim),x0);

    [Ahat_opt(:,:,sim),err_opt(sim,:)] = ...
        estimate_A_ls(x_opt(:,:,sim),u_opt(:,:,sim),B,theta_true,ival);

    % ===== (c) Multisine-Eingang =========================================
    [u_ms,info] = multisine_generation(freq_band,1,T,numComp,du);
    scale = sqrt(sigma_u ./ var(u_ms));
    u_sine(:,:,sim) = u_ms .* scale;
    freq(sim,:)     = info(1).frequencies;

    x_sine(:,:,sim) = simulate_states(u_sine(:,:,sim),A,B,e(:,:,sim),x0);

    [Ahat_sine(:,:,sim),err_sine(sim,:)] = ...
        estimate_A_ls(x_sine(:,:,sim),u_sine(:,:,sim),B,theta_true,ival);
end

%% --------------------- 7) (optional) Gamma_theo -------------------------
Gamma_theo = build_Gamma(U_full_loc,A,B,k);  %#ok<NASGU>

%% --------------------- 8) PLOT – mittl. Schätzfehler --------------------
figure
semilogy(ival,mean(err_wn ,1),'LineWidth',1); hold on
semilogy(ival,mean(err_opt,1),'LineWidth',1)
semilogy(ival,mean(err_sine,1),'LineWidth',1); hold off
xlabel('T'), ylabel('||Â_T – A^*||_F'), grid on
legend('White noise','Oracle','Multisine')
title('Least-Squares-Schätzung nur für A (B bekannt)')

%% ========================================================================
%% ----------------------- HILFSFUNKTIONEN --------------------------------
function x_meas = simulate_states(u,A,B,e,x0)
% Evolviert x_{k+1}=A x_k + B u_k  und liefert gemessene Zustände
T = size(u,1); n = size(A,1);
x_true = x0(:);                  % Startwert
x_meas = zeros(T,n);
for k = 1:T
    x_meas(k,:) = (x_true + e(k,:).').';     % Messung mit Rauschen
    x_true      = A*x_true + B*u(k,:).';
end
end
% -------------------------------------------------------------------------
function [A_hat_store,err] = estimate_A_ls(x,u,B,A_true_vec,ival)
% Least-Squares-Schätzung für A mit bekanntem B
[n,~]  = size(x); numSt = numel(ival);
A_hat_store = zeros(numel(A_true_vec),numSt);
err         = zeros(1,numSt);

for k = 1:numSt
    N = ival(k);
    Xk   = x(1:N-1,:).';            % n × (N-1)
    Uk   = u(1:N-1,:).';            % du×(N-1)
    Xkp1 = x(2:N  ,:).';            % n × (N-1)

    A_hat = (Xkp1 - B*Uk) * pinv(Xk);
    A_hat_store(:,k) = A_hat(:);
    err(k) = norm(A_hat(:) - A_true_vec);
end
end
% -------------------------------------------------------------------------
function [u_opt,U_full] = design_input(A,B,k,sigma_u)
% Heuristisches „Oracle“-Input-Design
d_u = size(B,2); m = k/2-1;
x0  = randn(2*d_u*m,1)/sqrt(2);
obj = @(x) -min_eig(assemble_U(x,d_u,k),A,B,k);
con = @(x)  energy_con(assemble_U(x,d_u,k),sigma_u,k);
opts = optimoptions('fmincon','Alg','interior-point','Display','iter', ...
    'MaxIter',4000,'MaxFunEvals',4000);
x_opt  = fmincon(obj,x0,[],[],[],[],[],[],con,opts);
U_full = assemble_U(x_opt,d_u,k);
u_opt  = ifft(U_full,[],2,'symmetric').';
end
% -------------------------------------------------------------------------
function lam = min_eig(U,A,B,k)
n = size(A,1); omega = (2*pi/k)*(0:k-1); Gamma = zeros(n);
for j = 1:k
    z = exp(1i*omega(j));
    H = (z*eye(n)-A)\B;
    Gamma = Gamma + real(H*U(:,j)*U(:,j)'*H');
end
lam = min(eig(Gamma));
end
% -------------------------------------------------------------------------
function U = assemble_U(x,d_u,k)
m = k/2-1; x = reshape(x,2*d_u,m);
Upos = x(1:d_u,:) + 1i*x(d_u+1:end,:);
U = zeros(d_u,k);
U(:,2:m+1)     = Upos;
U(:,k/2+2:end) = conj(fliplr(Upos));
end
% -------------------------------------------------------------------------
function [c,ceq] = energy_con(U,sigma_u,k)
c   = sum(abs(U(:)).^2) - sigma_u*k^2;
ceq = [];
end
% -------------------------------------------------------------------------
function [y,info] = multisine_generation(freqBand,fs,T,nComp,nSig)
f0   = fs/T;
kMin = max(1,ceil(freqBand(1)/f0));
kMax = floor(freqBand(2)/f0);
if kMax<kMin,  error('Band ausserhalb Bereich'); end
if nComp>(kMax-kMin+1), error('Zu viele Komponenten'); end

info = repmat(struct('frequencies',[],'magnitude',[],'phase',[]),nSig,1);
magTemplate   = ones(1,nComp)/nComp;
phaseTemplate = pi*(0:nComp-1).*(1:nComp)/nComp;

t = (0:T-1)/fs;  Y = zeros(nSig,T);
for s = 1:nSig
    bins = datasample(kMin:kMax,nComp,'Replace',false,'Weights',1./(kMin:kMax));
    bins = sort(bins); freqs = bins*f0;

    info(s).frequencies = freqs;
    info(s).magnitude   = magTemplate;
    info(s).phase       = phaseTemplate;

    for j = 1:nComp
        Y(s,:) = Y(s,:) + magTemplate(j)* ...
            cos(2*pi*freqs(j)*t + phaseTemplate(j));
    end
    Y(s,:) = Y(s,:)/max(abs(Y(s,:)));
end
y = Y.';                         % (T × nSig)
end
% -------------------------------------------------------------------------
function Gamma = build_Gamma(U_df,A,B,k)
d_u = size(B,2); omega = (2*pi/k)*(0:k-1);
U = reshape(U_df,d_u,k); n = size(A,1);
Gamma = zeros(n);
for j = 1:k
    z = exp(1i*omega(j));
    H = (z*eye(n)-A)\B;
    Gamma = Gamma + real(H*U(:,j)*U(:,j)'*H');
end
Gamma = Gamma/(k^2);
end
% -------------------------------------------------------------------------
