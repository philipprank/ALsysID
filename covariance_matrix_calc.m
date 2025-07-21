clear; clc; close all
%% Parameter -------------------------------------------------------------
k = 50;
T = 3000;
% A = [0.9 -0.2; 0 0];
% B = [0; 1];
A = [0.8 -0.6 0.2 -0.8; 1 0 0 0; 0 0 0 0; 0 0 1 0];
B = [0; 0; 1; 0];

%% 1) Periodisches Anregungssignal u_t ------------------------------------
rng(22);                              % Reproduzierbarkeit
U = zeros(k,1);
phi = 2*pi*rand(k/2-1,1);             % Zufallsphasen
U(2:k/2) = exp(1j*phi);
U(k/2+2:end)= conj(flipud(U(2:k/2)));     % Hermitesymmetrie  ⇒  reelles Zeitsignal
u = ifft(U);              % genau k-periodisch, reell

%% 2) Simulation + empirische Kovarianz ----------------------------------
% Zustände:  x_t = [y_t ; y_{t-1} ; u_{t-1} ; u_{t-2}]
x  = zeros(4,1);                      % Startwert: alles 0
S  = zeros(4);                        % Akkumulator  Σ x_t x_t^T

for t = 1:T
    % --- a) Kovarianz-Akkumulation (x_t basiert auf Vergangenheit) ------
    S = S + x*x.';                  % reine Transposition  (kein Conj.)

    % --- b) Eingabe u_t für diesen Zeitschritt -------------------------
    u_t = u(mod(t-1,k)+1);            % Index t (MATLAB 1-basiert)

    % --- c) Zustandsschritt  x_{t+1} = A2·x_t + B2·u_t ------------------
    x = A*x + B*u_t;
end
Gamma_emp = S/T;                    % empirischer Schätzer  \hat Γ_T

%% 3) Theoretische Kovarianz  (Frequenzsumme) -----------------------------
Gamma_theo = zeros(4);
I4 = eye(4);
for ell = 0:k-1
    if ell==0 || ell==k/2  
        continue;
    end
    z = exp(1i*2*pi*ell/k);
    H = (z*I4-A)\B;
    Gamma_theo = Gamma_theo + real(H*U(ell+1)*U(ell+1)'*H');
end
Gamma_theo = Gamma_theo/(k^2);

%% 4) Ausgabe -------------------------------------------------------------
disp('Theoretische Γ_k (4×4):');
disp(Gamma_theo);

disp('Empirische  Γ_T (4×4):');
disp(Gamma_emp);

disp('Betragsfehler |Γ_emp - Γ_theo|:');
disp(abs(Gamma_emp - Gamma_theo));

fprintf('Max. Abweichung: %.3e\n', max(abs(Gamma_emp(:) - Gamma_theo(:))));