%% Hauptskript: Systemidentifikation mit optimal gestalteten Eingängen von Anfang an
clear; 
close all; 
clc;
rng(44)

%% Parameter und Systemeinstellungen
n = 4;         % Systemdimension
m = 2;         % Eingangsdimension
sigma_u = 1;   % Eingangsvarianz
sigma_w = 0.05;% Rauschvarianz
T = 800;       % Gesamtzahl der Beobachtungen
k = T;         % Periodenlänge: k = T
update_interval = 10;
num_sim = 3 ;   % Anzahl Simulationsdurchläufe

%% Erzeugung zufälliger, stabiler Systemmatrizen
isStable = false;
while ~isStable
    A = randn(n);
    eigenvalues = eig(A);
    if max(abs(eigenvalues)) < 1
        isStable = true;
    else
        A = A / (1.1 * max(abs(eigenvalues))); % Normierung zur Stabilisierung
    end
end
B = randn(n, m);

%% Bestimme das optimale Eingangssignal (optimal input design)
u_df_all = cell(num_sim, 1);
final_costs = zeros(num_sim, 1);
optimal_inputs = cell(num_sim, 1);

for sim = 1:num_sim
    % Optimierung der Fourierkoeffizienten:
    [u_opt, u_df_all{sim}, final_costs(sim)] = update_inputs(A, B, sigma_u, n, m, T, sim);
    % u_opt hat nun die Größe T x m und wird als optimales Eingangssignal genutzt.
    optimal_inputs{sim} = u_opt;
end

%% Simulation mit optimalem Eingang von Anfang an
err_all_sim_opt = zeros(num_sim, T/update_interval);
for sim = 1:num_sim
    % Verwende das optimierte Eingangssignal (Größe: T x m)
    u = optimal_inputs{sim};
    
    % Initialisiere den Zustandsvektor
    x = zeros(n, T+1);
    
    errors_opt = [];
    for t = 1:T
        % Zustandsupdate: Systemdynamik plus Rauschen
        x(:, t+1) = A*x(:, t) + B*u(t, :)' + sigma_w * randn(n, 1);
        if mod(t, update_interval) == 0
            X = x(:, 1:t);
            Y = x(:, 2:t+1);
            A_hat = Y * pinv(X);  % Schätzung der Systemmatrix mittels Kleinste-Quadrate
            error_norm = norm(A - A_hat, 2);
            errors_opt = [errors_opt, error_norm];
        end
    end
    err_all_sim_opt(sim, :) = errors_opt;
end

%% Berechnung des Durchschnittsfehlers über alle Simulationen
avg_error_norm_opt = mean(err_all_sim_opt, 1);

%% Plot: Durchschnittlicher Schätzfehler über die Beobachtungen
figure;
plot(update_interval:update_interval:T, avg_error_norm_opt, 'o-', 'LineWidth', 2);
xlabel('Beobachtungen T');
ylabel('Durchschnittlicher Schätzfehler \(\|\hat{A} - A\|_2\)', 'Interpreter', 'latex');
title('Durchschnittlicher Schätzfehler (Optimaler Eingang, k = T)');
grid on;

%% --- FUNKTIONEN ---

%% Optimierungsfunktion zur Bestimmung der Fourierkoeffizienten und Zeitbereichseingabe
function [u_opt, u_df_mat, final_cost] = update_inputs(A, B, gamma, n, m, k, sim)
    % Um eine reelle, periodische Eingangsfolge zu erhalten, optimieren wir
    % über Fourierkoeffizienten u_df, die idealerweise symmetrisch sein sollten.
    rng(42 + sim);  % Reproduzierbare Zufallsstartwerte
    u_df0 = unifrnd(-10, 10, m, k);
    u0 = u_df0(:); % Vektorisierung
    
    % Frequenzen (nur für Zielfunktion benötigt)
    omega = (2*pi/k) * (0:k-1);
    cov_matrix = zeros(n); % Dummy-Matrix für Vorab-Daten

    % Zielfunktion: Maximierung des kleinsten Eigenwerts der Fisher-Informationsmatrix
    objective = @(u_vec) objective_fun(u_vec, A, B, omega, k, n, m, cov_matrix);

    % Nichtlineare Nebenbedingungen: Energie-Constraint + evtl. Symmetrie
    nonlcon = @(u_vec) nonlcon_all(u_vec, gamma, k, m);
    
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
        'MaxIterations', 3000, 'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [u_vec_opt, final_cost] = fmincon(objective, u0, [], [], [], [], [], [], nonlcon, options);
    
    % Reshape optimierte Fourierkoeffizienten in m x k Matrix
    u_df_mat = reshape(u_vec_opt, m, k);
    
    % Manuelle IDFT zur Rekonstruktion des Eingangssignals im Zeitbereich
    u_opt = zeros(k, m);
    for iChan = 1:m
        for nIdx = 1:k
            val = 0;
            for ell = 0:(k-1)
                alpha = 2*pi * ell * (nIdx-1) / k;
                val = val + u_df_mat(iChan, ell+1) * exp(1j * alpha);
            end
            u_opt(nIdx, iChan) = real(val) / k;
        end
    end
end

%% Zielfunktion: Negatives Minimum des Eigenwerts der Fisher-Informationsmatrix
function f = objective_fun(u_vec, A, B, omega, k, n, m, cov_matrix)
    u_df = reshape(u_vec, m, k);
    Gamma_u = zeros(n);
    for i = 1:k
        freq_gain = (exp(1j*omega(i))*eye(n) - A) \ B;
        Gamma_u = Gamma_u + real(freq_gain * (u_df(:, i)*u_df(:, i)') * freq_gain');
    end
    Gamma_u = Gamma_u + cov_matrix;
    min_eig_val = min(eig(Gamma_u));
    f = -min_eig_val;  % Maximierung des kleinsten Eigenwerts (Minimierung von -lambda_min)
end

%% Nichtlineare Nebenbedingungen: Energie-Constraint und (optionale) Symmetrie
function [c, ceq] = nonlcon_all(u_vec, gamma, k, m)
    u_df = reshape(u_vec, m, k);
    energy = sum(u_df(:).^2);
    c = energy - gamma * k;  
    ceq = [];
    for i = 1:m
        for j = 2:(k/2)
            ceq = [ceq; u_df(i, j) - u_df(i, k - j + 2)];
        end
    end
end
