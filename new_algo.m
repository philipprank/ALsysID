clear;
clc;
close all;
rng(24)
%% Parameters
% ARX model
p = 2;
q = 2;
d_y = 2;
d_u = 2;
n = p*d_y + q*d_u;

sigma_u = 1;           % Standard deviation input noise
sigma_e = sqrt(0.05);  % Standard deviation measurement noise
T = 1000;              % Number of observations
T0 = 100;              % Initial random excitation
k0 = 50;               % Number of frequency components
num_sim = 100;         % Number of simulations
update_interval = 50;

%% System Initialization
isStable = false;
while ~isStable
    a1 = randn(d_y); % zufällige Matrix
    a2 = randn(d_y);
    % Erstelle Companion–Matrix für die Ausgänge:
    A_comp = [ -a1, -a2; eye(d_y), zeros(d_y) ];
    eigenvalues = eig(A_comp);
    if max(abs(eigenvalues)) < 1
        isStable = true;
    else
        scale = 1.1 * max(abs(eigenvalues));
        a1 = a1 / scale;
        a2 = a2 / scale;
    end
end
% b–Matrizen (hier ohne Stabilitätsbedingungen)
b1 = randn(d_y, d_u);
b2 = randn(d_y, d_u);

% Wahres ARX–Parameter–Matrix (Beachte: in der Regression taucht -a_i auf)
theta_true = [-a1; -a2; b1; b2];  % Dimension: (p*d_y+q*d_u) x d_y

%% Zustandsraummodell für die Kovariaten (Companion–Form des ARX–Modells)
% Der Kovariaten–Vektor lautet: x_t = [ y_{t-1}; y_{t-2}; u_{t-1}; u_{t-2} ]
% Wir zerlegen den Zustand in den Ausgangs–Teil und den Eingangs–Teil:
n_y = p*d_y;
n_u = q*d_u;

% A11: Companion–Matrix für die Ausgänge (für p = 2)
A11 = [-a1, -a2; eye(d_y), zeros(d_y)];
% A12: Einfluss der vergangenen Eingänge in die Ausgangsgleichung;
% nur die erste Blockzeile erhält die b–Matrizen
A12 = [[b1, b2]; zeros(d_y*(p-1), q*d_u)];
% A22: Verschiebungsmatrix für die Eingänge (für q = 2)
A22 = [zeros(d_u, d_u), zeros(d_u, d_u); eye(d_u), zeros(d_u, d_u)];

% Vollständige Zustandsmatrix A für die Kovariaten–Entwicklung:
A_arx = [A11, A12; zeros(n_u, n_y), A22];

% Eingangsmatrix B für das Kombinierte Rauschen und den gesteuerten Eingang:
% v_t = [ ε_t; u_t ] mit Dimension (d_y + d_u)
% B1: ε_t wirkt nur in den ersten d_y–Zeilen des Ausgangsteils
B1 = [sigma_e * eye(d_y); zeros(n_y - d_y, d_y) ];
% B2: u_t wird direkt in den Eingangs–Teil injiziert (als neuester Eingang)
B2 = [ eye(d_u); zeros(n_u - d_u, d_u) ];
B_arx = [ B1, zeros(n_y, d_u);
    zeros(n_u, d_y), B2 ];

%% Schätzung mit Zufallseingang (ARX)
err_all_sim_rand = zeros(num_sim, T/update_interval);

for sim = 1:num_sim
    % Zufällige Eingänge
    u_rand = sqrt(sigma_u)*randn(T, d_u);
    x = zeros(n, T+1);  % Kovariaten–Vektor
    y = zeros(T, d_y);  % Beobachtete Ausgänge (y_t = erste d_y Elemente von x_{t+1})
    errors_random = [];

    for t = 1:T
        % v_t = [ε_t; u_t]
        epsilon_t = sigma_e * randn(d_y, 1);
        u_t = u_rand(t, :)';
        v_t = [epsilon_t; u_t];
        x(:,t+1) = A_arx*x(:, t) + B_arx*v_t;
        y(t, :) = x(1:d_y, t+1)';  % Ausgabe: die ersten d_y Elemente

        if mod(t, update_interval) == 0
            % Least-Squares–Schätzung:
            % Verwende alle bisherigen Daten: x_t als Kovariaten und y_t als Ausgabe.
            X_ls = x(:, 1:t);      % (n x t)
            Y_ls = y(1:t, :)';      % (d_y x t)
            theta_hat = pinv(X_ls * X_ls') * (X_ls * Y_ls');
            error_norm = norm(theta_true - theta_hat, 'fro');
            errors_random = [errors_random, error_norm];
        end
    end
    err_all_sim_rand(sim, :) = errors_random;
end

%% Schätzung mit aktiv gestaltetem Eingang (ARX)
errors_all_simulations_active = zeros(num_sim, T/update_interval);
u_df_all = cell(num_sim, 1);
final_costs = zeros(num_sim, 1);

for sim = 1:num_sim
    % Zunächst zufällige Excitation für T0 Zeitschritte
    u_rand_0 = sqrt(sigma_u) * randn(T0, d_u);
    x = zeros(n, T0+1);
    y = zeros(T0, d_y);
    errors_active = [];

    for t = 1:T0
        epsilon_t = sigma_e * randn(d_y, 1);
        u_t = u_rand_0(t, :)';
        v_t = [ epsilon_t; u_t ];
        x(:, t+1) = A_arx * x(:, t) + B_arx * v_t;
        y(t, :) = x(1:d_y, t+1)';
        if mod(t, update_interval) == 0
            X_ls = x(:, 1:t);
            Y_ls = y(1:t, :)';
            theta_hat = pinv(X_ls * X_ls') * (X_ls * Y_ls');
            error_norm = norm(theta_true - theta_hat, 'fro');
            errors_active = [errors_active, error_norm];
        end
    end

    % Aktiv gestaltetes Eingangssignal mittels Optimierung (wie im LTI–Fall)
    % Hier wird nun d_y als zusätzlicher Parameter übergeben.
    [u_opt, u_df_all{sim}, final_costs(sim)] = update_inputs(A_arx, B_arx, sigma_u, k0, n, d_u, d_y, k0, sim);
    % Anwenden des designten Eingangssignals für t = T0+1 bis T
    for t = T0+1:T
        epsilon_t = sigma_e * randn(d_y, 1);
        % Periodische Anwendung: u_opt hat Periode k0
        u_t = u_opt(mod(t-1, k0)+1, :)';
        v_t = [ epsilon_t; u_t ];
        x(:, t+1) = A_arx * x(:, t) + B_arx * v_t;
        y(t, :) = x(1:d_y, t+1)';
        if mod(t, update_interval) == 0
            X_ls = x(:, 1:t);
            Y_ls = y(1:t, :)';
            theta_hat = pinv(X_ls * X_ls') * (X_ls * Y_ls');
            error_norm = norm(theta_true - theta_hat, 'fro');
            errors_active = [errors_active, error_norm];
        end
    end

    % Falls die Länge der Fehlervektoren variiert, erfolgt hier eine Auffüllung
    if length(errors_active) ~= size(errors_all_simulations_active, 2)
        errors_active_padded = zeros(1, size(errors_all_simulations_active, 2));
        errors_active_padded(1:length(errors_active)) = errors_active;
        errors_all_simulations_active(sim, :) = errors_active_padded;
    else
        errors_all_simulations_active(sim, :) = errors_active;
    end
end

%% Mittelwert der Schätzfehler über alle Simulationen
mean_errors_random = mean(err_all_sim_rand, 1);
mean_errors_active = mean(errors_all_simulations_active, 1);

%% Plots

% Schätzfehler
figure;
plot(update_interval:update_interval:T, mean_errors_random, 'LineWidth', 1.5, 'DisplayName', 'Zufallseingang');
hold on;
plot(update_interval:update_interval:T, mean_errors_active, 'LineWidth', 1.5, 'DisplayName', 'Designtes Eingangssignal');
xlabel('Iteration $T$', 'Interpreter', 'latex');
ylabel('$\|\hat{\theta}_{T} - \theta^{\star}\|_{F}$', 'Interpreter', 'latex');
legend('Location', 'northeast');
grid on;

% Finaler Wert der Zielfunktion (aktives Design)
figure;
plot(1:num_sim, final_costs, 'o-', 'LineWidth', 1.5);
xlabel('Simulation $M$', 'Interpreter', 'latex');
ylabel('$\lambda_{\textrm{min}}(\Gamma)$', 'Interpreter', 'latex');
title('Finaler Zielfunktionswert', 'Interpreter', 'latex');
grid on;

% Optimale Fourier–Koeffizienten
figure;
for sim = 1:min(2, num_sim)
    stem(u_df_all{sim}(1, :), 'DisplayName', ['Simulation ' num2str(sim)]);
    hold on;
end
xlabel('Index der Fourier–Koeffizienten', 'Interpreter', 'latex');
ylabel('$\check{u}$', 'Interpreter', 'latex');
title('Optimale Fourier–Koeffizienten', 'Interpreter', 'latex');
legend('show');
grid on;

%% Optimierungsfunktion
function [u_opt, u_df, final_cost] = update_inputs(A, B, gamma, k, n, d_u, d_y, k0, sim)
% Diese Funktion optimiert das Eingangssignal für das
% ARX–Zustandsraummodell.
rng(42 + sim);
cov_matrix = zeros(n, n);
omega = (2 * pi / k) * (0:k-1); % Frequenzen
% Startwert für Fourier–Koeffizienten: Dimension d_u x k
u_df0 = unifrnd(-10,10, d_u, k);

% Zielfunktion: Maximiere den kleinsten Eigenwert (bzw. minimiere dessen Negative)
objective = @(U) -compute_min_eigenvalue(U, A, B, omega, k, cov_matrix, d_u, d_y);

nonlcon = @(u_df) energy_constraint(u_df, gamma, k);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',6000,'MaxFunctionEvaluations',6000);
[u_df, final_cost] = fmincon(objective, u_df0, [], [], [], [], [], [], nonlcon, options);

% Erzeuge periodisches Signal aus den Fourier–Koeffizienten (IFFT)
t = (1:k)';
u_opt = zeros(k, d_u);
for i = 1:d_u
    u_tmp = real(ifft(u_df(i,:), k))*sqrt(k);
    u_opt(:,i) = u_tmp / std(u_tmp,1);
end
end

%% Funktion zur Berechnung des minimalen Eigenwerts
function min_eig = compute_min_eigenvalue(U, A, B, omega, k, cov_matrix, d_u, d_y)
% U: Fourier–Koeffizienten (d_u x k)
% Berechne Gamma_u als Summe über die Frequenzkomponenten
[n, ~] = size(A);
U = reshape(U, d_u, k);
Gamma_u = zeros(n);
for i = 1:k
    % Annahme: B = [B1, B2] wobei B1 die ersten d_y Spalten sind.
    % Somit entspricht der u–Einfluss dem Block B2:
    B_u = B(:, d_y+1:end);
    freq_gain = (exp(1j*omega(i))*eye(n) - A) \ B_u;
    Gamma_u = Gamma_u + real(freq_gain * U(:, i) * U(:, i)' * freq_gain');
end
Gamma_u = Gamma_u + cov_matrix;
min_eig = min(eig(Gamma_u));
end

%% Energie–Nebenbedingung für die Fourier–Koeffizienten
function [c, ceq] = energy_constraint(u_df, gamma, k)
% Energiebedingung: trace(u_df'*u_df) ≤ gamma * k
energy = trace(u_df' * u_df);
c = energy - gamma * k;
ceq = [];
end

% %%
% % Holen der Fourier-Koeffizienten für die gegebene Simulation
% u_df = u_df_all{sim};  % u_df ist eine d_u x k Matrix
% k = size(u_df, 2);     % Anzahl der Frequenzkomponenten (k)
% d_u = size(u_df, 1);   % Dimension des Eingangssignals (d_u)
%
% % Zeitpunkte (angenommen, k Zeitpunkte für die Rekonstruktion)
% t = (1:k)';
%
% % Initialisiere das rekonstruierten Signal (d_u x k)
% u_reconstructed = zeros(k, d_u);
%
% % Rekonstruktion des Signals für jedes Eingangssignal (jede Dimension d_u)
% for i = 1:d_u
%     % Rekonstruiere das Signal für die i-te Dimension des Eingangssignals
%     u_reconstructed(:, i) = real(ifft(u_df(i, :), k));  % IFFT für die i-te Dimension
% end
%
% % Das rekonstruierten Signal ist jetzt in u_reconstructed
% % Du kannst es jetzt plotten oder weiterverarbeiten
%
% % Plotten des rekonstruierten Signals
% figure;
% plot(t, u_reconstructed);
% xlabel('Zeit');
% ylabel('Amplitude');
% title('Rekonstruiertes Signal');
% legend(arrayfun(@(i) sprintf('Input %d', i), 1:d_u, 'UniformOutput', false));