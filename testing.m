%% Hauptskript: Systemidentifikation mit zufälligen und aktiv gestalteten Eingängen
clear; close all; clc;
rng(44)

%% Parameter und Systemeinstellungen
n = 4;         % Systemdimension
m = 2;         % Eingangsdimension
sigma_u = 1;   % Eingangsvarianz
sigma_w = 0.05;% Rauschvarianz
T = 800;       % Gesamtzahl der Beobachtungen
T0 = 100;      % Anzahl der Beobachtungen für zufällige Anregung
k = 300;       % Periodenlänge bzw. Anzahl Frequenzkomponenten (muss gerade sein)
update_interval = 10;
num_sim = 2;   % Anzahl Simulationsdurchläufe
epsilon = randn(n,T,num_sim);
u_rand = sqrt(sigma_u/m) * randn(T,m,num_sim);
%% Erzeugung zufälliger, stabiler Systemmatrizen
isStable = false;
while ~isStable
    A = randn(n);
    eigenvalues = eig(A);
    if max(abs(eigenvalues)) < 1
        isStable = true;
    else
        A = A/(1.1*max(abs(eigenvalues))); % Normierung zur Stabilisierung
    end
end
B = randn(n,m);

%% Speicher für Ergebnisse
u_df_all = cell(num_sim, 1);
final_costs = zeros(num_sim, 1);

%% Simulation mit zufälligem Eingang (Baseline)
err_all_sim_rand = zeros(num_sim, T/update_interval);
for sim = 1:num_sim
    % Zufallseingang: Skaliert so, dass pro Zeitschritt die Varianz ~ sigma_u ist.
    x = zeros(n, T+1);
    errors_random = [];
    for t = 1:T
        x(:, t+1) = A*x(:, t) + B*u_rand(t,:,sim)' + sigma_w*epsilon(:,t,sim);
        if mod(t, update_interval) == 0
            X = x(:, 1:t);
            Y = x(:, 2:t+1);
            A_hat = Y * pinv(X);  % Kleinste-Quadrate-Schätzer
            error_norm = norm(A - A_hat, 2);
            errors_random = [errors_random, error_norm];
        end
    end
    err_all_sim_rand(sim, :) = errors_random;
end

%% Simulation mit aktiv gestalteten Eingängen
errors_all_sim_active = zeros(num_sim, T/update_interval);
for sim = 1:num_sim
    % Phase 1: Zufällige Anregung für T0 Schritte
    % u_rand_0 = sqrt(sigma_u/m) * randn(T0, m);
    u_rand_0 = u_rand(1:T0,:,sim);
    x = zeros(n, T0+1);
    errors_active = [];
    for t = 1:T0
        x(:, t+1) = A*x(:, t) + B*u_rand_0(t, :)' + sigma_w*epsilon(:,t,sim);
        if mod(t, update_interval) == 0
            X = x(:, 1:t);
            Y = x(:, 2:t+1);
            A_hat = Y * pinv(X);
            error_norm = norm(A - A_hat, 2);
            errors_active = [errors_active, error_norm];
        end
    end
    
    % Phase 2: Optimierung der Eingangs-Fourierkoeffizienten
    gamma = sigma_u; % Gewünschte Energie pro Zeitschritt im Zeitbereich
    [u_opt, u_df_all{sim}, final_costs(sim)] = update_inputs(A, B, gamma, n, m, k, sim);
    
    % Phase 3: Anwendung des optimierten periodischen Eingangs
    % u_opt ist ein k x m-Signal, das periodisch wiederholt wird.
    x = zeros(n, T+1);
    for t = T0+1:T
        % Periodische Anwendung: Index im Bereich 1:k
        idx = mod(t-1, k) + 1;
        x(:, t+1) = A*x(:, t) + B*u_opt(idx, :)' + sigma_w*epsilon(:,t);
        if mod(t, update_interval) == 0
            X = x(:, 1:t);
            Y = x(:, 2:t+1);
            A_hat = Y * pinv(X);
            error_norm = norm(A - A_hat, 2);
            errors_active = [errors_active, error_norm];
        end
    end
    
    % Falls Länge der Fehler-Vektoren nicht übereinstimmt, auffüllen:
    if length(errors_active) ~= size(errors_all_sim_active, 2)
        errors_padded = zeros(1, size(errors_all_sim_active, 2));
        errors_padded(1:length(errors_active)) = errors_active;
        errors_all_sim_active(sim, :) = errors_padded;
    else
        errors_all_sim_active(sim, :) = errors_active;
    end
end

%% Mittelwerte der Schätzfehler über die Simulationen
mean_errors_random = mean(err_all_sim_rand, 1);
mean_errors_active = mean(errors_all_sim_active, 1);

%% Plots
% Plot of the estimation error over iterations
figure;
plot(update_interval:update_interval:T, mean_errors_random, 'LineWidth', 1.5, ...
     'DisplayName', '$u_{rand}$'); % Randomly generated input signal
hold on;
plot(update_interval:update_interval:T, mean_errors_active, 'LineWidth', 1.5, ...
     'DisplayName', '$u_{opt}$'); % Optimized input signal
xlabel('Observations $T$', 'Interpreter', 'latex'); % Label for the x-axis (iterations)
ylabel('$\|\hat{\theta}_T - \theta\|_2$', 'Interpreter', 'latex'); % Label for the y-axis (estimation error)
legend('Location', 'northeast'); % Place the legend in the top-right corner
grid on; % Enable grid for better visualization

% Plot of the final values of the objective function (active design)
figure;
plot(1:num_sim, final_costs, 'o-', 'LineWidth', 1.5);
xlabel('Simulation $M$', 'Interpreter', 'latex'); % Label for the x-axis (simulation index)
ylabel('$\lambda_{\textrm{min}}(\Gamma)$', 'Interpreter', 'latex'); % Label for the y-axis (minimum eigenvalue of the Fisher information matrix)
title('Final Objective Function Value', 'Interpreter', 'latex'); % Title of the plot
grid on; % Enable grid

% Plot of the optimized Fourier coefficients (first row for each simulation)
figure;
for sim = 1
    stem(u_df_all{sim}(1,:), 'DisplayName', ['Simulation ' num2str(sim)]); % Fourier coefficients for the 1st input channel
    hold on;
end
xlabel('Frequency index', 'Interpreter', 'latex'); % Label for the x-axis (frequency index)
ylabel('$\check{u}$', 'Interpreter', 'latex'); % Label for the y-axis (Fourier coefficient values)
title('Optimized Fourier Coefficients (1st Input Channel)', 'Interpreter', 'latex'); % Title of the plot
legend('show'); % Show the legend
grid on; % Enable grid

% Spectrum of the optimized input signal
fft_signal = abs(fft(u_opt(:, 1))); % Compute the FFT of the first input channel
figure;
stem(0:k-1, fft_signal); % Plot the magnitude of the FFT
xlabel('Frequency index'); % Label for the x-axis (frequency index)


%% --- FUNKTIONEN ---

%% Optimierungsfunktion zur Bestimmung der Fourierkoeffizienten und Zeitbereichseingabe
function [u_opt, u_df_mat, final_cost] = update_inputs(A, B, gamma, n, m, k, sim)
    % Um eine reelle, periodische Eingangsfolge zu erhalten, optimieren wir
    % über Fourierkoeffizienten u_df, die idealerweise symmetrisch sein sollten.
    % Zur Reduktion der Variablenzahl parameterisieren wir u_df als reellen Vektor
    % der Länge (m*k) und erzwingen ggf. im nonlcon Symmetrie.

    rng(42 + sim);  % reproduzierbare Zufallsstartwerte
    u_df0 = unifrnd(-10, 10, m, k);
    u0 = u_df0(:); % Vektorisierung
    
    % Frequenzen (nur für Zielfunktion benötigt)
    omega = (2*pi/k) * (0:k-1);
    cov_matrix = zeros(n); % evtl. Vorab-Daten summieren, hier dummy
    
    % Zielfunktion: Wir maximieren das minimale Eigenwertmaß -> -lambda_min
    objective = @(u_vec) objective_fun(u_vec, A, B, omega, k, n, m, cov_matrix);
    
    % Nichtlineare Nebenbedingungen: Energie-Constraint + evtl. Symmetrie
    nonlcon = @(u_vec) nonlcon_all(u_vec, gamma, k, m);
    
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
        'MaxIterations', 3000, 'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);
    [u_vec_opt, final_cost] = fmincon(objective, u0, [], [], [], [], [], [], nonlcon, options);
    
    % Reshape optimierte Fourierkoeffizienten in m x k Matrix
    u_df_mat = reshape(u_vec_opt, m, k);
    
    % --- Hier die manuelle (I)DFT ohne ifft ---
    u_opt = zeros(k, m);
    for iChan = 1:m
        for nIdx = 1:k
            val = 0;
            for ell = 0:(k-1)
                alpha = 2*pi * ell * (nIdx-1) / k;
                % Komplexer Exponentialanteil:
                val = val + u_df_mat(iChan, ell+1) * exp(1j * alpha);
            end
            % Durch k teilen, Realteil nehmen
            u_opt(nIdx, iChan) = real(val) / k;
        end
    end
    
    % --- Alte Variante mit ifft (auskommentiert) ---
    % for iChan = 1:m
    %     u_opt(:, iChan) = real(ifft(u_df_mat(iChan, :), k));
    % end
end

%% Zielfunktion: Negatives Minimum des Eigenwerts der "Informations"-Matrix
function f = objective_fun(u_vec, A, B, omega, k, n, m, cov_matrix)
    u_df = reshape(u_vec, m, k);
    Gamma_u = zeros(n);
    for i = 1:k
        freq_gain = (exp(1j*omega(i))*eye(n) - A) \ B;
        % Aufsummierung der Beiträge
        Gamma_u = Gamma_u + real(freq_gain * (u_df(:, i)*u_df(:, i)') * freq_gain');
    end
    Gamma_u = Gamma_u + cov_matrix;  % ggf. Vorab-Daten
    min_eig_val = min(eig(Gamma_u));
    f = -min_eig_val;  % "Maximiere kleinsten Eigenwert" -> Minimieren von -lambda_min
end

%% Nichtlineare Nebenbedingungen
% (a) Energieconstraint: sum(u_df.^2) <= gamma * k
% (b) (Optional) Symmetrie: für j=2,...,k/2 => u_df(:, j)=u_df(:, k-j+2)
function [c, ceq] = nonlcon_all(u_vec, gamma, k, m)
    u_df = reshape(u_vec, m, k);
    energy = sum(u_df(:).^2);
    % Energieconstraint
    c = energy - gamma * k;  
    % Symmetrie (falls erwünscht; hier nur demonstrativ)
    ceq = [];
    for i = 1:m
        for j = 2:(k/2)
            ceq = [ceq; u_df(i, j) - u_df(i, k - j + 2)];
        end
    end
end
