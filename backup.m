clear;
clc;
close all;
rng(24)
%% System Matrices
% System 1
A1 = [0.9 -0.5; 0 0];
B1 = [0; 1];
% System 2
A2 = [0.9 -0.2 0.3; 0 0 0; 0 0 0];
B2 = [0 0; 1 0; 0 1];
% System 3
A3 = [0.9 -0.6 0.3; 0.4 0.5 0.1; 0 0 0];
B3 = [0; 0; 1];

%% Parameters
T = 2000;                  % Number of observations
num = 10000;               % Number of simulations
theta_true = [0.9; 0.2];   % True ARX parameters
sigma_u = 1;               % Standard deviation input noise
sigma_e = sqrt(0.05);      % Standard deviation measurement noise
error_norm = zeros(num, 1);

%% System Initialization
A = A1;
B = B1;
disp(['Eigenvalues of A ', mat2str(eig(A))]);

Gram_A = zeros(size(A));
Gram_B = zeros(size(A));
for s = 0:T-1
    Gram_A = Gram_A + A^s*(A^s)';
    Gram_B = Gram_B + (A^s*B)*(A^s*B)';
end
Gram = Gram_A + Gram_B;

k = 10;
Gram_A_k = zeros(size(A));
Gram_B_k = zeros(size(A));
for s = 0:k-1
    Gram_A_k = Gram_A_k + A^s*(A^s)';
    Gram_B_k = Gram_B_k + (A^s*B)*(A^s*B)';
end
Gram_k = Gram_A_k + Gram_B_k;


%% Simulation
for i = 1:num
    u = sigma_u*randn(T,1);
    epsilon = sigma_e*randn(T,1);
    y = zeros(T,1);

    % Simulation of ARX-Model
    for t = 2:T
        y(t) = theta_true(1)*y(t-1) + theta_true(2)*u(t-1) + epsilon(t);
    end

    % Parameter Estimation (OLS)
    X = [y(1:T-1), u(1:T-1)];
    Y = y(2:T);
    theta_hat = X\Y;

    % 2-Norm of Estimation Error
    error_norm(i) = norm(theta_hat - theta_true);
end

%% Analysis
% Average Error
mean_error = mean(error_norm);
disp(['Average error over ', num2str(num), ' trials: ', num2str(mean_error)]);

% Optional: Dispersion of the errors (Standard deviation)
std_error = std(error_norm);
disp(['Standard deviation of the error: ', num2str(std_error)]);

% Determine the 95% quantile of the error norm, i.e., 95% of all estimates are below this value
error_quantile_95 = prctile(error_norm, 95);
disp(['95% quantile of the error norm: ', num2str(error_quantile_95)]);

%% Histogram of the error norm from all trials with display of the 95% threshold
figure;
histogram(error_norm, 30);
xlabel('Error norm $\|\hat{\theta} - \theta\|_{2}$', 'Interpreter', 'latex');
ylabel('Frequency', 'Interpreter', 'latex');
title('Histogram of error norms from 1000 trials', 'Interpreter', 'latex');
hold on;
% Plot the 95% threshold as a vertical line
yl = ylim;
plot([error_quantile_95 error_quantile_95], yl, 'r--', 'LineWidth', 1.5);
%legend('Error norm', '$95\%$ quantile', 'Interpreter', 'latex');
xline(0.0290, 'b--', 'LineWidth', 1.5);
ax = gca;
ax.TickLabelInterpreter = 'latex';

%% Active Input Design
m = 1;
k = 10;
%freqs = [1 2 3 4 5];
freqs = 0:k-1;map
u = generate_periodic_input(m, k, freqs);
u_fft = fft(u, [], 2);
gamma_squared = norm(u, 'fro')^2 / k;

Gamma = compute_stationary_covariance(A, B, u_fft, gamma_squared);
disp(Gamma);

function Gamma = compute_stationary_covariance(A, B, u_fft, gamma_squared)
[n, ~] = size(B);
k = size(u_fft, 2);
Gamma = zeros(n, n);

for ell = 0:(k-1)
    z = exp(1i*2*pi*ell / k);
    M = inv(z*eye(n) - A);

    u_l = u_fft(:, ell + 1);
    term = M * B * (u_l * u_l') * B' * M';
    Gamma = Gamma + term;
end

Gamma = Gamma / (gamma_squared * k^2);
end

function u = generate_periodic_input(m, k, freqs)
t = 0:(k-1);
u = zeros(m, k);

for i = 1:m
    for f = freqs
        a_f = randn();
        b_f = randn();
        u(i, :) = u(i, :) + a_f * cos(2*pi*f*t/k) + b_f * sin(2*pi*f*t/k);
    end
end
end


[u_opt, u_df_all{sim}, final_costs(sim)] = update_inputs(A_arx, B_arx, sigma_u, k0, n, d_u, d_y, k0, sim);