function [beta_val, vals] = betaA(A, K)
% BETA_A  Approximates the quantity 
%    beta(A) = sup_{k >= 0} ||A^k||_2^(1/2 + rho(A)/2)
% by evaluating ||A^k||_2 for k = 0..K and taking the maximum.
%
% Usage:
%    [beta_val, vals] = betaA(A, K)
%
% Inputs:
%    A - square matrix
%    K - maximum exponent to check (nonnegative integer)
%
% Outputs:
%    beta_val - approximate value of beta(A)
%    vals     - array of values ||A^k||_2^(alpha) for k = 0..K
%
% Notes:
%  - If rho(A) > 1, we return Inf for both beta_val and vals, since the 
%    supremum diverges.
%  - If rho(A) <= 1, we compute A^k for k=0..K and record the values.
%  - Increase K if you need a finer approximation or if your matrix is
%    close to having rho(A) = 1.

% 1) Compute spectral radius of A
eigsA = eig(A);
rA   = max(abs(eigsA));  % spectral radius

% 2) Define the exponent alpha = 1/2 + rho(A)/2
alpha = 0.5 + 0.5*rA;

% 3) If the spectral radius is > 1, beta(A) is unbounded
if rA > 1
    beta_val = Inf;
    vals = Inf;
    return
end

% 4) Otherwise, evaluate ||A^k||_2^(alpha) for k = 0..K
vals = zeros(K+1, 1);
I = eye(size(A));  % Identity matrix
Ak = I;            % Start with A^0 = I

for kIdx = 0:K
    normAk = norm(Ak, 2);
    vals(kIdx+1) = normAk^alpha;
    
    % Update A^k -> A^(k+1) for next iteration
    Ak = Ak * A;
end

% 5) The supremum over 0..K is our approximation
beta_val = max(vals);
end


[max, all] = betaA(A_margin, 20);