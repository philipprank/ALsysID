function [sys_id,data_id] = opt_input_cvrlx(sys,options,e)
%%==================================%%
%%%     ACTIVE INPUT DESIGN        %%%
%%%     (CONVEX REALAXATION)       %%%

% Active input design (offline) in the time
% domain by convex relaxation based on
% "Input Design for System Identification
% via Convex Relaxation", Manchsester 2010
% Assumes a-priori knowledge of the true
% parameters and does all calculation before
% simulation (time-consuming for slow systems
% or long observations)

%%% Function needs CVX (cvxr.com/cvx/)! %%%

%% Initialization
MC = options.MC;
T = options.T;
H = options.H;
interval = options.interval;
estimation = options.estimation;

%% Input-Matrix Generation
Q = cell(p,p);
for i = 1:p
    for j = 1:p
        Q{2*i-1,2*j-1} = sumup(p,T,du,dy,@compcoeff_R,@compcoeff_R,i,j,G);
        Q{2*i-1,2*j} = sumup(p,T,du,dy,@compcoeff_S,@compcoeff_R,i,j,G);
        Q{2*i,2*j-1} = sumup(p,T,du,dy,@compcoeff_R,@compcoeff_S,i,j,G);
        Q{2*i,2*j} = sumup(p,T,du,dy,@compcoeff_S,@compcoeff_S,i,j,G);
    end
end

%% Optimization
Q_cvx = cell2mat(Q);
I = zeros(2*p);

cvx_begin
variable U(T,T) symmetric
expression I(2*p,2*p)

for i = 1:2*p
    for j = 1:2*p
        idx_row = T*(i-1)+1:T*i;
        idx_col = T*(j-1)+1:T*j;
        I(i,j) = trace(U*Q_cvx(idx_row,idx_col));
    end
end
maximize(lambda_min(I));
subject to
diag(U) <= ones(T,1);
% trace(U) <= 1
U == semidefinite(T);
cvx_end

%% Optimal Input Generation
D = chol(U, 'upper');
Xi  = randn(T,du);

u_opt = diag(1)*sign(D'*Xi);
u_opt = u_opt';
% [V, D] = eig(U);
% eigenvalues = diag(D);
% [~, maxIndex] = max(eigenvalues);
% u_opt = V(:, maxIndex);


toc
end

%% Functions
function S = compcoeff_S(z,du,T,~,~)
S = zeros(du,T);
S(1,z) = 1;
end

function R = compcoeff_R(z,dy,T,G,p)
R = zeros(dy,T);
if z-p+1 < 1
    R(1,1:z) = flip([0 G(1:z-1)]);
else
    R(1,z-p+1:z) = flip([0 G(1:end-1)]);
end
end

function Q_elem = sumup(p,T,du,dy,compcoeff1,compcoeff2,i,j,G)
Q_elem = zeros(T,T);
for t = p+1:T
    In1 = compcoeff1(t-i,du,T,G,p);
    In2 = compcoeff2(t-j,dy,T,G,p);
    Q_elem = Q_elem + In1'*In2;
end
end

%% Cost Function with Different Optimality Criteria
function J = calc_cost(I_f,I,options)
if strcmp(options.opt,'a-opt') == true
    J = trace(I_f + I);
elseif strcmp(options.opt,'d-opt') == true
    J = log(det(I_f + I));
elseif strcmp(options.opt,'e-opt') == true
    J = min(eig(I_f + I));
else
    J = 'error';
end
end



