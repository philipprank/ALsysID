function [sys_id,data_id] = opt_input_rcdhz(sys,options,e)
%%==================================%%
%%%     ACTIVE INPUT DESIGN        %%%
%%%     (RECEDING HORIZON)         %%%

% Active input design by optimizing online
% some cost function of the Fisher information
% using a receeding horizon approach.

% Three different cost functions based on:
% 1. A-optimal design
% 2. D-optimal design
% 3. E-optimal design
% -> calc_cost(...)

% Two cases:
% 1. Using true parameters
%   -> theta = theta_true
% 2. Using recursive estaimted parameters
%   -> theta will be updated online
%      (initial guess needed, e.g. by
%      prior identification (PEM) with
%      PRBS, state-space identification
%      or initial guess)
% -> fisher_pred(...)

% Two algorithms
% 1. oe (true): At defined timesteps (interval)
%               all data are used for PEM
% 2. recusriveoe (recursive): At every timestep the
%                 parameters are updated recursively
%
%%==================================%%

%% Initialization
MC = options.MC;
T = options.T;
H = options.H;
interval = options.interval;
estimation = options.estimation;
init_param = options.init_param;

na = length(sys.A) - 1;
nb = length(sys.B) - 1;
theta = [sys.B(2:end) sys.A(2:end)]; % free coefficients
H_eff = H - nb; % (effective) horizon length for input design
u_sym = sym('u_sym',[1 H]);
y_sym = sym('y_sym', [1 na]);
theta_sym = sym('theta_sym',[1 na+nb]);

I = zeros(na+nb,na+nb,T); % sum of past Fisher information
I_f = zeros(na+nb); % sum of future (predicted) Fisher information
J = zeros(size(2^H_eff,1),T); % cost function
u_opt_H = zeros(T,H_eff); % optimal (predicted) input
u_opt = zeros(T,1); % optimal (applied) input (only 1st input of u_opt_H is used at any time step)
y_opt = zeros(T,1); % output when u_opt is used

data_id = cell(MC,1);
sys_id = zeros(floor(T/interval-50/interval)+1,na+nb+2,MC);
t_eval = interval:interval:T; % time stepas where parameters get estimated/evaluated

% System given as output-error model (structure is assumed to be known)
oemodel = @(u,y) theta_sym*[u -y]';

%% Preallocation of Input Signals
% For direct search: Initialization of u with
% all possible values ('alg1'); avoid possible
% local maxima
if strcmp(options.alg,'alg1') == true
    M = 0:2^(H_eff)-1;
    u_all = dec2bin(M,H_eff) - '0';
    u_all(u_all == 0) = -1;
end

%% Predicted Fisher Information
F_pred = fisher_pred(H,H_eff,u_sym,y_sym,oemodel,na,nb,theta_sym);

%% Algorithm (Optimal Input Design, True Parameters for Optimization)
for i = 1:MC
    if strcmp(estimation,'true') == true
        init_sys = make_init(sys);
        u_opt(1:max(na,nb)) = idinput(max(na,nb)); % PRBS for first inputs, which can't be optimized
        idx = 1;
        % Optimization
        for j = max(na,nb)+1:T
            if strcmp(options.alg,'alg1') == true
                for k = 1:length(M)
                    u_H = u_all(k,:);
                    I_f(:,:,k) = F_pred(theta,[u_opt(j-nb:j-1,1)' u_H],y_opt(j-1:-1:j-na,1)');
                    J(k,j) = calc_cost(I_f(:,:,k),I(:,:,j-1),options);
                    if k == 1
                        u_opt_H(j,:) = u_H;
                        comp = J(k,j); % comp: variable for internal comparison
                    else
                        if J(k,j) > comp
                            u_opt_H(j,:) = u_H;
                            comp = J(k,j);
                        end
                    end
                end
                u_opt(j) = u_opt_H(j,1);

            elseif strcmp(options.alg,'alg2') == true
                lb = -1*ones(H_eff);
                ub = 1*ones(H_eff);
                % u_opt_H = fmincon(calc_cost,idinput(H_eff),[],[],[],[],lb,ub);
                % u_opt(j) = u_opt_H(j,1);
            end

            y_opt(j) = sim_system(sys,u_opt(j-1:-1:j-nb),y_opt(j-1:-1:j-na),e(j,i));
            I_upd = update_cost(sys,u_opt(j:-1:j-nb+1),y_opt(j:-1:j-na+1));
            I(:,:,j) = I(:,:,j-1) + I_upd;

            if any(t_eval(:) == j) && j >= 50
                data_id{i,1} = iddata(y_opt(1:j),u_opt(1:j),1);
                id_struct = oe(data_id{i,1},init_sys);
                sys_id(idx,:,i) = [id_struct.B id_struct.F];
                idx = idx + 1;
            end
        end

        %% Algorithm (Optimal Input Design, True Parameters for Optimization)
    elseif strcmp(estimation,'recursive') == true
        idx = 1;
        B0 = init_param.B;
        A0 = init_param.A;
        obj_oe = recursiveOE([nb na 1],B0,A0);
        theta = [B0(2:end) A0(2:end)];
        calc_covar = zeros(1,length(theta));
        calc_covar(theta == 0) = 10^(-7); % assumed known system structure
        calc_covar(theta ~= 0) = 0.1;
        obj_oe.InitialParameterCovariance = diag(calc_covar);
        u_opt(1:max(na,nb)) = idinput(max(na,nb)); % PRBS for first inputs, which can't be optimized
        for j = max(na,nb)+1:T
            if strcmp(options.alg,'alg1') == true
                for k = 1:length(M)
                    u_H = u_all(k,:);
                    I_f(:,:,k) = F_pred(theta,[u_opt(j-nb:j-1,1)' u_H],y_opt(j-1:-1:j-na,1)');
                    J(k,j) = calc_cost(I_f(:,:,k),I(:,:,j-1),options);
                    if k == 1
                        u_opt_H(j,:) = u_H;
                        comp = J(k,j); % comp: variable for internal comparison
                    else
                        if J(k,j) > comp
                            u_opt_H(j,:) = u_H;
                            comp = J(k,j);
                        end
                    end
                end
                u_opt(j) = u_opt_H(j,1);

            elseif strcmp(options.alg,'alg2') == true
                lb = -1*ones(H_eff);
                ub = 1*ones(H_eff);
                % u_opt_H = fmincon(calc_cost,idinput(H_eff),[],[],[],[],lb,ub);
                % u_opt(j) = u_opt_H(j,1);
            end

            y_opt(j) = sim_system(sys,u_opt(j-1:-1:j-nb),y_opt(j-1:-1:j-na),e(j,i));
            I_upd = update_cost(sys,u_opt(j:-1:j-nb+1),y_opt(j:-1:j-na+1));
            I(:,:,j) = I(:,:,j-1) + I_upd;
            [B_id, F_id, ~] = step(obj_oe,y_opt(j),u_opt(j));
            theta = [B_id(2:end) F_id(2:end)];

            if any(t_eval(:) == j) && j >= 50
                data_id{i,1} = iddata(y_opt(1:j),u_opt(1:j),1);
                sys_id(idx,:,i) = [B_id F_id];
                idx = idx + 1;
            end
        end
        reset(obj_oe)
    end
end


%% Initialize OE-Model ('full')
    function init_sys = make_init(sys)
        A = sys.A;
        B = sys.B;
        A([0 A(2:end)] ~= 0) = NaN;
        B([0 B(2:end)] ~= 0) = NaN;
        init_sys = idtf(B,A,1);
        init_sys.Structure.Numerator.Free = sys.Structure.B.Free;
        init_sys.Structure.Denominator.Free = sys.Structure.A.Free;
    end

%% Predicted Fisher Information
    function F_pred = fisher_pred(H,H_eff,u_sym,y_sym,oemodel,na,nb,theta_sym)
        % Compute the sum of the H predicted Fisher Information
        F_pred = zeros(na+nb);
        y_pred = sym('y_pred', [1 H]);
        grad = sym('grad', [na+nb H]);
        y = y_sym; % Initialization of y with past data (here: y_sym)

        for h = 1:H_eff
            u = u_sym(h:h+nb-1);
            y_pred(h) = oemodel(u,y);
            y = [y_pred(h) y(1:end-1)];
            grad(:,h) = gradient(y_pred(h),theta_sym);
            F_pred = F_pred + grad(:,h)*grad(:,h)';
        end

        % Symbolic expression to function handle (for faster computation)
        F_pred = matlabFunction(F_pred,"Vars",{theta_sym(1:end),u_sym(1:end),y_sym(1:end)});
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

%% Update Cost Function
    function I_upd = update_cost(sys,u,y)
        theta_a = sys.A;
        theta_a(theta_a ~= 0) = 1;
        theta_a = theta_a(2:end);
        theta_b = sys.B;
        theta_b(theta_b ~= 0) = 1;
        theta_b = theta_b(2:end);

        I_upd = (diag([theta_b -theta_a])*[u; y])*((diag([theta_b -theta_a])*[u; y])');
    end

end
