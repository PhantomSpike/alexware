function dexp_fit = fitdexp(k_h,input_params)
%dexp = fitdexp(k_h,time_bin_ms,Ainit_coeff,Binit_coeff,alphainit_coeff,betainit_coeff,method)
%This function fits a double exponential fit to the temporal kernel part of the strf - k_h
%The fit is done using minFunc with MSE loss term
%The double exponential has the form:
%k_hat(h) = A*exp(-h/alpha) - B*exp(-h/beta). A fast decaying positive exp
%and a slower decaying negative exp
%
% input params:
%   k_h -- the temporal kernel part of the strf k_h
%   input_params -- struct containing the input parameters used for
%   initialization of the fits
%       .A_init - the initialization val of the A param, the scalar of the
%           +ve exp
%       .B_init - the initialization val of the B param, the scalar of the
%           -ve exp
%       .alpha_init - the initialization val of the alpha param, the
%           time const of the +ve exp
%       .beta_init - the initialization val of the beta param, the
%           time const of the -ve exp
%       .init_method - the initialization method used for the parameters
%           'exp' - the vals are drawn from exponential distribution
%           'gaussian' - the vals are drawn from gaussian distribution 
%       time_bin_ms - the time bin used for the strf in ms


%
% output params:
%   dexp_fit -- struct containing the results of the fit
%       .A - the scalar A of the +ve exp
%       .B - the scalar B of the -ve exp
%       .alpha - the time const alpha of the +ve exp
%       .beta - the time const beta of the -ve exp
%       .k_h_hat - the best prediction k_h kernel fit (lowest MSE)
%
% Author: Aleksandar Ivanov
% Year: 2019
%% Params - Some other parameters used to run the fits
epsilon_AB = 1e-2; %This will be a small value we will use for the lower bound of A and B
lb_t_const_ms = 0.5; %This is the lowest sensible values for alba and beta
n_iter = 100; %The number of iterations to run before selecting
sigma = 1; %The standard deviation of the initialization
check_grad = 0; %Logical flag whether to check the gradient
interpolate = input_params.interpolate; %Logical flag whether to do interpolation
upsamp_factor = 100;
method = 'makima';
n_params = 4;
%% Unroll input parameters
time_bin_ms = input_params.time_bin_ms;
h = (0:1:numel(k_h)-1)'; %The history steps of the kernel in steps
h = time_bin_ms*h; %Convert to ms

if interpolate
    [k_h,h] = interpolate_kh(k_h,h,upsamp_factor,method);
end

h_max = max(h); %Find the max history step

A_init = input_params.A_init;
B_init = input_params.B_init;

if ~isempty(input_params.alpha_init)
    alpha_init = input_params.alpha_init;
else
    alpha_init = 20; %Initialize alpha to some sensible short value as we know that excitation is usually much shorter
end

if ~isempty(input_params.beta_init)
    beta_init = input_params.beta_init;
else
    beta_init = h_max; %Initialize beta to with the max time step so that we can try a range a of values later
end

if ~isempty(input_params.ub_A_B)
    ub_A_B = input_params.ub_A_B;
else
    ub_A_B = 3; %Initialize beta to with the max time step so that we can try a range a of values later
end




%% Specify params for fmincon
lb_A = epsilon_AB; lb_B = epsilon_AB; lb_alpha = lb_t_const_ms; lb_beta = lb_t_const_ms; % Define the lower bound for the fmincon. We want only positive values so 0
ub_A = ub_A_B; ub_B = ub_A_B; ub_alpha = h_max; ub_beta = h_max; %Define the upper bound for the fmincon. We want some sensible values for all params w/o introducing a bias
LB = [lb_A,lb_B,lb_alpha,lb_beta]; %Concatenate together
UB = [ub_A,ub_B,ub_alpha,ub_beta];
A_f = []; %These are other constraints that we don't need
b_f = [];
Aeq = [];
beq = [];
nonlcon = [];
opts = optimoptions('fmincon','SpecifyObjectiveGradient',true,'display','off'); %Here we specify some extra setting for fmincon

%% Run fmincon using many initializations with gaussian noise
vars_res = zeros(n_params,n_iter); %Here we will store the variables we get out of the different iterations
Loss = zeros(n_iter,1); %Here we store the MSE for every iteration

for n = 1:n_iter
    A = min(max(A_init + randn(1,1)*sigma,epsilon_AB),ub_A_B);
    B = min(max(B_init + randn(1,1)*sigma,epsilon_AB),ub_A_B);
    alpha = max(alpha_init*rand(1,1),lb_t_const_ms);
    beta = max(beta_init*rand(1,1),lb_t_const_ms);
    vars = [A,B,alpha,beta];
    if check_grad
        for j = 1:10
            d = checkgrad('double_exp_E_dE_dp', rand(n_params,1), 1e-5,rand(size(h)),rand(size(k_h)));
        end
        keyboard;
    end
    func = @(vars)double_exp_E_dE_dp(vars,h,k_h);
    [vars_res(:,n),Loss(n)] = fmincon(func,vars,A_f,b_f,Aeq,beq,LB,UB,nonlcon,opts);
end


%% Make best model and save results

%Find the fit with the lowest MSE
[best_Loss,best_fit_ix] = min(Loss);

A_r = vars_res(1,best_fit_ix);
B_r = vars_res(2,best_fit_ix);
alpha_r = vars_res(3,best_fit_ix);
beta_r = vars_res(4,best_fit_ix);

%Make the prediction k_h
exp_pos = A_r*exp(-h/alpha_r);
exp_neg = -B_r*exp(-h/beta_r);
k_h_hat = exp_pos + exp_neg;


%Save the final results to a struct
dexp_fit.A = A_r;
dexp_fit.B = B_r;
dexp_fit.alpha = alpha_r;
dexp_fit.beta = beta_r;

dexp_fit.k_h_hat = k_h_hat;
dexp_fit.k_h = k_h;
dexp_fit.h = h;
dexp_fit.exp_pos = exp_pos;
dexp_fit.exp_neg = exp_neg;

dexp_fit.all_fits = vars_res;
dexp_fit.best_Loss = best_Loss;
dexp_fit.best_fit_ix = best_fit_ix;
dexp_fit.Loss = Loss;
dexp_fit.units_t_const = 'ms';
dexp_fit.units_AB = 'AU';
dexp_fit.n_iter = n_iter;
dexp_fit.LB = LB;
dexp_fit.UB = UB;
