function [E,dE_dp] = double_exp_E_dE_dp(vars,h,k_h)
%[E,dE_dp] = double_exp_E_dE_dp(params,h,k_h)
%This function gives back the mean squared error/loss (MSE) associated with
%the dexp model prediction. It also provides the gradient w.r.t. to the
%four parameters in the model
%
%The double exponential has the form:
%k_hat(h) = A*exp(-h/alpha) - B*exp(-h/beta). A fast decaying positive exp
%and a slower decaying negative exp
%
% input params:
%   params -- struct containing the input parameters used for the model
%
%       params(1) = A - the val of the A param, the scalar of the
%           +ve exp
%       params(2) = B - the val of the B param, the scalar of the
%           -ve exp
%       params(3) = alpha - theval of the alpha param, the
%           time const of the +ve exp
%       params(4) = beta - the val of the beta param, the
%           time const of the -ve exp
%   h -- the vector of history steps
%   k_h -- the temporal k_h kernel which we are trying to fit

%
% output params:
%   E -- the MSE for these particular parameters
%   dE_dp -- a vector containing the partial derivatives of the MSE w.r.t. to the 4 parameters
%   in the model:
%       dE_dA
%       dE_dB
%       dE_alpha
%       dE_beta
%
% Author: Aleksandar Ivanov
% Year: 2019

plot_on = 0;

%Get the parameters out making sure to take abs as these values cannot be
%-ve e.g. because of the initialization
A = vars(1);
B = vars(2);
alpha = vars(3);
beta = vars(4);

h = h(:); k_h = k_h(:); %Force col vectors
h_max = length(h); %Find the number of data points in the model

% ---------  Compute Error -----------------
exp_pos = exp(-h/alpha); %The positive exp part of the fit 
exp_neg = exp(-h/beta); %The negative exp part of the fit
k_h_fit = A*exp_pos - B*exp_neg; %The full kernel fit (prediction)
D = k_h_fit-k_h; %Find the difference between the two
E = (1/h_max).*(D'*D); %Compute the MSE

% ---------  Compute Gradient (Partial derivatives) -----------------
if (nargout > 1)
    dE_dA = (2/h_max)*(D'*exp_pos);
    dE_dB = (-2/h_max)*(D'*exp_neg);
    dE_dalpha = (2/h_max)*(A/alpha^2)*(D'*(exp_pos.*h));
    dE_dbeta = (-2/h_max)*(B/beta^2)*(D'*(exp_neg.*h));
    
    dE_dp = [dE_dA; dE_dB; dE_dalpha; dE_dbeta]; %Combine the derivatives w.r.t. all the params together
end

% ---------  Compute Hessian -------------------
% if nargout > 2
% syms A1 B1 alpha1 beta1
% f = A1*exp(-a) + 2*z*x;
% hessian(f,[x,y,z])
% end

%Optionally plot the fits as they happen
if plot_on
    plot(k_h); hold on; plot(k_h_fit);hold off;drawnow;
end


