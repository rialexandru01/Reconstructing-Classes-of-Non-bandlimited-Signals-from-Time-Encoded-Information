function [s_m, n_vec_d, y_n_d, disc_d, c_m_n] = compute_channel_moments(n_vec_burst, L_phi, y_n_burst, P, alpha_vec, phi, t_phi)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : compute_true_moments.m
% -------------------------------------------------------------------------
% Compute the signal moments
%
% USAGE:
% [s_m, n_vec_d, y_n_d, disc_d, c_m_n] = compute_channel_moments(n_vec_burst, L_phi, y_n_burst, P, alpha_vec, phi, t_phi)
% INPUT:
%  - n_vec_burst     : Locations of output samples
%  - L_phi           : Length of the support of the sampling kernel
%  - y_n_burst       : Output samples
%  - P               : Order of the sampling kernel (exponential spline)
%  - alpha_vec       : Vector of frequencies, which the sampling kernel can reproduce
%  - phi             : Sampling kernel
%  - t_phi           : Support of the sampling kernel
% OUTPUT:
%  - s_m             : Signal moments
%  - n_vec_d         : Locations of output samples used in reconstruction of current burst
%  - y_n_d           : Output samples used in reconstruction of current burst
%  - disc_d          : Continuous interval defined by the current output samples
%  - c_m_n           : Coefficients used to reproduce exponentials in the continuous interval
%

%First 3 output samples
t_n = n_vec_burst(1,1);
t_n_2 = n_vec_burst(1,3);

disc_d = [t_n_2 - round(L_phi/2), t_n];
n_vec_d = [n_vec_burst(1,2), n_vec_burst(1,3)];
y_n_d = [y_n_burst(2,1), y_n_burst(3,1)]';

%Get the coefficients from the system of equations
no_kernels = 2;
c_m_n = get_c_m_tn_exp_from_system_eq(P, alpha_vec, disc_d, n_vec_d, no_kernels, phi, t_phi);
s_m = c_m_n * y_n_d;
end

