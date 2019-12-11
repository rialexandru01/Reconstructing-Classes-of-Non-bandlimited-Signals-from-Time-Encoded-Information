function [ c_m ] = get_c_m_n( alpha_vec, t_0, t_1, t_2, T_s)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : get_c_m_n.m
% -------------------------------------------------------------------------
% Compute the coefficients needed for reproduction of exponentials
%
% USAGE:
% [ c_m ] = get_c_m_n( alpha_vec, t_0, t_1, t_2)%
% INPUT:
%  - alpha_vec : Vector of frequencies, which the sampling kernel can
%  reproduce
%  - t_i       : Location of the i output sample
%  - T_s       : Resolution of the analog signals, when represented in
%  digital form
% OUTPUT:
%  - c_m       : Coefficients that reproduce exponentials
%


alpha_0 = alpha_vec(1,1);
alpha_1 = alpha_vec(1,2);

A = zeros(2,2);
c_m = zeros(2,2);

A(1,1) = exp(alpha_0*t_1)- exp(alpha_0*t_0); 
A(1,2) = exp(alpha_0*t_2)- exp(alpha_0*t_1); 
A(2,1) = exp(alpha_1*t_1)- exp(alpha_1*t_0); 
A(2,2) = exp(alpha_1*t_2)- exp(alpha_1*t_1); 
A = A.*(1/(alpha_0*(alpha_0-alpha_1)));

B = [1 0]';
c_m (1,:) = A\B;

B = [0 1]';
c_m (2,:) = A\B;

%Plot the reproduction of exponentials
% plot_exp_recon(alpha_vec, T_s, A, c_m );
end

