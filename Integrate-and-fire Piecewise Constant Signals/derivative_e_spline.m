function [ dphi_dt_direct, t_dphi_dt] = derivative_e_spline( phi, t_phi, L_phi, alpha_vec, T_s)
%DERIVATIVE_E_SPLINE Summary of this function goes here
%   Detailed explanation goes here
dphi_dt_direct = diff(phi)./T_s;
t_dphi_dt = t_phi(1:end-1);

%This code is for computation of the derivative of the first-order E-spline
% alpha_0 = alpha_vec(1,1);
% alpha_1 = alpha_vec(1,2);
% 
% t_dphi_dt_1 = t_phi(1:round(L_phi/2));
% t_dphi_dt_2 = t_phi(round(L_phi/2)+1:end);
% 
% dphi_dt = zeros(length(phi),1);
% dphi_dt(1:round(L_phi/2)) = (1/(alpha_0-alpha_1))*(alpha_0*exp(alpha_1-alpha_0-alpha_0*t_dphi_dt_1)-alpha_1*exp(-alpha_1+alpha_0-alpha_1*t_dphi_dt_1));
% 
% dphi_dt(round(L_phi/2)+1:end) = (1/(alpha_0-alpha_1))*(-alpha_0*exp(-alpha_0*t_dphi_dt_2)+alpha_1*exp(-alpha_1*t_dphi_dt_2));

%Plot the E-spline and its derivative, computed in two different ways
figure
plot(t_phi, phi, 'b');
hold on;
plot(t_dphi_dt, dphi_dt_direct, 'g', 'Linewidth', 2);
% plot(t_phi, dphi_dt, 'r');
Legend = cell(1,3);
Legend{1,1} = 'First-order E-spline';
Legend{1,2} = 'Derivative of E-spline using Matlab function';
Legend{1,3} = 'Derivative of E-spline using Mathematical formula';
lgd = legend(Legend);
lgd.FontSize = 10;

end

