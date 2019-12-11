function [  ] = plot_exp_recon(alpha_vec, T_s, A, c_m )
%PLOT_EXP_RECON Summary of this function goes here
%   Detailed explanation goes here

alpha_0 = alpha_vec(1,1);
alpha_1 = alpha_vec(1,2);
t = T_s:T_s:round(2/T_s)*T_s;
exp_0 = exp(-alpha_0*t);
exp_1 = exp(-alpha_1*t);

for m = 0:1
    exp_true = exp(-alpha_vec(1,m+1)*t);
    exp_estimated = c_m(m+1,1).*(A(1,1).*exp_0+A(2,1).*exp_1) + c_m(m+1,2).*(A(1,2).*exp_0+A(2,2).*exp_1);
    
    figure
    plot(t, real(exp_true)+0.00001, '--k');
    hold on;
    plot(t, real(exp_estimated),'r' );
    hold off;
    
    Legend = cell(1,2);
    Legend{1,1} = 'True exponential';
    Legend{1,2} = 'Reconstructed exponential';
    lgd = legend(Legend);
    lgd.FontSize = 10;
end

end

