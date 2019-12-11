function [  ] = plot_exp_recon(exp_true, exp_estimated, t_g)


plot(t_g, real(exp_true)+0.005, '--k');
hold on;
plot(t_g, real(exp_estimated),'r' );
hold off;


Legend = cell(1,2);
Legend{1,1} = 'True exponential';
Legend{1,2} = 'Reconstructed exponential';
lgd = legend(Legend);
lgd.FontSize = 10;

end

