function [ ] = plot_all_signals( t_sig, x, t_y, y, f, t_f, t_int)
%PLOT_ALL_SIGNALS Plot the input signal, the filtered input
%The output non-uniform time samples, and the output kernels

figure
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_f(1:end-1), f);
plot(t_y(1:end-1), y, 'b', 'Linewidth', 2);
% scatter(idx*T_s, zeros(1,length(idx)), 'or', 'filled');
%     scatter(idx+L_phi, zeros(1,length(idx)), 'or', 'filled');

hold off;

axis([0 t_int -2 3]);
xlabel('t[s]');

Legend = cell(1,5);
Legend{1,1} = 'Input signal';
Legend{1,2} = 'Filtered input';
Legend{1,3} = ['Integrated filtered input'];
Legend{1,4} = 'Output non-uniform time shifts';
Legend{1,5} = 'Non-uniformly shifted kernels';
lgd = legend(Legend);
lgd.FontSize = 10;


end

