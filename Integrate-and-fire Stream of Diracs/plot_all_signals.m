function [ ] = plot_all_signals( t_int,  t_y, y, f, t_f)
%PLOT_ALL_SIGNALS Plot the input signal, the filtered input
%The output non-uniform time samples, and the output kernels

figure
hold on;
plot(t_f, f, 'b');
plot(t_y, y, 'r', 'Linewidth', 2);
hold off;
xlabel('t[s]');
axis([0 t_int -1 1]);

Legend = cell(1,2);
Legend{1,1} = 'Input signal';
Legend{1,2} = ['Output of integrator'];
lgd = legend(Legend);
lgd.FontSize = 10;
% axis([0 t_int -1.2 1.2]);

end

