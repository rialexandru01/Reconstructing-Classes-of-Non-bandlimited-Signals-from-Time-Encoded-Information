function [ ] = plot_all_signals( t_int, t_sig, x, t_y1, y1, f1, t_f1, t_y2, y2, f2, t_f2)
%PLOT_ALL_SIGNALS Plot the input signal, the filtered input
%The output non-uniform time samples, and the output kernels
figure
subplot(2,1,1)
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_f1, f1);
plot(t_y1, y1, 'b', 'Linewidth', 2);
hold off;
xlabel('t[s]');
% axis([0 t_int -1.7 1.7]);
title('Channel 1: input signal, filtered input and output of integrator');
axis([0 5 -0.5 2]);

subplot(2,1,2)
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_f2, f2);
plot(t_y2, y2, 'b', 'Linewidth', 2);
hold off;
xlabel('t[s]');
title('Channel 2: input signal, filtered input and output of integrator');
axis([0 5 -0.5 2]);

Legend = cell(1,5);
Legend{1,1} = 'Input signal';
Legend{1,2} = 'Filtered input';
Legend{1,3} = ['Integrated filtered input'];
Legend{1,4} = 'Output non-uniform time shifts';
Legend{1,5} = 'Non-uniformly shifted kernels';
lgd = legend(Legend);
lgd.FontSize = 10;


end

