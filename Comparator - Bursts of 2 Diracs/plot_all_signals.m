function [ ] = plot_all_signals( t_sig, x, t_y1, y1, t_y2, y2, A, f_tem, idx1, idx2, T_s, t_tem1, phi1, t_tem2, phi2, L_phi1, L_phi2)
%PLOT_ALL_SIGNALS Plot the input signal, the filtered input
%The output non-uniform time samples, and the output kernels
figure
subplot(2,1,1)
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_y1, y1, 'b', 'Linewidth', 2);
plot(t_y1, A*cos(f_tem*2*pi*t_y1)',  'r', 'Linewidth', 1);
scatter(idx1*T_s, zeros(1,length(idx1)), 'or', 'filled');
for i = 1:length(t_tem1)
    t = t_tem1(1,i)-(L_phi1-1)*T_s:T_s:t_tem1(1,i);
    plot(t, phi1, '--k', 'Linewidth', 0.5);
end
%grid on;
hold off;
axis([0 8 -1.5 1.5]);
xlabel('t[s]');

Legend = cell(1,5);
Legend{1,1} = 'Input signal';
Legend{1,2} = 'Filtered input';
Legend{1,3} = ['Comparator', '''', 's reference signal'];
Legend{1,4} = 'Output non-uniform time shifts';
Legend{1,5} = 'Non-uniformly shifted kernels';
ldg = legend(Legend);
lgd.FontSize = 10;


subplot(2,1,2)
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_y2, y2, 'b', 'Linewidth', 2);
plot(t_y2, A*cos(f_tem*2*pi*t_y2)',  'r', 'Linewidth', 1);
scatter(idx2*T_s, zeros(1,length(idx2)), 'or', 'filled');
for i = 1:length(t_tem2)
    t = t_tem2(1,i)-(L_phi2-1)*T_s:T_s:t_tem2(1,i);
    plot(t, phi2, '--k', 'Linewidth', 0.5);
end
%grid on;
hold off;
axis([0 8 -1.5 1.5]);
xlabel('t[s]');


end

