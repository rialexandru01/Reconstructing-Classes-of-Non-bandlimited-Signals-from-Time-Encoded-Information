function [  ] = plot_shifted_kernels_and_diracs( t_k, a_k, t_n, y_n, t_n_burst, T_s, C_T, phi, t_phi )

figure
hold on;
stem(t_k, a_k, 'kx','fill','Markersize',5);
stem(t_n*T_s, y_n./abs(y_n)*C_T, '','fill','Markersize',3);

for n = 1:2
    t_h = t_n_burst(1,n)*T_s:T_s:t_n_burst(1,n+1)*T_s;
    h_box = ones(1,length(t_h));
    phi_t_n = T_s*conv(phi, h_box);
    t_start = t_phi(1) + t_h(1);
    t_end = t_phi(end) + t_h(end);
    t_g = t_start : T_s : t_end;
if n==1
    plot(t_g, phi_t_n, 'b-', 'Linewidth', 1);
else
    plot(t_g, phi_t_n, 'r-', 'Linewidth', 1);
end
end

Legend = cell(1,2);
Legend{1,1} = 'Input Diracs';
Legend{1,2} = 'Output samples';
lgd = legend(Legend);
lgd.FontSize = 10;
end

