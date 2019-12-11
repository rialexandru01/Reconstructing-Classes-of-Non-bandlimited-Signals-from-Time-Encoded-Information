figure
set(gcf, 'Position', [50 200 420 315])
stem(t_k(1:K), a_k(1:K), '^k', 'fill', 'LineWidth', 1, 'MarkerSize', 10)
axis([t_y(1) t_y(end) -.5 4])
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)

figure
K_new = 1;
set(gcf, 'Position', [500 50 420 315])
stem(tt_k(1:K_new), aa_k(1:K_new), '^k', 'fill', 'LineWidth', 1, 'MarkerSize', 10)
axis([t_y(1) t_y(end) -.5 4])
hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
set(hdl, 'FontSize', label_size)
% hdl = title('Reconstructed stream of Diracs', 'Interpreter', 'Latex');
% set(hdl, 'FontSize', label_size)
set(gca, 'FontSize', axis_size)
