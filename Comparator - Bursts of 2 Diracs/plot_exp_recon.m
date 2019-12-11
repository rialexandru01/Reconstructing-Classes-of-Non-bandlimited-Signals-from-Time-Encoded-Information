function [  ] = plot_exp_recon( P, alpha_vec, phi, t_phi, disc_d, n_vec_d, n_vec, T, T_s, c_m_n, t_k, a_k)
%PLOT_EXP_RECON Summary of this function goes here
%   Detailed explanation goes here

n_vec = sort(n_vec);
t1 = t_phi(1);
t2 = n_vec(end) * T_s * T + t_phi(end);
t     = (t1:T_s:t2)';

% Font sizes
label_size = 30;
axis_size  = 26;

for ith_m = 1:P+1
    alpha = alpha_vec(ith_m);   
    % Original exponential
    f = exp(alpha * t / T);
    
    % Estimated exponential
    f_rec = zeros(length(t), length(n_vec_d)+1);
    for i = 1:length(n_vec_d)
        %idx_phi are the indices in t_phi+n_vec(i)*T
        inc = n_vec_d(1, i)*T_s*T;
        tphi = t_phi + inc;
        %Find intersection
        %[~, idx_f, idx_phi] = intersect(t, tphi);
        %[idx_f, idx_phi] = range_intersection_simple(t, tphi);
        [~,idx_f, idx_phi] = intersect(t, tphi);
        f_rec(idx_f,i+1)    = c_m_n(ith_m, i) * phi(idx_phi);
        f_rec(:,1) = f_rec(:,1) + f_rec(:,i+1);
    end
    
    %Unused estimated exponential
    n_vec_unused = setdiff(n_vec, n_vec_d);
    f_rec_unused = zeros(length(t), length(n_vec_unused));
    for i = 1:length(n_vec_unused)
        %idx_phi are the indices in t_phi+n_vec(i)*T
        tphi = t_phi + n_vec_unused(1, i)*T_s*T;
        %Find intersection
        [~,idx_f, idx_phi] = intersect(t, tphi);
        f_rec_unused(idx_f,i)    = phi(idx_phi);
    end
    
    
    % Dynamic range of the reproduced polynomial
    dyn_range = max(real(f_rec(:,1))) - min(real(f_rec(:,1)));
    
    % Measure the error
    phi_support = length(phi);
    idx         = (1+phi_support:length(f)-phi_support)';
    rec_er      = f(idx) - f_rec(idx,1);
    MSE         = (rec_er' * rec_er) / length(rec_er);
    disp(['MSE = ' num2str(MSE)])

    figure
    set(gcf, 'Position', [50+100*ith_m 50 420 315])
    plot(t, real(f_rec(:,1)), 'k', 'LineWidth', 4)
    hold on
    
    plot(t, f, '--k', 'LineWidth', 2);
    
    for i = 1:length(n_vec_unused)
        plot(t, real(f_rec_unused(:,i)), '--r')
    end

    
    for i = 1:length(n_vec_d)
        plot(t, real(f_rec(:,i+1)), '-r')
    end

    hdl = xlabel('$$t$$ [s]', 'Interpreter', 'Latex');
    set(hdl, 'FontSize', label_size)
    set(gca, 'FontSize', axis_size)
    axis([t(1) t(end) min(real(f_rec(:,1)))-1*dyn_range max(real(f_rec(:,1)))+1.6*dyn_range])
    
    axis([0 8 -2 2]);
    scatter(disc_d*T_s*T, zeros(1, length(disc_d)), 'o', 'filled');
    
    stem(t_k, a_k, 'g');
    %grid on;
    hold off;
end

Legend = cell(1,3);
Legend{1,1} = 'Reconstructed exponential';
Legend{1,3} = 'Output non-uniform splines';
Legend{1,2} = 'True exponential';
lgd = legend(Legend);
lgd.FontSize = 10;

end

