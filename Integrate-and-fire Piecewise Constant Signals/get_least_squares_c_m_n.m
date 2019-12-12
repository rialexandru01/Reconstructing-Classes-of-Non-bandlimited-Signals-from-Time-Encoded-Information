function [ c_m, error_approximation] = get_least_squares_c_m_n( alpha_vec, h, t_h, t_n, T_s, L_phi,P)
%GET_LEAST_SQUARES_C_M_N 
%For some p in the right interval
%Integrate h (reversed phi) between t_n_0 and t_n_1
%t_g is chosen in a smooth interval of the modified kernel
t_g = t_n(end)-T_s*round(L_phi/(P)):T_s:t_n(1);

t_g = t_g';
no_spikes = length(t_n)-1;
length_g = length(t_g);
phi_t_n = zeros(length(t_g), no_spikes);

for m = 0:length(alpha_vec)-1
    %Exponential we want to reproduce
    g = exp(-alpha_vec(1,m+1)*t_g);
   
    for t = 1:length_g %compute the new kernel at each time t where we want to reproduce exponentials
        t_h_temp = t_h + t_g(t); %delay the new kernel by the time stored in t_g
        for n = 1: no_spikes %compute this for each delay, corresponding to the output spike time n
            t_start = find(t_h_temp == t_n(n));
            t_end = find(t_h_temp == t_n(n+1));
            phi_t_n(t,n) = T_s*trapz(h(t_start:t_end));%area between consecutive spikes, with spacing T_s
        end
    end


    B = [];
    A = [];

    for n = 1: no_spikes
        %Compute the inner product between g (the actual exponential) and
        %each shifted new kernel phi_t_n (corresponding to the area between
        %two consecutive spikes)
        phi_n = phi_t_n(:,n); %integral between spike n and spike n+1
        phi_t_n_rev = real(phi_n(end:-1:1));
        conv_g_phi_n = conv(g, phi_t_n_rev);
        
        B = vertcat(B, conv_g_phi_n);
        
        %Compute the inner products between phi_n and all the other new kernels
        A_row_n = [];
        for j = 1: no_spikes
            conv_phi_j_phi_n =  conv(phi_t_n(:,j), phi_t_n_rev);
            A_row_n = [A_row_n conv_phi_j_phi_n];
        end
        A = vertcat(A, A_row_n);
    end
    
    c_m(m+1,:) = A\B; %coefficients

    %Compute the error between the estimation and the true exponential
    exp_true = g;
    exp_estimated = c_m(m+1,1)*phi_t_n(:,1);
    for n = 2:no_spikes
        exp_estimated = exp_estimated+ c_m(m+1,n)*phi_t_n(:,n);
    end
    error_approximation = immse(real(exp_true), real(exp_estimated));

    %Plot the actual exponential, and its approximation
        
    %Plot the modified kernel.
%     figure  
%     subplot(1,2,1)
%     for n = 1: no_spikes
%         plot(t_g, phi_t_n(:,n));
%         hold on;
%     end
%     hold off;
%     
%     subplot(1,2,2)
%     plot_exp_recon(exp_true, exp_estimated, t_g);
end


end

