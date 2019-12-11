function [s_m, n_vec_d, y_n_d, disc_d, c_m_n] = compute_channel_moments(n_vec_burst, L_phi, T_s, y_n_burst, P, alpha_vec, phi, t_phi)
%COMPUTE_CHANNEL_MOMENTS This function computes the moments for a single
%input channel

%First 3 output samples
t_n = n_vec_burst(1,1);
t_n_1 = n_vec_burst(1,2);
t_n_2 = n_vec_burst(1,3);

%Find all discontinuities containing the non-zero samples corresponding
%to the current burst
discontinuities1 = find_discontinuities( n_vec_burst, L_phi, T_s );

%Compute the signal moments for first channel
disc_d = [t_n_2 - round(L_phi/2), t_n];
n_vec_d = [n_vec_burst(1,2), n_vec_burst(1,3)];
y_n_d = [y_n_burst(2,1), y_n_burst(3,1)]';
%Get the coefficients from the system of equations
no_kernels = 2;
[c_m_n, t_range] = get_c_m_tn_exp_from_system_eq(P, alpha_vec, disc_d, n_vec_d, no_kernels, phi, t_phi);
s_m = c_m_n * y_n_d;
end

