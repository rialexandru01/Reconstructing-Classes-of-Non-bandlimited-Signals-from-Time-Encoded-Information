function [ cum_error_shifts, cum_error_samples, y_n_est1, y_n_est2 ] = check_output_samples_using_estimated_diracs( ind_n_d1, ind_n_d2, n_vec_d1, n_vec_d2, y_n_d1, y_n_d2, tt_k1, tt_k2, aa_k1, aa_k2, T_s, t_sig, h1, t_h1, h2, t_h2, A, f_tem)
%CHECK_OUTPUT_SAMPLES_USING_ESTIMATED_DIRACS 
%Compute the output samples using the estimated Diracs, and the same
%acquisition model, based on the comparator
%Compute the error between the estimated output samples and the actual
%output samples


cum_error_shifts = Inf;
cum_error_samples = Inf;
y_n_est1 = [];
y_n_est2 = [];

if ~isempty(tt_k1) && tt_k1>0 && ~isempty(tt_k2) && tt_k2>0 
    %Estimated input
    x_est = zeros(size(t_sig));
    t_est1_index = round(tt_k1 / T_s);
    x_est(t_est1_index) = aa_k1;

    t_est2_index = round(tt_k2 / T_s);
    x_est(t_est2_index) = aa_k2;

    %Estimated filtered input
    [y_est1, t_y1] = filter_input( x_est, h1, t_sig, t_h1, T_s );
    [y_est2, t_y2] = filter_input( x_est, h2, t_sig, t_h2, T_s );

    n_vec_est1  = compare_filtered_input( A, f_tem, y_est1, t_y1, T_s );
    n_vec_est2  = compare_filtered_input( A, f_tem, y_est2, t_y2, T_s );

    %output values at the non-uniform time shifts
    y_n_est1 = y_est1(n_vec_est1');
    y_n_est2 = y_est2(n_vec_est2');

    %Output samples of estimated signal, at the actual non-uniform
    %times corresponding to the current burst
    n_vec_d1_est = n_vec_est1(ind_n_d1);
    n_vec_d2_est = n_vec_est2(ind_n_d2);

    y_n_d1_est = y_n_est1(ind_n_d1);
    y_n_d2_est = y_n_est2(ind_n_d2);

    error_shifts1 = sum(((n_vec_d1_est - n_vec_d1)*T_s).^2);
    error_samples1 = sum(((y_n_d1_est - y_n_d1)*T_s).^2);

    error_shifts2 = sum(((n_vec_d2_est - n_vec_d2)*T_s).^2);
    error_samples2 = sum(((y_n_d2_est - y_n_d2)*T_s).^2);

    cum_error_shifts = error_shifts1+error_shifts2;
    cum_error_samples = error_samples1 + error_samples2;

end


end

