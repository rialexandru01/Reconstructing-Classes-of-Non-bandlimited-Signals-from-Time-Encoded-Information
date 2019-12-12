function [ nonuniform_times_c, y_n_c ] = get_samples_burst( c, T_s, L_phi, y_n, t_n, P, tt_k, K_B)
%GET_SAMPLES_BURST Get a new set of output samples, for the next burst to
%estimate
if c>1        
    last_dirac = tt_k(K_B*(c-1));
    index_c = find(t_n>=last_dirac/T_s+L_phi);
    nonuniform_times_c = t_n(index_c);
    nonuniform_times_c = nonuniform_times_c(1:P+1);
    
    %Select the second and third output samples, in order to ensure correct
    %sequential retrieval
    y_n_c = y_n(index_c(2:P+1));
else
    nonuniform_times_c = t_n(1:P+1);
    y_n_c = y_n(2:P+1);
end

nonuniform_times_c_in_sec = nonuniform_times_c*T_s;
end

