function [ discontinuities ] = find_discontinuities( n_vec_burst, L_phi, T_s )
%FIND_DISCONTINUITIES This function splits the interval corresponding to
%the current burst of 2 Diracs into smaller smooth intervals

discontinuities = [];
n_vec_burst = sort(n_vec_burst);
for i = 1:length(n_vec_burst)-1
    n_vec_1 = n_vec_burst(i);
	n_vec_2 = n_vec_burst(i+1);
    if (n_vec_2-n_vec_1)*T_s>L_phi/2
        r = [n_vec_2-L_phi, n_vec_1];
        discontinuities = [discontinuities; r];
    elseif (n_vec_2-n_vec_1)*T_s<L_phi/2
        r = [n_vec_2-L_phi, n_vec_1-L_phi/2];
        discontinuities = [discontinuities; r];
        r = [n_vec_1-L_phi/2, n_vec_2-L_phi/2];
        discontinuities = [discontinuities; r];
        r = [n_vec_2-L_phi/2, n_vec_1];
        discontinuities = [discontinuities; r];
    end
end


end

