function [ aa_k ] = retrieve_amplitudes(alpha_vec, s_m_1, s_m_2, tt_k)
%RETRIEVE_AMPLITUDES Retrieve the input Dirac amplitudes, using the signal
%moments

A = zeros(2,2);
alpha_0 = alpha_vec(1,1);
alpha_1 = alpha_vec(1,2);
A(1,1) = exp(alpha_0*tt_k(1,1));
A(1,2) = exp(alpha_0*tt_k(2,1));
A(2,1) = exp(alpha_1*tt_k(1,1));
A(2,2) = exp(alpha_1*tt_k(2,1));

S = [s_m_1,s_m_2]';

aa_k = A\S; 

end

