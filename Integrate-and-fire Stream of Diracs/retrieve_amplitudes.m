function [ aa_k_current ] = retrieve_amplitudes(alpha_vec, s_m_1, tt_k_current)
%RETRIEVE_AMPLITUDES Retrieve the amplitude of the Dirac, using the signal
%moments

alpha_0 = alpha_vec(1,1);
A = exp(alpha_0*tt_k_current(1,1));
S = s_m_1;

aa_k_current = S/A; 

end

