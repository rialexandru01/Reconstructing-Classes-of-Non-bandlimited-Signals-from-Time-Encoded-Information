function [ x_noisy, e_n] = add_noise( x, sigma )
%ADD_NOISE Adds noise to the input signal x

e_n  = sigma .* randn(length(x),1);
if var(e_n)>sigma
    e_n  = sigma .* randn(length(x),1); %randn returns normally distributed random numbers
end
x_noisy = x + e_n';

P_x = x * x' / length(x);
P_n = e_n'*e_n /length(e_n);
noise_var = 10.* log10(P_x/P_n);
end

