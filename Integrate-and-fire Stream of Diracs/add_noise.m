function [ x_noisy, e_n, sigma] = add_noise( x, SNR, length_x )
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 12/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : add_noise.m
% -------------------------------------------------------------------------
% Add noise to the input signal based on SNR 
%
% USAGE:
% [ x_noisy, e_n, sigma] = add_noise( x, SNR, length_x )
% INPUT:
%  - x              : Input signal
%  - SNR            : Signal to noise ratio
%  - length_x       : Length of input signal

% OUTPUT:
%  - x_noisy        : Input signal corrupted by noise
%  - e_n            : White Gaussian noise added to input
%  - sigma          : Variance of the white additive Gaussian noise
%

P_x = x * x'/length_x;
% P_x = x * x';
% Standard deviation of the noise for each realisation and noise matrix
sigma = sqrt(P_x * 10.^(-SNR/10));
%sigma = inamp(1)*sqrt(10.^(-SNR/10)); % PSNR
e_n  = sigma .* randn(length(x),1);
if var(e_n)>sigma
    e_n  = sigma .* randn(length(x),1); %randn returns normally distributed random numbers
end
x_noisy = x + e_n';

P_e_n = e_n'*e_n/length(e_n);
snr_test = 10*log10(P_x/P_e_n);
end

