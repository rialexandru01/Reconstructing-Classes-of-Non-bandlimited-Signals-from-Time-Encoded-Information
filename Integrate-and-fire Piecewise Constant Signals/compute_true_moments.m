function [ s_m ] = compute_true_moments( a_k, t_k, alpha_vec)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : compute_true_moments.m
% -------------------------------------------------------------------------
% Compute the true signal moments, using the true locations and amplitudes
% of the input Diracs
%
% USAGE:
% [ s_m ] = compute_true_moments( a_k, t_k, alpha_vec)
% INPUT:
%  - a_k       : Correct amplitudes of all the Diracs
%  - t_k       : Correct time locations of all the Diracs
%  - alpha_vec : Vector of frequencies, which the sampling kernel can
%  reproduce
% OUTPUT:
%  - s_m       : Signal moments
%

s_m = zeros(2,1);

for m = 0:1
    for i = 1:length(a_k)
        s_m(m+1,1) = s_m(m+1,1) + a_k(1,i)*exp(-alpha_vec(1,m+1)*t_k(1,i));
    end
end

end

