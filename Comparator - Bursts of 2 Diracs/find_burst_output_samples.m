function [ y_n_burst, n_vec_burst, y_n_temp ] = find_burst_output_samples( y_n_temp, y_n, n_vec )
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 12/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : find_burst_output_samples.m
% -------------------------------------------------------------------------
% Find the output samples corresponding to the
%current burst
%
% USAGE:
% [ y_n_burst, n_vec_burst, y_n_temp ] = find_burst_output_samples( y_n_temp, y_n, n_vec )
% INPUT:
%  - y_n_temp                       : Updated output samples (set to zero up to current burst)
%  - y_n                            : Vector of unaltered output samples
%  - n_vec                          : Vector of locations of output samples                 
%
% OUTPUT:
%  - y_n_burst                      : Samples used to reconstruct the Diracs in the current burst
%  - n_vec_burst                    : Time locations of samples used to reconstruct the Diracs in the current burst
%  - y_n_temp                       : Updated output samples 
%

%Find the next burst of non-zero output samples
y_n_start = find(y_n_temp ~= 0, 1, 'first');

%Find the next zero sample after y_n_start (i.e. the point where the
%samples of the current burst end)
y_n_stop = y_n_start + find(y_n_temp(y_n_start:end) == 0, 1, 'first')-2;
y_n_burst = y_n(y_n_start:y_n_stop);
n_vec_burst = n_vec(y_n_start:y_n_stop);

%Set the output samples already used for estimation to zero
y_n_temp(1:y_n_stop) = 0;
end

