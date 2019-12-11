function [ y_n ] = remove_dirac_contribution( tt_k_burst, aa_k_burst, T_s, t_sig, y_n, t_n, h, t_h)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : remove_dirac_contribution.m
% -------------------------------------------------------------------------
% Remove the contribution of the last estimated burst of Diracs, from the
% output samples
%
% USAGE:
% [ y_n ] = remove_dirac_contribution( tt_k_burst, aa_k_burst, T_s, t_sig, y_n, t_n, h, t_h)%
% INPUT:
%  - tt_k_burst      : Vector of estimated time locations of the Diracs in
%  the current burst
%  - aa_k_burst      : Vector of estimated amplitudes of the Diracs in the
%  current burst
%  - T_s             : Resolution of the analog signals, when represented in
%  digital form
%  - t_sig           : Support of the input signal
%  - y_n             : Vector of amplitudes of output spikes
%  - t_n             : Vector of locations of output spikes
%  - h               : Impulse response of filter
%  - t_h             : Support of the filter h(t)
% OUTPUT:
%  - y_n             : Modified output samples 
%


%Filter the estimated Dirac with the same filter
itk_burst = round(tt_k_burst/T_s);
x_burst = zeros(size(t_sig));
x_burst(itk_burst) = aa_k_burst;
f1_burst = conv(x_burst, h);
t_start = t_sig(1) + t_h(1);
t_end = t_sig(end) + t_h(end);
t_f1_burst = t_start : T_s : t_end;

%Set the first sample to 0 since this sample was used for sure in
%estimation
y_n(1) = 0;

%Remove the contribution of the current estimated Dirac from the output
%samples
 for i = 2:length(t_n)
     t_start = find(t_f1_burst==t_n(i-1)*T_s);
     t_stop = find(t_f1_burst == t_n(i)*T_s);
     y_n(i) = y_n(i) - T_s*trapz(f1_burst(t_start:t_stop)); %subtract the contribution of the estimated Dirac between the two output spikes at t_start and t_stop
 end


 
end

