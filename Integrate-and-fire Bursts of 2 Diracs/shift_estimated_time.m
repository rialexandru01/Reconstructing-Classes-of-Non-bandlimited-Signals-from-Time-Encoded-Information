function [ tt_k_burst ] = shift_estimated_time( tt_k, uu_k, L_phi, T_s, lambda, t_n_0_channel1, t_n_0_channel2)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : shift_estimated_time.m
% -------------------------------------------------------------------------
% Adjust the location of the current estimated Dirac, such that it is
% located after tau_k+L (where tau_k is the location of the previous
% estimated Dirac), and before the first output sample we used to estimate
% this Dirac (which is located at t_n_0)
%
% USAGE:
% [ tt_k_burst ] = shift_estimated_time( tt_k, uu_k, L_phi, T_s, lambda, t_n_0_channel1, t_n_0_channel2)  
% INPUT:
%  - tt_k      : Vector of estimated time locations
%  - uu_k      : Roots of the annihilating filter (Prony's method)
%  - L_phi     : Length of the support of the sampling kernel
%  - T_s       : Sampling period (to represent an analog signal in digital
%                form)
%  - lambda    : Parameters used in the definition of the frequencies
%  reproduced by the sampling kernel
%  - t_n_0_channeli     : First output sample used to estimate the current 
%                         Dirac, for channel i
% OUTPUT:
%  - tt_k_burst       : Estimation of the current Dirac in the stream
%

angle_uu_k = mod(angle(uu_k),2*pi) ;
if ~isempty(tt_k)
    last_dirac = tt_k(end);
    
    %Find shift c_min to ensure estimated time is after the start time of
    %the interval (at a minimum length L_phi from previous burst)
    start_time = last_dirac + L_phi*T_s;
    c_min = ceil((start_time*lambda - angle_uu_k)/(2*pi));
    c_min = min(c_min);
    
    %Find maximum shift c_max to ensure estimated time is before the first
    %output samples of this burst (t_n_0_channel1 and t_n_0_channel2)
    end_time = min(t_n_0_channel1, t_n_0_channel2);
    c_max = floor((end_time*lambda - angle_uu_k)/(2*pi));
    c_max = max(c_max);
    
    %Shift the estimated time
    if c_max >= c_min
        tt_k_burst = sort((angle_uu_k + 2*pi*c_min)/lambda);
    end

else
    tt_k_burst = sort(angle_uu_k/lambda);
end

