function [y_n_burst, t_n_burst] = get_samples_burst(b, tt_k, T_s, min_separation, t_n, y_n, s, C_T)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : get_samples_burst.m
% -------------------------------------------------------------------------
% Retrieve the second and third output samples, occurring after the dirac
% indexed 'b'
%
% USAGE:
% [y_n_burst, t_n_burst] = get_samples_burst(b, tt_k, T_s, L_phi, t_n, y_n);    
%
% INPUT:
%  - b       : Index of current Dirac (e.g. 1 for the first Dirac in the
%  stream)
%  - tt_k    : Stream of estimated Diracs
%  - T_s     : Time resolution (when representing an analog signal in digital
%  form)
%  - L_phi   : Length of the support of the sampling kernel
%  - t_n     : Locations of output spikes
%  - y_n     : Value of output spikes (ideally equal to the trigger mark
%  C_T, however, due to the finite time resolution of the signal dictated by T_s, may
%  slightly differ from C_T
%  - C_T     : Value of the trigger mark of the integrate-and-fire TEM
%  - s       : (s-1) is the number of samples we discard in each burst
%
% OUTPUT:
%  - y_n_burst   : Samples which can be used for the estimation of the
%  Dirac burst indexed 'b'
%  - t_n_burst   : Locations of the samples y_n_burst

if b>1
        last_dirac = tt_k(2*(b-1),1);
        y_n_b = union(find(y_n>=C_T*0.1), find(y_n<=C_T*0.1));
        index_b = intersect(y_n_b,find(t_n>=last_dirac/T_s+min_separation));        
        t_n_burst = t_n(index_b);

        t_n_burst = t_n_burst(s:s+2);
        %Select the second and third output samples, in order to ensure correct
        %sequential retrieval
        y_n_burst = y_n(index_b(s+1:s+2));

else
%For the first burst, the first sample s is not influences by previous
%bursts, and hence can be used.
        t_n_burst = t_n(s-1:s+1);
        y_n_burst = y_n(s:s+1);
end

