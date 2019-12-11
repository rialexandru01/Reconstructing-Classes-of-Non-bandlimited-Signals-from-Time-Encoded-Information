function [y_n_current, t_n_current] = get_samples_current_dirac(b, tt_k, T_s, L_phi, t_n, y_n)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : get_samples_current_dirac.m
% -------------------------------------------------------------------------
% Retrieve the second and third output samples, occurring after the dirac
% indexed 'b'
%
% USAGE:
% [y_n_current, t_n_current] = get_samples_current_dirac(b, tt_k, T_s, L_phi, t_n, y_n);    
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
%
% OUTPUT:
%  - y_n_current   : Samples which can be used for the estimation of the
%  Dirac indexed 'b'
%  - t_n_current   : Locations of the samples y_n_current



if b>1
        last_dirac = tt_k((b-1),1);
        index_c = find(t_n>=last_dirac/T_s+L_phi); %get samples after last_dirac+L, where last_dirac is the location of the previously (last) estimated Dirac
        t_n_current = t_n(index_c); 
        
        %Select the second and third output samples, in order to ensure correct
        %sequential retrieval
        t_n_current = t_n_current(1:3);
        y_n_current = y_n(index_c(2:3));
else
        t_n_current = t_n(1:3);
        y_n_current = y_n(2:3);
end

end

