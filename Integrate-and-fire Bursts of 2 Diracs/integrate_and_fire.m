function [ t_n, y_n, y_with_reset] = integrate_and_fire( y, C_T )
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 11/12/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : integrate_and_fire.m
% -------------------------------------------------------------------------
% Reset the integrated filtered input, whenever the trigger mark C_T is
% reached
%
% USAGE:
% [ t_n, y_n ] = integrate_and_fire( y, C_T )
% INPUT:
%  - y         : Integrated signal
%  - C_T       : Trigger mark of the threshold detector
%
% OUTPUT:
%  - y_n      : Amplitude of output samples (close to C_T)
%  - t_n      : Location of the output spikes
%  - y        : The integrate filtered input, with spike-triggered resets


t_n = [];
y_n = [];
finished = 0;
t_n_last = 0;

while finished == 0
    t_n_new = find(y(t_n_last+1:end)>=C_T | y(t_n_last+1:end)<=-C_T,1,'first')+t_n_last;
    t_n = [t_n, t_n_new]; 
    y_n = [y_n, y(t_n_new)];
    if ~isempty(t_n_new)
        y(t_n_new+1:end) = y(t_n_new+1:end) - y(t_n_new);
    else
        finished = 1;
    end
    t_n_last = t_n_new;
end

y_with_reset = y;

