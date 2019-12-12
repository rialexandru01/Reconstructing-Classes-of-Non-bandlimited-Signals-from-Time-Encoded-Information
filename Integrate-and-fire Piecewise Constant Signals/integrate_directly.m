function [ y, t_y ] = integrate_directly( f, t_f, T_s )
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2019
%
% Date        : 31/10/2019
% Supervisor  : Professor Pier Luigi Dragotti
% Author      : Roxana Alexandru
%
% File        : integrate_directly.m
% -------------------------------------------------------------------------
% Integrate the filtered input
%
% USAGE:
% [ y, t_y ] = integrate_directly( f, t_f, T_s )
% INPUT:
%  - f         : filtered input
%  - t_f       : Support of the filtered input
%  - T_s       : Resolution of the analog signals, when represented in
%  digital form
%
% OUTPUT:
%  - y        : Integrated fltered signal
%  - t_y      : Support of integrated signal
%

t_y = t_f;

y = zeros(length(t_y)-1,1);
for i = 1:length(t_y)-1
    y(i,1) = T_s*trapz(f(1:i));
end

end

