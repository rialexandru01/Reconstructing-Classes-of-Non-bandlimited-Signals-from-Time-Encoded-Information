function [ y, t_y ] = filter_input( x, h, t_sig, t_h, T_s )
%FILTER_INPUT This function takes as input the continuous-time signal x(t)
%and outputs the filtered input signal, where the sampling kernel is h(t)=phi(-t)

y = conv(x, h);

t_0 = t_sig(1) + t_h(1);
t_f = t_sig(end) + t_h(end);
t_y = t_0 : T_s : t_f;

end

