function [ n_vec ] = compare_filtered_input( A, f_tem, y, t_y, T_s )
%COMPARE_FILTERED_INPUT This function outputs the non-uniform samples,
%obtained by comparing the filtered input y with a comparator's reference
%signal (a cosine function)

comp = y-A*cos(f_tem*2*pi*t_y)';
tem = ZeroX(t_y,comp);
n_vec = round(tem/T_s);
t_tem = n_vec *T_s;
[~,idx] = intersect(t_y, t_tem);
n_vec = idx';
end

