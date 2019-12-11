function [ y, t_y ] = get_output_of_a_channel( f, t_f, a_k, t_k, omega, T_s, L_phi )
%GET_OUTPUT_OF_A_CHANNEL SThis function integrates multiple Diracs

delta_integrate = L_phi; %we integrate from tau_1 to tau_1 + delta_integrate
%Discret time of the output
t_y = t_f;
y = zeros(length(t_y), 1);
for i = 1:length(a_k)
    a_i = a_k(i);
    t_i = t_k(i);
    [ y_i, t_y_i ] = integrate_no_fire( f, t_f, delta_integrate, a_i, t_i, omega, T_s, L_phi );
    y = y + y_i;
end


end

