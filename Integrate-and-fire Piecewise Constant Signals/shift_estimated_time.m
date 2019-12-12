function [ tt_k_burst ] = shift_estimated_time( tt_k, uu_k, L_phi, T_s, lambda, t_n_0_channel1, t_n_0_channel2)
%SHIFT_ESTIMATED_TIME

angle_uu_k = mod(angle(uu_k),2*pi) ;
% angle_uu_k = angle(uu_k);
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

