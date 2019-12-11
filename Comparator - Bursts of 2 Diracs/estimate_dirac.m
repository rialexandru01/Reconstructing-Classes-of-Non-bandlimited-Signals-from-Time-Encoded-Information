function [ tt_k ] = estimate_dirac( interval_start_time, interval_stop_time, angle_uu_k, lambda )
%ESTIMATE_DIRAC Estimate the location and amplitude of the Dirac, using the
%current continuous interval

if(angle_uu_k/lambda > interval_stop_time)
    %Find the smallest integer, such that the inter-Dirac separation
    %constraint is satisfied (larger than L)
    r = (angle_uu_k - lambda * interval_stop_time)/(2*pi); %the estimated location tt_k should be larger than the start of the support of all kernels in the interval
    r = ceil(r); 
    %Estimated Dirac
    tt_k    = sort((angle_uu_k - 2*r*pi) / lambda); %estimated location
elseif (angle_uu_k/lambda < interval_start_time)
    %Find the smallest integer, such that the inter-Dirac separation
    %constraint is satisfied (larger than L)
    r = ( lambda * interval_start_time -  angle_uu_k)/(2*pi); %the estimated location tt_k should be larger than the start of the support of all kernels in the interval
    r = ceil(r); 
    %Estimated Dirac
    tt_k    = sort((angle_uu_k + 2*r*pi) / lambda); %estimated location
else
    tt_k    = sort(angle_uu_k / lambda); %estimated location
end

end

