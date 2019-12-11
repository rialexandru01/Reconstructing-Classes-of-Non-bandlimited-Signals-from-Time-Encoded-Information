function [idx_t, idx_tphi] = range_intersection_simple(t, tphi)
%This function finds the intersection of two time intervals t and tphi,
%Assuming the interval tphi comes later compared to t.
%The output idx_t contains the indices in t which are part of the
%intersection.
%The output idx_tphi contains the indices in tphi which are part of the
%intersection.

if tphi(end,1) <= t(end,1) %If the intersection contains the whole tphi
    idx_t = find(t>=tphi(1,1) & t<=tphi(end,1));
    idx_tphi = find(tphi<=t(end,1));
elseif tphi(end,1) > t(end,1) %If tphi extends beyong t
    if tphi(1,1) <= t(end,1) %If there is some overlap between tphi and t
        idx_t = find(t>=tphi(1,1));
        idx_tphi = find(tphi<=t(end,1));
    else %If there is no overlap between tphi and t
        idx_t = [];
        idx_tphi = [];
    end
end

%The lengths of the resulting indices vectors should be the same, but in
%practice they may differ by one sample.
%Here we ensure the lengths are the same by eliminating a sample from the
%longer vector.
len_f = length(idx_t);
len_phi = length(idx_tphi); 
if len_f ~= len_phi
    if len_f > len_phi
        idx_t = idx_t(1:len_phi,1);
    elseif len_f < len_phi
        idx_tphi = idx_tphi(1:len_f,1);
    end
end

end

