function [ dirac_within_support ] = check_dirac_within_smooth_interval( disc_d1, disc_d2, tt_k1, tt_k2, T_s )
%CHECK_DIRAC_WITHIN_SUPPORT_OF_KERNELS Check that the estimated Dirac is
%within the support of all kernels in the interval we use for estimation

dirac_within_support = 1;
if tt_k1 > disc_d1(1,2)*T_s || tt_k2 > disc_d1(1,2)*T_s || tt_k1 < disc_d1(1,1)*T_s|| tt_k2 < disc_d1(1,1)*T_s
    dirac_within_support = 0;
elseif tt_k1 > disc_d2(1,2)*T_s || tt_k2 > disc_d2(1,2)*T_s || tt_k1 < disc_d2(1,1)*T_s|| tt_k2 < disc_d2(1,1)*T_s
    dirac_within_support = 0;
else
    dirac_within_support = 1;
end
end

