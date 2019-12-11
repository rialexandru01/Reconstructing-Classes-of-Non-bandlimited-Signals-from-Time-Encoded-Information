function [ c_m_n ] = get_c_m_tn_directly( alpha_0, alpha_1, t_1, t_2, disc_d, T_s)
%GET_C_M_TN_DIRECTLY This function computes the coefficients used in the
%reproduction of two exponentials
c_m_n = zeros(2,2);

a1 = 1/(alpha_0-alpha_1)*exp(alpha_0*t_1);
b1 = -1/(alpha_0-alpha_1)*exp(alpha_1*t_1);
a2 = 1/(alpha_0-alpha_1)*exp(alpha_0*t_2);
b2 = -1/(alpha_0-alpha_1)*exp(alpha_1*t_2);

beta1 = -b2/(a2*b1-b2*a1);
beta2 = b1/(a2*b1-b2*a1);

gamma1 = -a2/(b2*a1-a2*b1);
gamma2 = a1/(b2*a1-a2*b1);

c_m_n(1,1) = beta1;
c_m_n(1,2) = beta2;
c_m_n(2,1) = gamma1;
c_m_n(2,2) = gamma2;

%test
failed_recon1 = 0;
for t = disc_d(1,1)*T_s:T_s:disc_d(1,2)*T_s
    test1 = beta1*(a1*exp(-alpha_0*t)+b1*exp(-alpha_1*t)) + beta2*(a2*exp(-alpha_0*t)+b2*exp(-alpha_1*t));
    test1_ideal = exp(-alpha_0*t);
    if round(real(test1),4)  ~= round(real(test1_ideal),4) || round(imag(test1),4)  ~= round(imag(test1_ideal),4)
        failed_recon1 = 1;
        failed_test1 = real(test1);
        failed_test1_ideal  = real(test1_ideal);
    end
end

failed_recon2 = 0;
for t = disc_d(1,1)*T_s:T_s:disc_d(1,2)*T_s
    test2 = gamma1*(a1*exp(-alpha_0*t)+b1*exp(-alpha_1*t)) + gamma2*(a2*exp(-alpha_0*t)+b2*exp(-alpha_1*t));
    test2_ideal = exp(-alpha_1*t);
    if round(real(test2),4)  ~= round(real(test2_ideal),4) || round(imag(test2),4)  ~= round(imag(test2_ideal),4)
        failed_recon2 = 1;
        failed_test2 = real(test2);
        failed_test2_ideal  = real(test2_ideal);
    end
end
end

