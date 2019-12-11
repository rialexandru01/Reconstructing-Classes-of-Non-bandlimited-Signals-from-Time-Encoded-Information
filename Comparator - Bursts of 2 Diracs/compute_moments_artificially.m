function [ s_m ] = compute_moments_artificially( a_k, t_k, alpha_vec1, alpha_vec2)
%COMPUTE_MOMENTS_ARTIFICIALLY Summary of this function goes here
%   Detailed explanation goes here

s_m = zeros(4,1);
a_1 = a_k(1,1);
a_2 = a_k(1,2);
t_1 = t_k(1,1);
t_2 = t_k(1,2);
s_m(1,1) = a_1*exp(alpha_vec1(1,1)*t_1) + a_2*exp(alpha_vec1(1,1)*t_2);
s_m(2,1) = a_1*exp(alpha_vec2(1,1)*t_1) + a_2*exp(alpha_vec2(1,1)*t_2);
s_m(3,1) = a_1*exp(alpha_vec2(1,2)*t_1) + a_2*exp(alpha_vec2(1,2)*t_2);
s_m(4,1) = a_1*exp(alpha_vec1(1,2)*t_1) + a_2*exp(alpha_vec1(1,2)*t_2);


end

