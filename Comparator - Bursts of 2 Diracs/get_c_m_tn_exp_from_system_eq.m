function [c_m_n] = get_c_m_tn_exp_from_system_eq(P, alpha_m, disc_d, n_vec_d, no_kernels, phi, t_phi)
%N is the number of shifted kernels overlapping the interval between the
%two consecutive discontinuities in disc_d

% Rearrange the arguments (n row vector, alpha_m column vector)
disc_d = disc_d(:).';
T_s     = t_phi(2) - t_phi(1);

% Output matrix
c_m_n = zeros(P+1, no_kernels);

% Kernel's boundaries
lower_bound = t_phi(1) ;
upper_bound = t_phi(end) ;

%--------------------------------------------------------------------------
%------------------CREATE SYSTEM OF EQUATIONS------------------------------
%--------------------------------------------------------------------------
%Ax=B, find x

%Range of times t 
%t should take N values to solve the system of equations
%of N unknown coefficients
%t should be larger than the peak time of the first kernel
n_first = disc_d(1,1) ;
n_last = disc_d(1,2);
n_inc = floor((n_last - n_first)/(no_kernels * T_s*1));
%Times where we want to achieve perfect reconstruction
t_range = n_first + n_inc*T_s:n_inc*T_s:n_first+n_inc*T_s*no_kernels;
% t_range = [n_first+n_inc*T_s, n_last-n_inc*T_s];
for m = 0:P
    A = zeros(no_kernels,no_kernels);
    B = zeros(no_kernels,1);
    count = 1;
    for t = t_range
        %First, get the values of the kernel to form matrix A
        t_n = (t-n_vec_d)*T_s;
        idx_t_n = find(t_n >= lower_bound & t_n <= upper_bound);
        l_t = t_n(idx_t_n);

        % Find indices in phi which correspond to the times tn within the support
        % of the kernel phi. 
        %i.e. we need to find the times tn which are within the range [-1,1] and
        %the value of the kernel at those points
        idx   = round( (-lower_bound + l_t) / T_s) + 1;
        phi_l = phi(idx);
        
        A(count,idx_t_n) = phi_l;
        
        %Then, form matrix B
        B(count,1) = exp(-alpha_m(1,m+1)*t*T_s);
        
        count = count+1;

    end
    
    [~,colind] = rref(A);
    C = A(:, colind); 
    
    c_m_n(m+1,:) = A\B;
  
end


