ccc;
%%
%--------------------------------------------------------------------------
%----------------------------SET-UP----------------------------------------
%--------------------------------------------------------------------------

% Parameters
t_int = 20;                   % Temporal interval where Diracs can be located is [0 t_int]
K     = 6;                    % Total number of Diracs
K_B = 2;                      % Number of Diracs in a burst     
B = K/K_B;                    % Number of input bursts
P     = 4;                    % E-spline order
T_s   = 1 / 64;               % "continuous" time resolution
s = 2;                        % at most s-1 samples we cannot use in each burst
C_T = 0.005;                  % Threshold value of integrator
noise_var =0.01;              % Variance of input additive white Gaussian noise

% E-spline of order P (can reproduce P+1 exponentials)
omega_0 = -pi/2;
lambda = -2*omega_0/(2*P-1);

%Channel 1
m_1 = 0:1:P/2-1;
m_1 = [m_1 2*P-P/2:1:2*P-1];
omega_m   = omega_0 + lambda * m_1;
alpha_vec1 = 1j * omega_m;
[phi1, t_phi1] = generate_e_spline(alpha_vec1, T_s, 1, 'anticausal');
L_phi1        = length(phi1);
h1            = real(phi1(end:-1:1));
t_h1          = -t_phi1(end:-1:1);
%Derivative of E-spline
[ dphi_dt_direct1, t_dphi_dt1] = derivative_e_spline( phi1, t_phi1, L_phi1, alpha_vec1, T_s );
h_derivative1 = real(dphi_dt_direct1(end:-1:1));
t_h_derivative1 = t_h1;


%Channel 2
m_2 = P/2:1:2*P-P/2-1;
omega_m   = omega_0 + lambda * m_2;
alpha_vec2 = 1j * omega_m;
[phi2, t_phi2] = generate_e_spline(alpha_vec2, T_s, 1, 'anticausal');
L_phi2        = length(phi2);
h2            = real(phi2(end:-1:1));
t_h2          = -t_phi2(end:-1:1);
%Derivative of E-spline
[ dphi_dt_direct2, t_dphi_dt2] = derivative_e_spline( phi2, t_phi2, L_phi2, alpha_vec2, T_s );
h_derivative2 = real(dphi_dt_direct2(end:-1:1));
t_h_derivative2 = t_h2;



%Minimum separation between consecutive bursts
min_separation = max(L_phi1, L_phi2);

%%
%--------------------------------------------------------------------------
%----------------------------INPUT ----------------------------------------
%--------------------------------------------------------------------------
% Stream of Diracs - input
t_sig = (0 : T_s : t_int)';
itk = [60  80 500 520 1000 1030];%EXP 1
t_k = t_sig(itk).';
a_k = [1 0.5 -1 -1 1.5 0.2] ; %EXP 1

% Generate the piecewise constant signal x(t)
x = create_piecewise_signal( t_sig, itk, a_k );
%Find the Diracs corresponding to the discontinuities in the piecewise
%signal x(t)
x_derivative = zeros(size(t_sig));
x_derivative(itk(1:K)) = a_k(1:K);


% Add noise to the input
[x_noisy, noise, SNR] = add_noise( x', noise_var);
% x = x_noisy;

%%
%--------------------------------------------------------------------------
%----------------------------OUTPUT----------------------------------------
%--------------------------------------------------------------------------
%CHANNEL 1
% Compute the continuous time signal f(t) = x(t) * h(t)
%This is equivalent to the inner product between x and the E-spline phi
f1 = conv(x, h_derivative1)*T_s;
t_start = t_sig(1) + t_h_derivative1(1);
t_end = t_sig(end) + t_h_derivative1(end);
t_f1 = t_start : T_s : t_end;
%Integrate the filtered signal
[ y1, t_y1 ] = integrate_directly( f1, t_f1, T_s );

%CHANNEL 2
% Compute the continuous time signal f(t) = x(t) * h(t)
%This is equivalent to the inner product between x and the E-spline phi
f2 = conv(x, h_derivative2)*T_s;
t_start = t_sig(1) + t_h_derivative2(1);
t_end = t_sig(end) + t_h_derivative2(end);
t_f2 = t_start : T_s : t_end;
%Integrate the filtered signal
[ y2, t_y2 ] = integrate_directly( f2, t_f2, T_s );

%%
%--------------------------------------------------------------------------
%--------------------------NON-UNIFORM SAMPLES-----------------------------
%--------------------------------------------------------------------------
[ t_n1, y_n1, y_with_feedback1 ] = integrate_and_fire( y1, C_T );
[ t_n2, y_n2, y_with_feedback2 ] = integrate_and_fire( y2, C_T );

y_n1_temp = y_n1;
y_n2_temp = y_n2;

plot_all_signals( t_sig, x, t_y1, y_with_feedback1, f1, t_f1, t_int); 
plot_all_signals( t_sig, x, t_y2, y_with_feedback2, f2, t_f2, t_int); 

%%
%--------------------------------------------------------------------------
%--------------------------SEQUENTIAL INPUT ESTIMATION---------------------
%--------------------------------------------------------------------------

%Initialize the estimated time locations t_est_k, and amplitudes a_est_k
tt_k = [];
aa_k = [];
nonuniform_times_tot1 = [];
nonuniform_times_tot2 = [];
y_n_tot1 = [];
y_n_tot2 = [];
for c = 1:B %FOR EACH BURST
    
    %Get samples for channel 1
    [ nonuniform_times_c1, y_n_c1 ] = get_samples_burst( c, T_s, L_phi1, y_n1, t_n1, P, tt_k, K_B);
    nonuniform_times_tot1 = [nonuniform_times_tot1, nonuniform_times_c1(2:end)];
    y_n_tot1 = [y_n_tot1, y_n_c1];
    
    %Get samples for channel 1
    [ nonuniform_times_c2, y_n_c2 ] = get_samples_burst( c, T_s, L_phi2, y_n2, t_n2, P, tt_k, K_B);
    nonuniform_times_tot2 = [nonuniform_times_tot2, nonuniform_times_c2(2:end)];
    y_n_tot2 = [y_n_tot2, y_n_c2];

    %Least-squares coefficients
    [c_ls_m1, error_approximation_dirac1] = get_least_squares_c_m_n( alpha_vec1, h1, t_h1, nonuniform_times_c1*T_s, T_s, L_phi1, P);
    [c_ls_m2, error_approximation_dirac2] = get_least_squares_c_m_n( alpha_vec2, h2, t_h2, nonuniform_times_c2*T_s, T_s, L_phi2, P);

    %Compute the signal moments
    s_m1 = c_ls_m1*y_n_c1';
    s_m2 = c_ls_m2*y_n_c2';
    
    %s_m_true = compute_true_moments( a_k, t_k, alpha_vec);
    s_m1_1  = s_m1(1,1);
    s_m2_1  = s_m2(1,1);

    % Compute signal moments
    s_m(1:P/2,1) = s_m1(1:P/2)';
    s_m(2*P-P/2+1:1:2*P,1) = s_m1(P/2+1:1:P)';
    s_m(P/2+1:1:2*P-P/2,1) = s_m2'; 
    
    
    S = toeplitz(s_m(3:end), s_m(3:-1:1));
    [~,~,V] = svd(S);
    h_prony       = V(:,end);
    uu_k    = roots(h_prony);
    angle_uu_k = mod(-angle(uu_k),2*pi) ;
    %Dirac's location
    if c >1
        t_start =  tt_k(K_B*(c-1))+L_phi1*T_s;
    else
        t_start = 0;
    end
     tt_k_burst= shift_estimated_time( tt_k, uu_k, min_separation, T_s, lambda, nonuniform_times_c1(1,1)*T_s, nonuniform_times_c2(1,1)*T_s);

    %Dirac's amplitude
    aa_k_burst = -real(retrieve_amplitudes(alpha_vec1,s_m1(1,1),s_m1(2,1),K_B,tt_k_burst));
    tt_k = vertcat(tt_k, tt_k_burst);
    aa_k = vertcat(aa_k, aa_k_burst); 
end


%Reconstructed piecewise constant signal
x_est = create_piecewise_signal( t_sig, round(tt_k'/T_s), aa_k' );




%%
%--------------------------------------------------------------------------
%----------------------PLOT INPUT AND ESTIMATION---------------------------
%--------------------------------------------------------------------------
figure 
subplot(2,2,1)
plot(t_sig, x, 'k')
xlabel('(a)');
 axis([0 t_int -1 1.7]);

subplot(2,2,2)
plot(t_dphi_dt1, dphi_dt_direct1, 'k');
% axis([-4 0 -0.6 0.6]);
xlabel('(b)');

subplot(2,2,3)
stem(t_n1*T_s, y_n1./abs(y_n1)*C_T, 'k.');
%  stem(nonuniform_times_tot*T_s, y_n_tot./abs(y_n_tot)*C_T, 'k.');
xlabel('(c)');
% axis([0 t_int -C_T C_T]);

subplot(2,2,4)
plot(t_sig, x_est, 'k')
xlabel('(d)');
axis([0 t_int -1 1.7]);


figure 
subplot(1,3,1)
plot(t_sig, x, 'Linewidth', 2);
hold on;
plot(t_f1(1:end-1), f1);
hold off;
axis([0 t_int -1 1.7]);
xlabel('(a)')
Legend = cell(1,2);
Legend{1,1} = 'Input signal';
Legend{1,2} = 'Filtered input';
lgd = legend(Legend);

subplot(1,3,2)
% stem(t_n1*T_s, y_n1./abs(y_n1)*C_T, 'k.');
stem(nonuniform_times_tot1*T_s, y_n_tot1./abs(y_n_tot1)*C_T, 'k.')
xlabel('(b)');
axis([1 2 -C_T C_T]);

subplot(1,3,3)
plot(t_sig, x_est, 'k')
xlabel('(c)');
axis([0 t_int -1 1.7]);

