ccc;
%%
%--------------------------------------------------------------------------
%----------------------------SET-UP----------------------------------------
%--------------------------------------------------------------------------

% Parameters
t_int = 15;                   % Temporal interval where Diracs can be located is [0 t_int]
K     = 4;                    % Total number of Diracs
K_B = 2;                      % Number of Diracs in a burst     
B = K/K_B;                    % Number of input bursts
P     = 2*K_B - 1;            % E-spline order
T_s   = 1 / 64;               % "continuous" time resolution
s = 2;                        % at most s-1 samples we cannot use in each burst
C_T = 0.1;                    % Threshold value of integrator
SNR = 50;                     % Signal-to-noise ratio

% E-spline of order P (can reproduce P+1 exponentials)
omega_0 = -pi/3;
lambda = -2*omega_0/P ;

%Channel 1
m = [0,3];
omega_m   = omega_0 + lambda * m;
alpha_vec1 = 1j * omega_m;
[phi1, t_phi1] = generate_e_spline(alpha_vec1, T_s, 1, 'anticausal');
L_phi1        = length(phi1);
h1            = real(phi1(end:-1:1));
t_h1          = -t_phi1(end:-1:1);

%Channel 2
m = [1,2];
omega_m   = omega_0 + lambda * m;
alpha_vec2 = 1j * omega_m;
[phi2, t_phi2] = generate_e_spline(alpha_vec2, T_s, 1, 'anticausal');
L_phi2        = length(phi2);
h2            = real(phi2(end:-1:1));
t_h2          = -t_phi2(end:-1:1);



%Minimum separation between consecutive bursts
min_separation = max(L_phi1, L_phi2);

%%
%--------------------------------------------------------------------------
%----------------------------INPUT ----------------------------------------
%--------------------------------------------------------------------------
% Stream of Diracs - input
t_sig = (0 : T_s : t_int)';
itk = [60 80 700 725];
t_k = t_sig(itk).';
a_k = [1 1 0.7 1 ]*1.1 ; 

% Generate the continuous-time signal x(t)
x = zeros(size(t_sig));
x(itk(1:K)) = a_k(1:K);
length_x = itk(end)+L_phi1;

% Add noise to the input
[x_noisy, noise, sigma_noise] = add_noise( x', SNR, length_x );
% x = x_noisy;

%%
%--------------------------------------------------------------------------
%----------------------------OUTPUT----------------------------------------
%--------------------------------------------------------------------------
%CHANNEL 1
% Compute the continuous time signal f(t) = x(t) * h(t)
%This is equivalent to the inner product between x and the E-spline phi
f1 = conv(x, h1);
t_start = t_sig(1) + t_h1(1);
t_end = t_sig(end) + t_h1(end);
t_f1 = t_start : T_s : t_end;
%Integrate the filtered signal
[ y1, t_y1 ] = integrate_directly( f1, t_f1, T_s );

%CHANNEL 2
% Compute the continuous time signal f(t) = x(t) * h(t)
%This is equivalent to the inner product between x and the E-spline phi
f2 = conv(x, h2);
t_start = t_sig(1) + t_h2(1);
t_end = t_sig(end) + t_h2(end);
t_f2 = t_start : T_s : t_end;
%Integrate the filtered signal
[ y2, t_y2 ] = integrate_directly( f2, t_f2, T_s );

%%
%--------------------------------------------------------------------------
%--------------------------NON-UNIFORM SAMPLES-----------------------------
%--------------------------------------------------------------------------
[ t_n1, y_n1, y_with_feedback1] = integrate_and_fire( y1, C_T );
[ t_n2, y_n2, y_with_feedback2] = integrate_and_fire( y2, C_T );

y_n1_temp = y_n1;
y_n2_temp = y_n2;

plot_all_signals( t_int,t_sig, x, t_y1, y_with_feedback1, f1, t_f1, t_y2, y_with_feedback2, f2, t_f2);    

%%
%--------------------------------------------------------------------------
%--------------------------SEQUENTIAL INPUT ESTIMATION---------------------
%--------------------------------------------------------------------------

%Initialize the estimated time locations t_est_k, and amplitudes a_est_k
t_est_k = [];
a_est_k = [];

for b = 1:B

    %Get the 2 samples in this burst, corresponding to the areas between output
    %spikes s and s+1, and output spikes s+1 and s+2
    [y_n1_burst, t_n1_burst] = get_samples_burst(b, t_est_k, T_s, min_separation, t_n1, y_n1, s, C_T, K_B);
    [y_n2_burst, t_n2_burst] = get_samples_burst(b, t_est_k, T_s, min_separation, t_n2, y_n2, s, C_T, K_B);

    %Get the coefficients that reproduce exponentials
    %CHANNEL 1
    t_n_0 = t_n1_burst(1,1)*T_s;
    t_n_1 = t_n1_burst(1,2)*T_s;
    t_n_2 = t_n1_burst(1,3)*T_s;
    [ c_m1 ] = get_c_m_n( alpha_vec1, t_n_0, t_n_1, t_n_2, T_s);
    %Save the first sample of this burst, of this channel
    t_n_0_channel1 = t_n_0;

    %CHANNEL 2
    t_n_0 = t_n2_burst(1,1)*T_s;
    t_n_1 = t_n2_burst(1,2)*T_s;
    t_n_2 = t_n2_burst(1,3)*T_s;
    [ c_m2 ] = get_c_m_n( alpha_vec2, t_n_0, t_n_1, t_n_2, T_s);
    %Save the first sample of this burst, of this channel
    t_n_0_channel2 = t_n_0;

    %Compute the signal moments
    %CHANNEL 1
    s_m1 = c_m1*y_n1_burst';
    s_m1_true = compute_true_moments( a_k(2*b-1:2*b), t_k(2*b-1:2*b), alpha_vec1);
    
    %CHANNEL 2
    s_m2 = c_m2*y_n2_burst';
    s_m2_true = compute_true_moments( a_k(2*b-1:2*b), t_k(2*b-1:2*b), alpha_vec2);
    
    %Overall moments
    s_m = [s_m1(1,1),s_m2(1,1),s_m2(2,1),s_m1(2,1)]';

    % Locate the Diracs using Prony's method
    S = toeplitz(s_m(K_B+1:end), s_m(K_B+1:-1:1));
    [~,~,V] = svd(S);
    h_prony       = V(:,end);
    uu_k    = roots(h_prony);
    angle_uu_k = mod(angle(uu_k),2*pi) ;
    
    %Estimate the Diracs' locations
    tt_k_burst = shift_estimated_time( t_est_k, uu_k, min_separation, T_s, lambda, t_n_0_channel1, t_n_0_channel2);
    
    %Estimate the Diracs' amplitudes
    aa_k_burst = real(retrieve_amplitudes(alpha_vec1, s_m1(1,1), s_m1(2,1), tt_k_burst));
    t_est_k = vertcat(t_est_k, tt_k_burst);
    a_est_k = vertcat(a_est_k, aa_k_burst);
    
    %Remove the contribution of the current estimation
    %CHANNEL 1
    y_n1 = remove_dirac_contribution( tt_k_burst, aa_k_burst, T_s, t_sig, y_n1, t_n1, h1, t_h1);
    
    %CHANNEL 2
    y_n2 = remove_dirac_contribution( tt_k_burst, aa_k_burst, T_s, t_sig, y_n2, t_n2, h2, t_h2);
end

%%
%--------------------------------------------------------------------------
%----------------------PLOT INPUT AND ESTIMATION---------------------------
%--------------------------------------------------------------------------
figure 
subplot(2,2,1)
stem(t_k,a_k, 'k', 'fill','Markersize',2)
axis([0 t_int -1.6 1.6]);
title('Input signal');
xlabel('(a)');

subplot(2,2,2)
plot(t_f1, f1, '-k');
hold on;
plot(t_y1, y_with_feedback1, 'b', 'Linewidth', 2);
hold off;
axis([0 t_int -2 2]);
Legend = cell(1,2);
Legend{1,2} = 'Output of integrator';
Legend{1,1} = 'Filtered input';
lgd = legend(Legend);
lgd.FontSize = 10;
title('Channel 1: filtered input and output of integrator');
xlabel('(b)');

subplot(2,2,3)
stem(t_n1*T_s, y_n1_temp, 'k','fill','Markersize',3); %note samples are not exactly equal to C_T due to finite time resolution determined by T_s
axis([0 t_int -C_T*(1+0.1) C_T*(1+0.1)]);
title('Channel 1: output non-uniform samples');
xlabel('(c)');

subplot(2,2,4)
stem(t_est_k, a_est_k, 'k', 'fill','Markersize',2)
axis([0 t_int -1.6 1.6]);
title('Reconstruction of the input signal');
xlabel('(d)');