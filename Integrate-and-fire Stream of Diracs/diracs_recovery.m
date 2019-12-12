ccc;
%%
%--------------------------------------------------------------------------
%----------------------------SET-UP----------------------------------------
%--------------------------------------------------------------------------

% Parameters
t_int = 10;             % Temporal interval where Diracs can be located
K     = 2;              % Number of Diracs
P     = 1;              % E-spline order
T_s   = 1 / 64;         % "continuous" time resolution
SNR   = 20;             % Signal-to-noise ratio

% E-spline of order P (can reproduce P exponentials)
m = 0:P;
omega_0 = -pi/3;
lambda = -2*omega_0;
omega_m   = omega_0 + lambda * m;
alpha_vec = 1j * omega_m;
[phi, t_phi] = generate_e_spline(alpha_vec, T_s, 1, 'anticausal');
L_phi        = length(phi);
h            = real(phi(end:-1:1));
t_h          = -t_phi(end:-1:1);

% Threshold value of integrate-and-fire TEM
C_T = 0.05;

%%
%--------------------------------------------------------------------------
%----------------------------INPUT ----------------------------------------
%--------------------------------------------------------------------------
% Stream of Diracs - input
t_sig = (0 : T_s : t_int)';
itk = [95 250];  %time locations of the Diracs
t_k = t_sig(itk).'; 
a_k = [1 1]; %amplitudes of the Diracs

% Generate the continuous-time signal x(t)
x = zeros(size(t_sig));
x(itk(1:K)) = a_k(1:K);
length_x = itk(end)+L_phi;

% Add noise to the input signal
[x_noisy, noise, sigma_noise] = add_noise( x', SNR, length_x );
% x = x_noisy;

%%
%--------------------------------------------------------------------------
%----------------------------OUTPUT----------------------------------------
%--------------------------------------------------------------------------
% Compute the continuous time signal f(t) = x(t) * h(t)
%This is equivalent to the inner product between x(t) and the E-splinephi(t)
f = conv(x, h);
t_start = t_sig(1) + t_h(1);
t_end = t_sig(end) + t_h(end);
t_f = t_start : T_s : t_end;

%Integrate the filtered signal
[ y, t_y ] = integrate_directly( f, t_f, T_s );

%%
%--------------------------------------------------------------------------
%--------------------------NON-UNIFORM SAMPLES-----------------------------
%--------------------------------------------------------------------------
%Reset the integrated signal, after each spike
[ t_n, y_n, y_with_reset] = integrate_and_fire( y, C_T);

plot_all_signals( t_int, t_y, y_with_reset, f, t_f);    

%%
%--------------------------------------------------------------------------
%--------------------------SEQUENTIAL INPUT ESTIMATION---------------------
%--------------------------------------------------------------------------

tt_k = [];
aa_k = [];

for b = 1:K %Sequentially retrieve the Diracs
    
    [y_n_current, t_n_current] = get_samples_current_dirac(b, tt_k, T_s, L_phi, t_n, y_n);    

    %Get the coefficients that reproduce exponentials
    t_n_0 = t_n_current(1,1)*T_s;
    t_n_1 = t_n_current(1,2)*T_s;
    t_n_2 = t_n_current(1,3)*T_s;
   
%     [ c_m ] = get_c_m_n( alpha_vec, t_n_0, t_n_1, t_n_2); %Compute the coefficients for a spline of order P=1
    [ c_m ] = get_least_squares_c_m_n( alpha_vec, h, t_h, t_n_current*T_s, T_s, L_phi); %Compute least square coefficients, by solving a system of equations

    %Compute the signal moments
    s_m1 = c_m*y_n_current';
    s_m = vertcat(s_m1(2,1), s_m1(1,1));
    
    % Locate the Diracs using Prony's method
    S = toeplitz(s_m(2:end), s_m(2:-1:1));
    [~,~,V] = svd(S);
    h_prony = V(:,end);
    uu_k    = roots(h_prony);
    
    %Dirac's location
    tt_k_current = shift_estimated_time( tt_k, uu_k, L_phi, T_s, lambda, t_n_0);
    
    %Dirac's amplitude
    aa_k_current = real(retrieve_amplitudes(alpha_vec, s_m(1,1), tt_k_current));
    tt_k = vertcat(tt_k, tt_k_current);
    aa_k = vertcat(aa_k, aa_k_current);

    plot_shifted_kernels_and_diracs( t_k, a_k, t_n, y_n, t_n_current, T_s, C_T, phi, t_phi )

end

%%
%--------------------------------------------------------------------------
%----------------------PLOT INPUT AND ESTIMATION---------------------------
%--------------------------------------------------------------------------
figure 
subplot(2,2,1)
stem(t_k,a_k,'k','fill','Markersize',3)
axis([0 t_int+L_phi*T_s -2 2]);
title('Input signal');
xlabel('(a)');

subplot(2,2,2)
plot(t_f, f, '-k');
hold on;
plot(t_y, y_with_reset, 'b', 'Linewidth', 2);
hold off;
axis([0 t_int+L_phi*T_s -1.2 1.7]);
Legend = cell(1,2);
Legend{1,2} = 'Output of integrator';
Legend{1,1} = 'Filtered input';
lgd = legend(Legend);
lgd.FontSize = 10;
title('Filtered input and output of integrator')
xlabel('(b)');

subplot(2,2,3)
stem(t_n*T_s, y_n./abs(y_n)*C_T, 'k','fill','Markersize',3);
axis([0 t_int+L_phi*T_s -C_T-0.1 C_T+0.1]);
title('Output spikes')
xlabel('(c)');

subplot(2,2,4)
stem(tt_k,aa_k,'k','fill','Markersize',3)
axis([0 t_int+L_phi*T_s -2 2]);
title('Reconstructed signal');
xlabel('(d)');

