clearvars;
close all;

%% Part 1 a)

syms t T t0 real
syms k
syms ck(k,to,T)
assume(k, 'integer');
assume(t0, 'positive');
assume(T, 'positive');

ck(k,to,T) = 1/T *  int((cos(2*pi*t/t0))^2 * exp(-1i*2*pi*k*t/T), t, -t0/4, t0/4);
ck(k,to,T) = simplify(ck(k,to,T));
% ck = (2*T^2*sin((sym(pi)*k*t0)/(2*T)))/(k*sym(pi)*(4*T^2 - k^2*t0^2))

% see attached pdf for steps

% Final Expression:
% c_k(k, t0, T) = (t0 / (4 * T)) * sinc(k * t0 / 2T) - (t0 / (4 * T)) * (k^2 * sinc(k * t0 / 2T)) / (k^2 - (4 * T^2 / t0^2))

%% Part 1 b)
x_t = @(t, t0) (abs(t) <= t0/4) .* (cos(2*pi*t/t0).^2);
c_k = @(k, t0, T) (t0/(4*T)) * sinc(k*t0/(2*T)) - (t0/(4*T)) * (k.^2 .* sinc(k*t0/(2*T))) ./ (k.^2 - (4*T^2/t0^2));

%% Part1 c)

t0 = 2.7e-3;
T = 6e-3;

fo = 1/T;   % Fundamental Frequency = Frequency resolution
N = 50;     % Starting value for N to iterate
SNR = 0;    % Starting value of SNR to enter the while loop

while (SNR < 96)    
    N = N + 1;                  % Increase N by 1
    k = -N:N;                   % k vector from -N to N
    ck = c_k(k, t0, T);        % Evaluate ck coefficients

    fmax = N * fo;             % The max frequency is N times the frequency resolution
    dt = 1/(10 * fmax);        % Oversampling by 10 => Sampling frequency is 10*2*fmax
    t = -T/2:dt:T/2;           % Time vector with time resolution dt
    
    % Compute the reconstruction using just 2N+1 coefficients
    xr = zeros(size(t)); 
    for i = 1:length(k)
        xr = xr + ck(i) * exp(1i * 2 * pi * k(i) * t / T);
    end
    xr = real(xr);              % Ensures that the signal is real
    
    xt = x_t(t, t0);            % Exact signal x(t)
    Pxt = mean(xt.^2);          % Power of x(t)
    
    e = xt - xr;                % Error Signal
    Pe = mean(e.^2);            % Power of error signal
    
    SNR = 10 * log10(Pxt / Pe); % SNR
end

figure; 
plot(t, xt, t, xr, '.'); 
grid on;
xlabel('Time (s)'); 
ylabel('x(t) signal'); 
legend('Original', 'Reconstruction');

disp(['There are ', num2str(2 * N + 1), ' coefficients']);
disp(['Maximum frequency is ', num2str(fmax), 'Hz'])

%% Part1 d)
% fs > 2*N*fo

M = 2 * N + 1;
fs = M / T;
dt = 1 / fs;               % Time step
t = linspace(-T/2, T/2 - dt, M);  % -T/2 to T/2
xt = x_t(t,t0);

k = -2 * N : 2 * N;        % Frequency index vector
ck_numerical = zeros(1, length(k));

% Compute the numerical Fourier coefficients
for i = 1:length(k)
    integrand = xt .* exp(-1i * 2 * pi * k(i) * t / T);  % Integrand for each k
    ck_numerical(i) = sum(integrand) / M;                % Numerical Fourier coefficient
end

% Plot analytical vs numerical Fourier coefficients
figure;
plot(k * fo, 20 * log10(abs(c_k(k, t0, T))), 'b');  % Analytical Fourier coefficients
hold on;
plot(k * fo, 20 * log10(abs(ck_numerical)), 'r');   % Numerical Fourier coefficients
grid on;
xlabel('Frequency (Hz)');
ylabel('Fourier Coefficients |c_k| (dB)');
legend('Analytical', 'Numerical');


%% Part 2 b)
X_f = @(f,t0) (1/t0)*sinc(f*t0/2)./((2/t0)^2-f.^2);

%% Part2 c)

E_x = sum(xt.^2) * dt;
disp(['Energy of x(t) is ', num2str(E_x)]);

%% Part2 d)
% Ex =  Intgral of  |X(f)|^2
df = 50;    % Frequency resolution step
f_max = 20400;   % Initial value for maximum frequency
SNR = 0;         % Initialize SNR to enter the loop
f_max = 0;       % Start at zero

while SNR < 96
    f_max = f_max + df;  % Increment fmax
    f = -f_max : df : f_max; % Frequency vector
    X_f_vals = X_f(f, t0);
    E_f = trapz(abs(X_f_vals).^2) * df;  % Compute energy in frequency domain
    SNR = 10 * log10(E_x / (E_x - E_f));
end

disp(['Bandwidth of the signal is ', num2str(2 * f_max), 'Hz']);


%% e)
f_max = ceil(f_max/df) * df;
f = -fmax : df : fmax;
x_approx = zeros(1, length(t));  %Reconstruction signal

% Perform the Inverse Fourier Transform using the computed bandwidth
for i = 1:length(t)
    integrand = X_f(f, t0) .* exp(1i * 2 * pi * f * t(i));
    x_approx(i) = trapz(integrand) * df;
end
x_approx = real(x_approx);  % Discard imaginary part due to numerical errors

figure;
plot(t, xt, 'b', t, x_approx, 'r--');
grid on;
xlabel('Time (s)');
ylabel('x(t) and reconstructed x(t)');
legend('Original x(t)', 'Reconstructed x(t)');
title('Comparison of Original and Reconstructed Signal');