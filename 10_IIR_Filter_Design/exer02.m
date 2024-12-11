%% Low-Pass Filter Design
clearvars;
close all;
clc;

% Load the signal
load('din.mat');

% Filter Specifications
fp = 1.8e6;     % Passband frequency (Hz)
fs = 6e6;       % Stopband frequency (Hz)
Ap = 1;         % Passband ripple (dB)
As = 60;        % Stopband attenuation (dB)
fm = 32e6;      % Sampling frequency (Hz)

wp = 2*pi*fp;
ws = 2*pi*fs;

% Butterworth Filter Design
% Frequency pre-distortion
wap_butter = 2*fm*tan(wp/2/fm);
was_butter = 2*fm*tan(ws/2/fm);

% Filter Order and Cutoff Calculation
e2_butter = 10^(0.1*Ap) - 1;
n_butter = ceil(log10(sqrt((10^(0.1*As)-1)/e2_butter)) / log10(was_butter/wap_butter));
wc_butter = wap_butter / (e2_butter^(1/(2*n_butter)));

[Z_butter, P_butter, K_butter] = buttap(n_butter);
[B_butter, A_butter] = zp2tf(Z_butter, P_butter, K_butter);
[b_butter, a_butter] = bilinear(B_butter, A_butter, fm);

% Chebyshev Type I Filter Design
eps_cheby1 = sqrt(10^(0.1*Ap) - 1);
n_cheby1 = ceil(acosh(sqrt((10^(0.1*As) - 1)/(eps_cheby1^2))) / acosh(was_butter/wap_butter));

[Z_cheby1, P_cheby1, K_cheby1] = cheb1ap(n_cheby1, eps_cheby1);
[B_cheby1, A_cheby1] = zp2tf(Z_cheby1, P_cheby1, K_cheby1);
[b_cheby1, a_cheby1] = bilinear(B_cheby1, A_cheby1, fm);

% Chebyshev Type II Filter Design
n_cheby2 = ceil(acosh(sqrt((10^(0.1*As) - 1) / eps_cheby1^2)) / acosh(was_butter/wap_butter));

[Z_cheby2, P_cheby2, K_cheby2] = cheb2ap(n_cheby2, As);
[B_cheby2, A_cheby2] = zp2tf(Z_cheby2, P_cheby2, K_cheby2);
[b_cheby2, a_cheby2] = bilinear(B_cheby2, A_cheby2, fm);

% Frequency Response Analysis
nfft = 1024;
[h_butter, w] = freqz(b_butter, a_butter, nfft, fm);
[h_cheby1, ~] = freqz(b_cheby1, a_cheby1, nfft, fm);
[h_cheby2, ~] = freqz(b_cheby2, a_cheby2, nfft, fm);

% Magnitude Response
figure;
plot(w, 20*log10(abs(h_butter)), 'b', 'LineWidth', 1.5);
hold on;
plot(w, 20*log10(abs(h_cheby1)), 'r', 'LineWidth', 1.5);
plot(w, 20*log10(abs(h_cheby2)), 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response of Low-Pass Filters');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Group Delay
[gd_butter, w_gd] = grpdelay(b_butter, a_butter, nfft, fm);
[gd_cheby1, ~] = grpdelay(b_cheby1, a_cheby1, nfft, fm);
[gd_cheby2, ~] = grpdelay(b_cheby2, a_cheby2, nfft, fm);

figure;
plot(w_gd, gd_butter*1e6, 'b', 'LineWidth', 1.5);
hold on;
plot(w_gd, gd_cheby1*1e6, 'r', 'LineWidth', 1.5);
plot(w_gd, gd_cheby2*1e6, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Group Delay (µs)');
title('Group Delay of Low-Pass Filters');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Filter the signal
y_butter = filter(b_butter, a_butter, din);
y_cheby1 = filter(b_cheby1, a_cheby1, din);
y_cheby2 = filter(b_cheby2, a_cheby2, din);

% Time axis
t = (0:length(din)-1)/fm;

% Plot filtered signals
figure;
plot(t*1e6, y_butter, 'b', 'LineWidth', 1.5);
hold on;
plot(t*1e6, y_cheby1, 'r', 'LineWidth', 1.5);
plot(t*1e6, y_cheby2, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)');
ylabel('Amplitude');
title('Filtered GMSK Signal');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Print filter orders
fprintf('Filter Orders:\n');
fprintf('Butterworth: %d\n', n_butter);
fprintf('Chebyshev I: %d\n', n_cheby1);
fprintf('Chebyshev II: %d\n', n_cheby2);