%% Band-Pass Filter Design
clearvars;
close all;
clc;

% Bandpass Filter Specifications
fp1 = 12e6;     % First passband frequency (Hz)
fp2 = 25e6;     % Second passband frequency (Hz)
fs1 = 3e6;      % First stopband frequency (Hz)
fs2 = 30e6;     % Second stopband frequency (Hz)
Ap = 1;         % Passband ripple (dB)
As = 60;        % Stopband attenuation (dB)
fm = 80e6;      % Sampling frequency (Hz)

% rad/s
wp1 = 2*pi*fp1;
wp2 = 2*pi*fp2;
ws1 = 2*pi*fs1;
ws2 = 2*pi*fs2;

% Butterworth Filter Design
% Frequency pre-distortion
wap1_butter = 2*fm*tan(wp1/2/fm);
wap2_butter = 2*fm*tan(wp2/2/fm);
was1_butter = 2*fm*tan(ws1/2/fm);
was2_butter = 2*fm*tan(ws2/2/fm);

% Low-Pass Equivalent
wx2_butter = wap1_butter*wap2_butter;
if (was1_butter*was2_butter < wx2_butter)
    was1_butter = wx2_butter/was2_butter;
else
    was2_butter = wx2_butter/was1_butter;
end

wlpp_butter = wap2_butter - wap1_butter;
wlps_butter = was2_butter - was1_butter;

% Filter Order and Cutoff Calculation
e2_butter = 10^(0.1*Ap) - 1;
n_butter = ceil(log10(sqrt((10^(0.1*As)-1)/e2_butter)) / log10(wlps_butter/wlpp_butter));
wc_butter = wlpp_butter / (e2_butter^(1/(2*n_butter)));

[Z_butter, P_butter, K_butter] = buttap(n_butter);
[B_butter, A_butter] = zp2tf(Z_butter, P_butter, K_butter);
[B1_butter, A1_butter] = lp2bp(B_butter, A_butter, sqrt(wx2_butter), wc_butter);
[b_butter, a_butter] = bilinear(B1_butter, A1_butter, fm);

% Chebyshev Type I Filter Design
wap1_cheby1 = 2*fm*tan(wp1/2/fm);
wap2_cheby1 = 2*fm*tan(wp2/2/fm);
was1_cheby1 = 2*fm*tan(ws1/2/fm);
was2_cheby1 = 2*fm*tan(ws2/2/fm);

wx2_cheby1 = wap1_cheby1*wap2_cheby1;
if (was1_cheby1*was2_cheby1 < wx2_cheby1)
    was1_cheby1 = wx2_cheby1/was2_cheby1;
else
    was2_cheby1 = wx2_cheby1/was1_cheby1;
end

wlpp_cheby1 = wap2_cheby1 - wap1_cheby1;
wlps_cheby1 = was2_cheby1 - was1_cheby1;

eps_cheby1 = sqrt(10^(0.1*Ap) - 1);
n_cheby1 = ceil(acosh(sqrt((10^(0.1*As) - 1)/(eps_cheby1^2))) / acosh(wlps_cheby1/wlpp_cheby1));

[Z_cheby1, P_cheby1, K_cheby1] = cheb1ap(n_cheby1, eps_cheby1);
[B_cheby1, A_cheby1] = zp2tf(Z_cheby1, P_cheby1, K_cheby1);
[B1_cheby1, A1_cheby1] = lp2bp(B_cheby1, A_cheby1, sqrt(wx2_cheby1), wlpp_cheby1);
[b_cheby1, a_cheby1] = bilinear(B1_cheby1, A1_cheby1, fm);

% Chebyshev Type II Filter Design
wap1_cheby2 = 2*fm*tan(wp1/2/fm);
wap2_cheby2 = 2*fm*tan(wp2/2/fm);
was1_cheby2 = 2*fm*tan(ws1/2/fm);
was2_cheby2 = 2*fm*tan(ws2/2/fm);

wx2_cheby2 = wap1_cheby2*wap2_cheby2;
if (was1_cheby2*was2_cheby2 < wx2_cheby2)
    was1_cheby2 = wx2_cheby2/was2_cheby2;
else
    was2_cheby2 = wx2_cheby2/was1_cheby2;
end

wlpp_cheby2 = wap2_cheby2 - wap1_cheby2;
wlps_cheby2 = was2_cheby2 - was1_cheby2;

n_cheby2 = ceil(log10(10^(0.1*As) - 1) / (2*log10(wlps_cheby2/wlpp_cheby2)));

[Z_cheby2, P_cheby2, K_cheby2] = cheb2ap(n_cheby2, As);
[B_cheby2, A_cheby2] = zp2tf(Z_cheby2, P_cheby2, K_cheby2);
[B1_cheby2, A1_cheby2] = lp2bp(B_cheby2, A_cheby2, sqrt(wx2_cheby2), wlpp_cheby2);
[b_cheby2, a_cheby2] = bilinear(B1_cheby2, A1_cheby2, fm);

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
title('Magnitude Response of Band-Pass Filters');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Phase Response
figure;
plot(w, unwrap(angle(h_butter))*180/pi, 'b', 'LineWidth', 1.5);
hold on;
plot(w, unwrap(angle(h_cheby1))*180/pi, 'r', 'LineWidth', 1.5);
plot(w, unwrap(angle(h_cheby2))*180/pi, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Phase Response of Band-Pass Filters');
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
title('Group Delay of Band-Pass Filters');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Impulse Response
figure;
[h_imp_butter] = impz(b_butter, a_butter);
[h_imp_cheby1] = impz(b_cheby1, a_cheby1);
[h_imp_cheby2] = impz(b_cheby2, a_cheby2);

plot(h_imp_butter, 'b', 'LineWidth', 1.5);
hold on;
plot(h_imp_cheby1, 'r', 'LineWidth', 1.5);
plot(h_imp_cheby2, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Impulse Response of Band-Pass Filters');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

% Random Signal
N = 1000;
t = (0:N-1)/fm;
x = randn(1, N);

y_butter = filter(b_butter, a_butter, x);
y_cheby1 = filter(b_cheby1, a_cheby1, x);
y_cheby2 = filter(b_cheby2, a_cheby2, x);

% Plot
figure;
plot(t*1e6, y_butter, 'b', 'LineWidth', 1.5);
hold on;
plot(t*1e6, y_cheby1, 'r', 'LineWidth', 1.5);
plot(t*1e6, y_cheby2, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)');
ylabel('Amplitude');
title('Filtered Random Signals');
legend('Butterworth', 'Chebyshev I', 'Chebyshev II');

fprintf('Filter Orders:\n');
fprintf('Butterworth: %d\n', n_butter);
fprintf('Chebyshev I: %d\n', n_cheby1);
fprintf('Chebyshev II: %d\n', n_cheby2);

