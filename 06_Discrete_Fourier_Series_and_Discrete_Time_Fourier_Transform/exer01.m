clearvars;
close all;

%% Exercise 1 a)
load din.dat;
N = length(din);
fs = 4096;

ck = zeros(size(din));
x = din;

n = 0:N-1;
for i=1:N
    k = i - 1;
    ck(i) = 1/N * sum(x.' .* exp(-1i * 2 * pi * k * n / N));  % DFT calculation
end

% Apply FFT shift to the coefficients
ck = [ck(N/2+1:N); ck(1:N/2)];

% Frequency vector (centered at zero using fftshift)
k = 0:N-1;  % Frequency index
fk = k*fs/N - fs/2;  % Frequency vector (Hz)

% Plot the power spectrum
figure;
stem(fk, 20*log10(abs(ck))); 
grid on;
xlabel('Frequency (Hz)');
ylabel('Discrete Fourier Series Coeffs |c_k| (dB)');
title('Power Spectrum of the Discrete Signal');

%% Exercise 1 b)
delta_f = fs / N;
fprintf('Frequency Resolution: %.2f Hz\n', delta_f);

%% Exercise 1 c)
% Zero out frequencies beyond 500 Hz
threshold_freq = 500;  % Frequency threshold in Hz
valid_indices = abs(fk) <= threshold_freq;  % Logical array for frequencies within [-500, 500] Hz
ck(~valid_indices) = 0;  % Set coefficients outside 500 Hz to zero

% Compute the inverse DFT (IDFT) to get the time-domain signal back
ck_unshifted = [ck(N/2+1:N); ck(1:N/2)];  % Unshift the coefficients
x_reconstructed = real(N * ifft(ck_unshifted));  % Inverse DFT, scaled by N

% Plot the original and reconstructed signals
t = (0:N-1)/fs;  % Time vector

figure;
plot(t, din, 'b', 'LineWidth', 1.5);  % Original signal in blue
hold on;
plot(t, x_reconstructed, 'r--', 'LineWidth', 1.5);  % Reconstructed signal in red dashed line
hold off;

title('Original and Reconstructed Signal (after zeroing ck beyond 500 Hz)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Reconstructed Signal');
grid on;

%% Exercise 1 d)
% Compute the error signal
error_signal = din - x_reconstructed;

% Calculate the Signal-to-Error Ratio (SER)
ser = 10 * log10(sum(din.^2) / sum(error_signal.^2));
fprintf('Signal-to-Error Ratio (SER): %.2f dB\n', ser);

figure;
plot(t, error_signal, 'k', 'LineWidth', 1.5);  % Error signal in black
title('Error Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Exercise 1 e)
% The minimum required sampling rate is 2*500Hz = 1000 Hz
% Calculate the new data length based on the new sampling rate
% T = N / fs
% T = 2458 / 4096 ≈ 0.6 seconds

% Using the new sampling rate f_min = 1000 Hz, the new number of samples N_min would be:
%
% N_min = T × f_min
% N_min = 0.6 × 1000 = 600 samples