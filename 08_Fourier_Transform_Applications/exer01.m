clearvars;
close all;
clc;

load ecg.mat;

% Signal parameters
Fs = 360;  % Sampling frequency in Hz
N = length(sig);  % Signal length
f = (0:N-1)*(Fs/N);  % Frequency vector
f = f - Fs/2;  % Shift frequency range to [-Fs/2, Fs/2]
    
% Calculate FFT of the signal
sig_fft = fft(sig);
sig_fft_shifted = fftshift(sig_fft);

% 1. Remove frequencies below 0.5 Hz
cutoff_low = 0.5;
low_freq_indices = find(abs(f) < cutoff_low);
sig_fft_filtered = sig_fft_shifted;
sig_fft_filtered(low_freq_indices) = 0;

% 2. Identify and plot 60 Hz harmonics in spectrum
figure(1);
plot(f, abs(sig_fft_filtered));
title('Original Signal Spectrum with removed frequency');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
% Add markers for 60 Hz harmonics
hold on;
plot([60 60], ylim, 'r--');
plot([120 120], ylim, 'r--');
plot([180 180], ylim, 'r--');
plot([-60 -60], ylim, 'r--');
plot([-120 -120], ylim, 'r--');
plot([-180 -180], ylim, 'r--');
legend('Signal Spectrum', '60 Hz and Harmonics');
hold off;

% 3. Remove 60 Hz and harmonics
bandwidth = 0.4;
harmonics = [-180 -120 -60 60 120 180];

% Remove harmonics by directly zeroing out frequency components
for harm = harmonics
    % Find indices for positive and negative frequencies
    harm_indices_pos = find(abs(f - harm) < bandwidth);
    harm_indices_neg = find(abs(f + harm) < bandwidth);
    
    sig_fft_filtered(harm_indices_pos) = 0;
    sig_fft_filtered(harm_indices_neg) = 0;
end

% Convert back to time domain
sig_filtered = real(ifft(ifftshift(sig_fft_filtered)));

% Plot
figure(2);
subplot(2,1,1);
plot(tm, sig);
title('Original ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(tm, sig_filtered);
title('Filtered ECG Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;