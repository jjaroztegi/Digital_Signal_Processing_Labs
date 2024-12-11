clearvars;
close all;
clc;


% Parameters
Nw = 1024;          % Width of each portion (samples)
Sw = Nw/8;          % Overlap size
alpha = 0.8;        % Filter coefficient for smoothing spectra

[y, Fs] = audioread('audio_file.wav');

y1 = y(:,1);         % Left
y2 = y(:,2);         % Right

L = floor((length(y) - Sw)/(Nw - Sw));  % Number of portions

yw1 = zeros(Nw, L);
yw2 = zeros(Nw, L);
Yw1 = zeros(Nw, L);
Yw2 = zeros(Nw, L);
Yw_filtered1 = zeros(Nw/2 + 1, L);
Yw_filtered2 = zeros(Nw/2 + 1, L);

% Hanning window
window = hanning(Nw);

% 1. Divide signal and apply window
for l = 0:L-1
    start_idx = l*(Nw-Sw) + 1;
    end_idx = start_idx + Nw - 1;
    
    if end_idx <= length(y)
        yw1(:,l+1) = y1(start_idx:end_idx) .* window;
        yw2(:,l+1) = y2(start_idx:end_idx) .* window;
    end
end

% 2. Calculate FFT
for l = 1:L
    Yw1(:,l) = fft(yw1(:,l));
    Yw2(:,l) = fft(yw2(:,l));
    
    Yw_filtered1(:,l) = abs(Yw1(1:Nw/2+1,l));
    Yw_filtered2(:,l) = abs(Yw2(1:Nw/2+1,l));
end

% 3. Filter
filtered_spectrum1 = Yw_filtered1(:,1);
filtered_spectrum2 = Yw_filtered2(:,1);

for l = 2:L
    filtered_spectrum1 = alpha * filtered_spectrum1 + (1-alpha) * Yw_filtered1(:,l);
    filtered_spectrum2 = alpha * filtered_spectrum2 + (1-alpha) * Yw_filtered2(:,l);

    Yw_filtered1(:,l) = filtered_spectrum1;
    Yw_filtered2(:,l) = filtered_spectrum2;
end

% 4. Create animation
figure('Name', 'Spectral Analysis', 'NumberTitle', 'off');

freq_axis = linspace(0, Fs/2, Nw/2 + 1);

for l = 1:L
    plot(freq_axis, 20*log10(Yw_filtered1(:,l)), 'b-', ...
         freq_axis, 20*log10(Yw_filtered2(:,l)), 'r-');
    
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(sprintf('Spectrum of Portion %d of %d (Left and Right channels)', l, L));
    legend('Left Channel', 'Right Channel');
    
    ylim([-80 40]);
    grid on;
    
    pause(0.05);

    set(gca, 'XTick', 0:1000:Fs/2);
    xlabel('Frequency (Hz)');

end
