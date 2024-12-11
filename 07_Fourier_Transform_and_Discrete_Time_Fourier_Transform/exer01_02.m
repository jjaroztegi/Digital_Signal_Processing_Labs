clearvars;
close all;
clc;

%% exer01
% Define syms and parameter
syms t f t0
assume(t0, 'positive');
t0 = 1e-3; 

% Part 1a: Define the signal x(t) and compute its Fourier Transform X(f)
x_t = exp(-t^2 / t0^2);
X_f = int(x_t * exp(-1i * 2 * pi * f * t), t, -inf, inf); % FT
X_f = simplify(X_f);
% X_f = (pi^(1/2)*exp(-(f^2*pi^2)/1000000))/1000
% X_f = t0*sqrt(sym(pi))*exp((-f^2*t0^2*sym(pi)^2))

disp('Fourier Transform X(f):');
disp(X_f);

% Part 1b: Define anonymous functions for x(t) and X(f)
x_anom = @(t0, t) exp(-t.^2 / t0^2);
X_anom = @(t0, f) (t0*sqrt(sym(pi))*exp((-f^2*t0^2*sym(pi)^2)));


% Part 1c: Compute the energy of the signal x(t) analytically
energy_x = int(x_t^2, t, -inf, inf);
energy_x = simplify(energy_x);
energy_x = eval(energy_x);

% Display the energy
disp('Energy of the signal x(t):');
disp(energy_x);

% Part 1d: Determine fmax such that energy for f > fmax is 120 dB below total energy
total_energy = double(energy_x);
threshold_energy = total_energy * 10^(-120/10); % Convert dB to linear scale

f0 = 0;                         % Starting frequency
df = 50;                        % Step size for frequency increment
energy_above_f0 = total_energy; % Start with total energy

while energy_above_f0 > threshold_energy
    energy_above_f0 = double(int(abs(X_anom(t0, f)).^2, f, f0, inf));
    f0 = f0 + df;
end

fmax = f0 - df;

disp('fmax (120 dB below total energy):');
disp(fmax);

% Part 1e: Select fs
fs = 2 * fmax;
disp('Sampling frequency fs:');
disp(fs);

% Part 1f: SNR = 60 dB
threshold_energy_60dB = total_energy * 10^(-60/10);
f0 = 0;
energy_above_f0 = total_energy;

while energy_above_f0 > threshold_energy_60dB
    energy_above_f0 = double(int(abs(X_anom(t0, f)).^2, f, f0, inf));
    f0 = f0 + df;
end

fmax_60dB = f0 - df;

disp('fmax (60 dB below total energy):');
disp(fmax_60dB);

% Select appropriate sampling frequency fs for SNR = 60 dB
fs_60dB = 2 * fmax_60dB;
disp('Sampling frequency fs (60 dB below total energy):');
disp(fs_60dB);

%% exer02

% interval [-3t0, 3t0]
n_samples = ceil(6 * t0 * fs_60dB); % Total number of samples, symmetric around 0
t_sample = linspace(-3*t0, 3*t0, n_samples);
x_n = x_anom(t0, t_sample);

disp('Number of sampling points (n_samples):');
disp(n_samples);

% DTFT of x[n]
% F = [-1/2, 1/2]
F = linspace(-0.5, 0.5, 1000);

X_F = zeros(size(F));
for k = 1:length(F)
    X_F(k) = sum(x_n .* exp(-1i * 2 * pi * F(k) * (0:n_samples-1)));
end

% Normalize
X_F = X_F / max(abs(X_F));


% Part1 X(f) over the same frequency range for comparison
% SNR 60dB
f_range60 = F * fs_60dB;
X_anom = @(t0, f) (t0*sqrt(sym(pi))*exp((-f.^2*t0.^2*sym(pi).^2)));

X_f_values = double(X_anom(t0, f_range60));
% Normalize
X_f_values = X_f_values / max(abs(X_f_values));


% SNR 120dB
f_range120 = F * fs;

X_f_values120 = double(X_anom(t0, f_range120));
% Normalize
X_f_values120 = X_f_values120 / max(abs(X_f_values120));

figure;
plot(F, abs(X_F), 'LineWidth', 2);
hold on;
plot(F, abs(X_f_values), 'r--', 'LineWidth', 2);
hold on;
plot(F, abs(X_f_values120), 'g--', 'LineWidth', 2);

xlabel('Frequency (Hz)');
ylabel('|X(f)| and |X(F)|');
title('Comparison of DTFT and FT');
legend('DTFT |X(F)|', 'FT |X(f)| 60dB', 'FT |X(f)| 120dB');
grid on;

% Part 2d:
% X(f) is the Continuous Fourier Transform of x(t), representing all
% frequency components In contrast, X(F) is the DTFT of the sampled version
% of x(t) and is limited to the frequency range [-fs/2, fs/2]. Sampling
% restricts X(F) to the Nyquist range, which causes it to differ from X(f)
% due to the finite frequency range.

% The FT curve for 120 dB SNR appears narrower in frequency, while the 60
% dB FT curve is wider due to the larger fmax. This shows that stricter SNR
% thresholds require sampling across a broader frequency range to maintain
% high fidelity in representing the original continuous signal.