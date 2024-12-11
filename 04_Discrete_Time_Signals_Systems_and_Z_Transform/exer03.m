%% Exercise 3
clearvars;
close all;
clc;

student_id = 149372;
[original_signal, fs, noisy_signal] = degrade_fragment(student_id);

% Cross-correlate
correlation = xcorr(original_signal, noisy_signal);

% Find the index of the maximum correlation value
[~, max_corr_position] = max(correlation);

% Calculate the estimated starting sample position of the fragment
fragment_start_sample = max_corr_position - length(original_signal);

% Convert the sample position to time (in seconds)
fragment_start_time = fragment_start_sample / fs;

% Display the estimated starting time of the fragment
fprintf('Estimated fragment start time: %.5f seconds\n', fragment_start_time);