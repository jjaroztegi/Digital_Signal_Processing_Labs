%% Exercise 01
% Parameters
n = 0:500; 
% x1[n] = cos(0.11*pi*n)
x1 = @(n) cos(0.11*pi*n);
% x2[n] = cos(781*n/2260)
x2 = @(n) cos(781*n/2260);

% Periodicity for x1
[num, den] = rat(0.11/2);  
N1 = den;

% Periodicity for x2
% 781/2261 is an irrational multiple of pi

% Plot x1[n] for one period
figure;
subplot(3, 1, 1);
stem(n(1:N1), x1(1:N1), 'r', 'LineWidth', 1.5);
title('x1[n] = cos(0.11\pi n) for one period');
xlabel('n');
ylabel('x1[n]');
grid on;

% Plot x2[n] 
subplot(3, 1, 2);
stem(n(1:N1), x2(1:N1), 'b', 'LineWidth', 1.5);
title('x2[n] = cos(781n/2260) for 200 samples');
xlabel('n');
ylabel('x2[n]');
grid on;

% Plot x1[n] - x2[n]
large_n = 0:1000;
x1_large = @(large_n) cos(0.11*pi*large_n);
x2_large = @(large_n) cos(781*large_n/2260);
resta = @(large_n) x1_large(large_n) - x2_large(large_n);

subplot(3, 1, 3);
stem(large_n, resta(large_n), 'k', 'LineWidth', 1.5);
title('x1[n] - x2[n]');
xlabel('n');
ylabel('x1[n] - x2[n]');
grid on;

% The difference shows how the two signals diverge over time.
