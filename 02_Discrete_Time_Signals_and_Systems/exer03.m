%% Exercise 03

% Define the delta function
d = @(n) (n == 0) * 1.0;  % Delta Function
n = -10:10; 

% 1. Characteristic polynomial and roots
% y[n] - 0.4*y[n-1] + 0.2*y[n-2] = x[n] - x[n-1]
% y[n] - 0.4*y[n-1] + 0.2*y[n-2] = 0
A = [1 -0.4 0.2];   % Coefficients y[n]
B = [1 -1];         % Coefficients x[n]

% Calculate roots
roots = roots(A);
disp('Roots of the characteristic polynomial:');
disp(roots);

% 2. Mathematical expression of the impulse response 
z1 = roots(1);
z2 = roots(2);

% Initialize h
h1 = zeros(1, length(n));
h1(n==0) = 1;   % Null initial conditions

% Apply the difference equation to find h[n]
for i = 3:length(n)  
    h1(i) = 0.4 * h1(i-1) - 0.2 * h1(i-2) + d(n(i)) - d(n(i-1));
end

% 3. Impulse response using filter()
h = filter(B, A, d(n));

% 4. Plot and compare both impulse responses
figure;
hold on;
stem(n, h1, 'filled', 'DisplayName', 'h1 (Manual Calculation)');
stem(n, h, 'DisplayName', 'h (Using filter)');
xlabel('n');
ylabel('h[n]');
title('Impulse Response of the System');
legend show;
grid on;
