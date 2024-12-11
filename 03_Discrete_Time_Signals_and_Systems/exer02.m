% Exercise 2

syms n z;

% Step Function
u = heaviside(n);

% Define the impulse response h[n] = 0.6^n * (u[n + 5] - u[n - 5])
h = 0.6^n * (heaviside(n + 5) - heaviside(n - 5));

% Z-transform of h[n]
H_z = ztrans(h, n, z);

% Define the input signal x[n] = cos(2*pi*n/5) * u[n]
x1 = cos(2 * pi * n / 5) * heaviside(n);
% Z-transform of x1[n]
X1_z = ztrans(x1, n, z);

% Manual convolution for y1[n]
Y1_z = X1_z * H_z;
y1 = simplify(iztrans(Y1_z, z, n));

% Define the input signal x2[n] = cos(2*pi*n/5)
x2 = cos(2 * pi * n / 5);
% Z-transform of x2[n]
X2_z = ztrans(x2, n, z);

% Manual convolution for y2[n]
Y2_z = X2_z * H_z;
y2 = simplify(iztrans(Y2_z, z, n));

% Define range for n
n_values = -100:100;

% Evaluate y1 and y2 for n_values
y1_values = double(subs(y1, n, n_values));
y2_values = double(subs(y2, n, n_values));

% Evaluate h[n] for n_values
h_values = double(subs(h, n, n_values));

% conv() for y1[n] and y2[n]
conv_y1 = conv(h_values, y1_values);
conv_y2 = conv(h_values, y2_values);


% Define new n range for convolution results
n_conv_y1 = -100:(length(conv_y1) - 100 - 1);
n_conv_y2 = -100:(length(conv_y2) - 100 - 1);

figure;
% Plot results
subplot(2, 2, 1);
stem(n_values, y1_values, 'filled');
title('Output for x1[n] = cos(2\pi n/5) * u[n]');
xlabel('n');
ylabel('y1[n]');
grid on;

subplot(2, 2, 2);
stem(n_conv_y1, conv_y1, 'filled');
title('Conv() of y1[n]');
xlabel('n');
ylabel('Conv(y1[n])');
grid on;

subplot(2, 2, 3);
stem(n_values, y2_values, 'filled');
title('Output for x2[n] = cos(2\pi n/5)');
xlabel('n');
ylabel('y2[n]');
grid on; 

subplot(2, 2, 4);
stem(n_conv_y2, conv_y2, 'filled');
title('Conv() of y2[n]');
xlabel('n');
ylabel('Conv(y2[n])');
grid on;

% In this case,  h[n]  is finite within a limited range and decays exponentially, 
% indicating that the system is BIBO stable.

% A system is causal if the output at any time depends only on the present and 
% past input values:
% Since  h[n]  is zero for  n < -5 , it means that the system does not respond to 
% future inputs, making the system causal.

% The system is classified as an FIR (Finite Impulse Response) system since the impulse 
% response  h[n]  is non-zero for a finite number of samples 

% Determine transient and stationary parts
transient_duration = 20; 
stationary_start = transient_duration + 1;

% Plot Transient Response
figure;
subplot(2, 1, 1);
stem(n_values(1:transient_duration), y1_values(1:transient_duration), 'filled');
title('Transient Response of y1[n]');
xlabel('n');
ylabel('y1[n]');
grid on;

% Plot Stationary Response
subplot(2, 1, 2);
stem(n_values(stationary_start:end), y1_values(stationary_start:end), 'filled');
title('Stationary Response of y1[n]');
xlabel('n');
ylabel('y1[n]');
grid on;

% y_1[n]  reflects a transient response that gradually stabilizes.
% y_2[n] exhibits oscillatory behavior because it is not being influenced by the 
% unit step function