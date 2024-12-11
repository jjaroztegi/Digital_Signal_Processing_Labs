%% Exercise 1
clearvars;
close all;
clc;

% 1. Analytically response y[n] to a unit step u[n].

% Transfer function H(z)
% H(z)=(1 - z^-1 / 1 + 0.8*z^-1 + 0.01*z^-2)
H = @(z) (z^2 - z) / (z^2 + 0.8*z + 0.01);

% Input signal X(z) for unit step
% X(z) = 1 / (1 - z^-1)
X = @(z) z / (z - 1);

% Convolution in the z-domain
% Y(z) = X(z) * H(z)
Y = @(z) (z(z^2 - z)) / ((z - 1)*(z^2 + 0.8*z + 0.01));

% Coefficients for the convolution
x1 = [1];              
h1 = [1 -1];           
x2 = [1 -1];           
h2 = [1 0.8 0.01];     

% Compute residues
[res, pk, ks] = residuez(conv(x1, h1), conv(x2, h2));

% Inverse Z-transform for causal response
n = -30:30;
resul = 0;
u = @(n) (n >= 0);       % Unit step function
d = @(n) (n == 0);       % Delta function

% Compute the response y[n]
for i = 1:length(res)
    resul = resul + (res(i) .* pk(i).^n .* u(n));
end

y1 = @(n) resul;

% 2. Plot and compare with filter() function.
% y[n] + 0.8 y[n−1] + 0.001 y[n−2] = 1 x[n] - 1 x[n−1]
B = [1 -1]; % Coeffs that multiply x[n], x[n−1], ...
A = [1 0.8 0.01]; % Coeffs that multiply y[n], y[n−1], ... 
y2 = @(n) filter(B,A,u(n)); % Relaxed => zero initial conditions

figure
stem(n,y1(n),'o')
hold on 
stem(n,y2(n),'x')
title('y[n]')
xlabel('n')
ylabel('y[n]')
legend('manual','filter() ')
grid on

% 3. Initial conditions are set to y[-1] = 1 and y[-2] = 2

a=1;
b=2;
init=filtic(B,A,a,b);
y2_init= @(n) filter(B,A,u(n),init);

figure
stem(n,y2_init(n))
title('y[n] not relaxed')
xlabel('n')
ylabel('y[n]')

% 4. Calculate the impulse response h[n] of the system analytically and using the filter function.

% Coefficients of the quadratic equation
a = A(1);
b = A(2);
c = A(3);

% Calculate the discriminant
D = b^2 - 4*a*c;

% Calculate the roots of the quadratic equation
p1 = (-b + sqrt(D)) / (2 * a);
p2 = (-b - sqrt(D)) / (2 * a);

% We know it is causal thanks to the difference equation
% h[n] = A_k * p_k^n * u[n]
% A_k = (z - p_k) * H(z) / z | evaluated at z = p_k

% Calculate coefficients A1 and A2 based on the roots p1 and p2
A1 = (p1 - 1) / (p1 - p2);
A2 = (p2 - 1) / (p2 - p1);

h = @(n) (A1 .* p1.^n .* u(n) + A2 .* p2.^n .* u(n));

yh = filter(B, A, d(n));

figure
stem(n, h(n), 'o')
hold on
stem(n, yh, 'x')
title('Impulse Response h[n]')
xlabel('n')
ylabel('h[n]')
legend('manual','filter() ')
grid on

% because h[n] is absolutely summable, system is stable











