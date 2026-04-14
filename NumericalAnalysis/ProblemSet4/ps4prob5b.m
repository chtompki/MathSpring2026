% Forward Euler for y' = sin(t) + y, y(0) = 0 on [0,1]
clear; clc;
f = @(t,y) sin(t) + y;
yexact = @(t) 0.5*(exp(t) - sin(t) - cos(t));

a = 0;
b = 1;
Nh = 2;
h = (b-a)/Nh;

t = zeros(Nh+1,1);
y = zeros(Nh+1,1);

% initial condition
t(1) = 0;
y(1) = 0;

% Forward Euler method
for n = 1:Nh
    t(n+1) = t(n) + h;
    y(n+1) = y(n) + h*f(t(n),y(n));
end

% exact value and error at t=1
y_exact_at_1 = yexact(1);
error_at_1 = abs(y(end) - y_exact_at_1);

fprintf('Forward Euler approximation at t=1: %.10f\n', y(end));
fprintf('Exact value at t=1: %.10f\n', y_exact_at_1);
fprintf('Absolute error at t=1: %.10f\n', error_at_1);

