clear; clc; close all;

% step size and grid
h = pi/4;
t = 0:h:2*pi;
N = length(t)-1;

% exact solution
tt = linspace(0,2*pi,1000);
y_exact = exp(sin(tt));

% initialize arrays
y_fe = zeros(size(t));   % forward Euler
y_be = zeros(size(t));   % backward Euler

% initial condition
y_fe(1) = 1;
y_be(1) = 1;

% forward Euler
for n = 1:N
    y_fe(n+1) = y_fe(n) + h*y_fe(n)*cos(t(n));
end

% backward Euler
for n = 1:N
    y_be(n+1) = y_be(n) / (1 - h*cos(t(n+1)));
end

% plot solutions
figure;
plot(tt, y_exact, 'k-', 'LineWidth', 2); hold on;
plot(t, y_fe, 'bo--', 'LineWidth', 1.5, 'MarkerSize', 7);
plot(t, y_be, 'rs--', 'LineWidth', 1.5, 'MarkerSize', 7);
grid on;
xlabel('t');
ylabel('y');
legend('Exact: e^{sin(t)}', 'Forward Euler', 'Backward Euler', 'Location', 'best');

% actual forward Euler error at mesh points
err_fe = abs(exp(sin(t)) - y_fe);
max_err_fe = max(err_fe);

% theoretical bound from part (c)
bound_fe = h * exp(2*pi) * (exp(2*pi) - 1);

fprintf('Maximum actual forward Euler error = %.10f\n', max_err_fe);
fprintf('Theoretical forward Euler error bound = %.10f\n', bound_fe);