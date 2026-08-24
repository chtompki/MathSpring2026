%% Problem 2: second-order convergence of four numerical ODE methods
% Global error is max_n |y_n - y_exact(t_n)| on 0 <= t <= 5.

clear; clc; close all;

T = 5;
Nvalues = [25, 50, 100, 200, 400, 800];
hvalues = T ./ Nvalues;

methodNames = { ...
    'Optimal RK2 (Ralston)', ...
    'Adams-Bashforth 2', ...
    'Adams-Moulton 1', ...
    'Heun predictor-corrector'};

errors = zeros(numel(methodNames), numel(hvalues));
rates  = NaN(size(errors));

for m = 1:numel(methodNames)
    for j = 1:numel(hvalues)
        h = hvalues(j);

        switch m
            case 1
                [t, y] = solveOneStep(@ralstonStep, T, h);
            case 2
                [t, y] = solveAB2(T, h);
            case 3
                [t, y] = solveOneStep(@adamsMoultonStep, T, h);
            case 4
                [t, y] = solveOneStep(@heunStep, T, h);
        end

        errors(m,j) = max(abs(y - yExact(t)));
        if j > 1
            rates(m,j) = log(errors(m,j-1)/errors(m,j))/log(2);
        end
    end
end

%% Print convergence tables
for m = 1:numel(methodNames)
    fprintf('\n%s\n', methodNames{m});
    fprintf('%10s %16s %10s\n', 'h', 'max error', 'rate');
    fprintf('%10.6f %16.8e %10s\n', hvalues(1), errors(m,1), '--');
    for j = 2:numel(hvalues)
        fprintf('%10.6f %16.8e %10.4f\n', ...
            hvalues(j), errors(m,j), rates(m,j));
    end
end

%% Plot every method on the same log-log axes
figure('Color', 'w');
hold on;
styles = {'o-', 's-', 'd-', '^-'};

for m = 1:numel(methodNames)
    loglog(hvalues, errors(m,:), styles{m}, ...
        'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', methodNames{m});
end

% A line proportional to h^2 for comparison.
reference = errors(1,end) .* (hvalues./hvalues(end)).^2;
loglog(hvalues, reference, 'k--', 'LineWidth', 1.3, ...
    'DisplayName', 'reference h^2');

grid on;
xlabel('Step size h');
ylabel('Maximum nodal error');
title('Second-order convergence on 0 \leq t \leq 5');
legend('Location', 'northwest');
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'problem2_convergence_matlab.png', 'Resolution', 200);

%% Local functions
function value = rhs(t, y)
    value = 4*t*sqrt(1 + y.^2)./y;
end

function value = yExact(t)
    value = sqrt((2*t.^2 + sqrt(2)).^2 - 1);
end

function yNext = ralstonStep(t, y, h)
    % Optimal two-stage, second-order Runge-Kutta method.
    k1 = rhs(t, y);
    k2 = rhs(t + 2*h/3, y + 2*h*k1/3);
    yNext = y + h*(k1/4 + 3*k2/4);
end

function yNext = heunStep(t, y, h)
    % Euler predictor followed by one trapezoidal correction.
    predictor = y + h*rhs(t, y);
    yNext = y + (h/2)*(rhs(t, y) + rhs(t + h, predictor));
end

function yNext = adamsMoultonStep(t, y, h)
    % One-step Adams-Moulton (implicit trapezoidal rule).
    % Newton iteration solves
    % z - y - h/2*(f(t,y) + f(t+h,z)) = 0.
    fn = rhs(t, y);
    z = y + h*fn;                  % Euler initial guess
    tolerance = 1e-13;
    maximumIterations = 20;

    for iteration = 1:maximumIterations
        residual = z - y - (h/2)*(fn + rhs(t + h, z));

        % Since d/dz[sqrt(1+z^2)/z] = -1/(z^2*sqrt(1+z^2)),
        % this is the derivative of the residual with respect to z.
        jacobian = 1 + 2*h*(t + h)/(z^2*sqrt(1 + z^2));
        correction = residual/jacobian;
        z = z - correction;

        if abs(correction) <= tolerance*max(1, abs(z))
            yNext = z;
            return;
        end
    end

    error('Newton iteration did not converge.');
end

function [t, y] = solveOneStep(stepFunction, T, h)
    N = round(T/h);
    t = (0:N)*h;
    y = zeros(size(t));
    y(1) = 1;

    for n = 1:N
        y(n+1) = stepFunction(t(n), y(n), h);
    end
end

function [t, y] = solveAB2(T, h)
    N = round(T/h);
    t = (0:N)*h;
    y = zeros(size(t));
    y(1) = 1;

    % A second-order startup is required to retain AB2's global order.
    y(2) = ralstonStep(t(1), y(1), h);

    for n = 2:N
        fn = rhs(t(n), y(n));
        fnMinus1 = rhs(t(n-1), y(n-1));
        y(n+1) = y(n) + h*(3*fn - fnMinus1)/2;
    end
end

%===========OUTPUT=======================
%Optimal RK2 (Ralston)
%         h        max error       rate
%  0.200000   8.92440542e-03         --
%  0.100000   2.08787877e-03     2.0957
%  0.050000   5.09660918e-04     2.0344
%  0.025000   1.26166968e-04     2.0142
%  0.012500   3.14027146e-05     2.0064
%  0.006250   7.83433508e-06     2.0030
%
%Adams-Bashforth 2
%         h        max error       rate
%  0.200000   2.88534227e-02         --
%  0.100000   7.22165611e-03     1.9983
%  0.050000   1.80363844e-03     2.0014
%  0.025000   4.50467567e-04     2.0014
%  0.012500   1.12546757e-04     2.0009
%  0.006250   2.81276299e-05     2.0005
%
%Adams-Moulton 1
%         h        max error       rate
%  0.200000   5.87539115e-03         --
%  0.100000   1.44663578e-03     2.0220
%  0.050000   3.60330581e-04     2.0053
%  0.025000   9.00004500e-05     2.0013
%  0.012500   2.24949879e-05     2.0003
%  0.006250   5.62376734e-06     2.0000
%
%Heun predictor-corrector
%         h        max error       rate
%  0.200000   9.20480095e-03         --
%  0.100000   2.10736388e-03     2.1269
%  0.050000   5.11123795e-04     2.0437
%  0.025000   1.26238053e-04     2.0175
%  0.012500   3.13904354e-05     2.0077
%  0.006250   7.82788897e-06     2.0036
%>>