function [p,absError,relError] = fixpt(a,y,n)
% Input - a, the number to sqaure root
%       - y, the initial guess
%       - n, the number of iterations to run
    p = y; % Initialize p with the initial guess
    absError = zeros(1,n);
    relError = zeros(1,n);
    for k = 1:n
        p_new = 0.5 * (p + a / p); % Update p using the fixed-point iteration formula
        absError(k) = abs(p_new - p); % Calculate absolute error
        relError(k) = absError(k) / abs(p_new); % Calculate relative error
        if absError < 1e-10 % Check for convergence
            break;
        end
        p = p_new; % Update p for the next iteration
    end
end

format long
[p,absError,relError] = fixpt(2,30,6)
loglog(1:6,absError, '-k')
hold on
loglog(1:6,relError,'--k')
hold off
legend("Absolute Error", "Relative Error")