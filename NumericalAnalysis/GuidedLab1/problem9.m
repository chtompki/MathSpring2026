function [p,absError,relError] = fixpt(a,y,n)
% Input - a, the number to sqaure root
%       - y, the initial guess
%       - n, the number of iterations to run
    p = y; % Initialize p with the initial guess
    absError = zeros(1,n);
    relError = zeros(1,n);
    for k = 1:n
        p_new = p + p^2 - a; % Update p using the fixed-point iteration formula
        absError(k) = abs(p_new - sqrt(a)); % Calculate absolute error
        relError(k) = absError(k) / abs(p_new); % Calculate relative error
        if absError < 1e-10 % Check for convergence
            break;
        end
        p = p_new; % Update p for the next iteration
    end
end

format long
[p,absError,relError] = fixpt(2,1.5,6)
loglog(1:6,absError, '-k')
hold on
loglog(1:6,relError,'--k')

[p,absError18,relError18] = fixpt(2,1.8,6)
loglog(1:6,absError18, '-r')
loglog(1:6,relError18,'--r')

format long
[p,absError12,relError12] = fixpt(2,1.2,6)
loglog(1:6,absError12, '-b')
loglog(1:6,relError12,'--b')
hold off
legend("Absolute Error 1.5", "Relative Error 1.5", "Absolute Error 1.8", "Relative Error 1.8", "Absolute Error 1.2", "Relative Error 1.2")`