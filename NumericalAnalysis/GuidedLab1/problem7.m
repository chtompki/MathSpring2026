function [p,absError,relError] = fixpt(y,n)
    p = y; % Initialize p with the initial guess
    absError = zeros(1,n);
    relError = zeros(1,n);
    for k = 1:n
        p_new = p^2/(2+p); % Update p using the fixed-point iteration formula
        absError(k) = abs(p_new - 0); % Calculate absolute error
        relError(k) = absError(k) / abs(p_new); % Calculate relative error
        if absError < 1e-10 % Check for convergence
            break;
        end
        p = p_new; % Update p for the next iteration
    end
end

format long
[p,absError,relError] = fixpt(0.344444444444444,10)
loglog(1:10,absError, '-k')
hold on
loglog(1:10,relError,'--k')
hold off
legend("Absolute Error", "Relative Error")