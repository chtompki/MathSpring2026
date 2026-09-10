function [t,w] = rk4(f,a,b,alpha,N)
%RK4 Classical fourth-order Runge-Kutta method.
%
%   [t,w] = rk4(f,a,b,alpha,N)
%
% approximates the solution of the initial value problem
%
%       y' = f(t,y),    a <= t <= b,
%       y(a) = alpha,
%
% using the classical fourth-order Runge-Kutta method with
% N equally spaced steps.
%
% INPUTS:
%   f      - function handle for the right-hand side f(t,y)
%   a      - left endpoint of the interval
%   b      - right endpoint of the interval
%   alpha  - initial condition (scalar or column vector)
%   N      - number of time steps
%
% OUTPUTS:
%   t      - (N+1)-by-1 vector of time values
%   w      - (N+1)-by-n matrix of approximations, where
%            each row w(i,:) approximates y(t(i))^T
%
% Example for a scalar equation:
%   f = @(t,y) y - t.^2 + 1;
%   [t,w] = rk4(f,0,2,0.5,10);
%
% Example for a system:
%   f = @(t,y) [y(2); -y(1)];
%   [t,w] = rk4(f,0,2*pi,[0;1],100);

    h = (b-a)/N;

    % Store the state internally as a column vector so that the same
    % code works for scalar equations and systems.
    alpha = alpha(:);
    n = length(alpha);

    t = linspace(a,b,N+1).';
    w = zeros(N+1,n);
    w(1,:) = alpha.';

    for i = 1:N
        wi = w(i,:).';

        k1 = h*f(t(i), wi);
        k2 = h*f(t(i) + h/2, wi + k1/2);
        k3 = h*f(t(i) + h/2, wi + k2/2);
        k4 = h*f(t(i) + h,   wi + k3);

        wi_next = wi + (k1 + 2*k2 + 2*k3 + k4)/6;
        w(i+1,:) = wi_next.';
    end
end
