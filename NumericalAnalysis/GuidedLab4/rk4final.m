function w = rk4final(f,a,b,alpha,N)
%RK4FINAL Classical fourth-order Runge-Kutta method storing only
%the current solution.
%
%   w = rk4final(f,a,b,alpha,N)
%
% returns the approximation at t = b.

    h = (b-a)/N;
    w = alpha(:);
    t = a;

    for i = 1:N
        k1 = h*f(t,w);
        k2 = h*f(t+h/2,w+k1/2);
        k3 = h*f(t+h/2,w+k2/2);
        k4 = h*f(t+h,w+k3);

        w = w + (k1+2*k2+2*k3+k4)/6;

        % Reconstruct t from the index to reduce drift from repeated
        % floating-point addition.
        t = a + i*h;
    end
end
