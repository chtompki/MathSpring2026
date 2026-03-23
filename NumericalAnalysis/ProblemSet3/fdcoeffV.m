function c = fdcoeffV(n, dx)
%FDCOEFFV  Finite difference coefficients for the n-th derivative
%
% Input:
%   n  = derivative order
%   dx = vector of offsets dx_j = x_j - xbar
%
% Output:
%   c  = coefficients such that
%        f^(n)(xbar) ≈ sum_j c(j)*f(xbar + dx(j))

    dx = dx(:);          % make column vector
    k = length(dx);

    if n >= k
        error('Need at least k > n stencil points.');
    end

    A = zeros(k,k);
    b = zeros(k,1);

    % Match Taylor-series moments
    for m = 0:k-1
        A(m+1,:) = (dx').^m;
        if m == n
            b(m+1) = factorial(n);
        end
    end

    c = A\b;
end