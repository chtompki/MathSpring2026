function ystar = naturalCubicSplines(x,y,xstar)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Natural cubic spline interpolator
    %
    % Inputs:
    %   x     - vector of x data points
    %   y     - vector of y data points
    %   xstar - point(s) where spline is to be evaluated
    %
    % Output:
    %   ystar - spline value(s) at xstar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(x) < 2 | length(y) < 2
        error('Input vectors are not long enough, please make them of length two or more and be of the same size')
    end
    if length(x) ~= length(y)
        error('Input vectors are not the same length')
    end
    
    % Make sure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % n = number of subintervals
    n = length(x) - 1;
    
    % Initialize A and rvec
    A = zeros(n+1,n+1);
    rvec = zeros(n+1,1);
    
    % Calculate vector of h values
    h = diff(x);
    
    % Natural spline boundary conditions
    A(1,1) = 1;
    A(n+1,n+1) = 1;
    
    % Populate rest of matrix A and RHS vector
    for i = 2:n
        A(i,i-1) = h(i-1);
        A(i,i)   = 2*(h(i-1) + h(i));
        A(i,i+1) = h(i);
        rvec(i)  = 3*((y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1));
    end
    
    % Solve for c coefficients
    c = A \ rvec;
    
    % Compute a coefficients
    a = zeros(n+1,1);
    for j = 1:n+1
        a(j) = y(j);
    end
    
    % Compute b and d coefficients
    b = zeros(n,1);
    d = zeros(n,1);
    
    for j = 1:n
        d(j) = (c(j+1) - c(j)) / (3*h(j));
        b(j) = (a(j+1) - a(j))/h(j) - (h(j)/3)*(2*c(j) + c(j+1));
    end
    
    % Build piecewise polynomial coefficients
    coeffs = [];
    for k = 1:n
        coeffs(k,:) = [d(k) c(k) b(k) a(k)];
    end
    
    breaks = x.';
    pp = mkpp(breaks, coeffs);
    ystar = ppval(pp, xstar);
end