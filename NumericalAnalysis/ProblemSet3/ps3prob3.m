function p = nevilleInterpolation(x, y, xbar)
    n = length(x);
    Q = zeros(n, n);

    % First column
    for i = 1:n
        Q(i,1) = y(i);
    end

    % Build Neville table
    for j = 2:n
        for i = 1:n-j+1
            Q(i,j) = ((xbar - x(i)) * Q(i+1,j-1) - (xbar - x(i+j-1)) * Q(i,j-1)) / (x(i+j-1) - x(i));
        end
    end

    p = Q(1,n);

    % Optional: display the table
    disp('Neville table:');
    disp(Q);
end

function p = newtonDividedDiff(x, y, xbar)
    n = length(x);
    DD = zeros(n, n);

    % First column: function values
    DD(:,1) = y(:);

    % Build divided difference table
    for j = 2:n
        for i = 1:n-j+1
            DD(i,j) = (DD(i+1,j-1) - DD(i,j-1)) / (x(i+j-1) - x(i));
        end
    end

    % Evaluate Newton polynomial at xbar
    p = DD(1,1);
    prodterm = 1;

    for j = 2:n
        prodterm = prodterm * (xbar - x(j-1));
        p = p + DD(1,j) * prodterm;
    end

    % Display divided difference table if desired
    disp('Divided difference table:')
    disp(DD)
end

x = [0 1 2];
y = [2 -1 4];
xbar=1.3;

ybar = nevilleInterpolation(x,y,xbar)
ybarNewton = newtonDividedDiff(x, y, xbar)