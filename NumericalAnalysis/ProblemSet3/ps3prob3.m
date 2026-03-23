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

x = [0 1 2];
y = [2 -1 4];
xbar=1.3;

ybar = nevilleInterpolation(x,y,xbar)

