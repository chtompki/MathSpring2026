function x = qrSolve(Q, R, b)
    % x = qrSolve(Q, R, b)
    % Q - an n by n matrix that is the left hand side of a QR decomposition
    % R - an n by n matrix that is the right hand side of the QR
    % b - the vector on the right hand side of our equation Xx=b or QRx = b
    % The goal here is to solve our linear set of equations relying upon a
    % QR decomposition. This is equivalent to writing an algorithm
    % for solving an upper triangular equation because our equation is also
    % seen as Q'QRx = Q'b => Rx = Q'b
    % by Rob Tompkins, 20260215
    n = length(b);
    V = Q'*b;
    x = zeros(n,1);
    x(n) = V(n)/R(n,n);
    for k=n-1:-1:1
        x(k) = (V(k) - R(k, k+1:n) * x(k+1:n)) / R(k, k);
    end
end
