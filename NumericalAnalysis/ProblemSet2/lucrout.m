function [L,U] = lucrout(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('A must be square');
    end
    L = zeros(n,n);
    U = eye(n);
    for j = 1:n
        for i = j:n
            L(i,j) = A(i,j) - L(i,1:j-1)*U(1:j-1,j);
        end
        if L(j,j) == 0
            error('Zero pivot at j=%d (need pivoting).', j);
        end
        for k = j+1:n
            U(j,k) = (A(j,k) - L(j,1:j-1)*U(1:j-1,k)) / L(j,j);
        end
    end
end
