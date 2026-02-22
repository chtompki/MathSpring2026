function [L,U] = lucrout(A)
    n = length(A);
    L = zeros(size(A));
    U = zeros(size(A));
    L(:,1) = A(:,1);
    U(:,1) = A(:,1)/L(1,1);
    U(1,1) = 1;
    for k = 2:n
        for j = 2:n
            for i = j:n
                L(i,j) = A(i,j) - dot(L(i,1:j-1),U(1:j-1,j));
            end
            U(k,j) = (A(k,i) - dot(L(k,1:k-1),U(1:k-1,j))) / L(k,k);
        end
    end
end
