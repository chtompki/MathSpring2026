function L = cholesky(A)
    n = size(A,1);
    if size(A,2) ~= n
        error('A must be square.');
    end
    if norm(A - A.', 'fro') > 1e-12
        error('A must be symmetric (within tolerance).');
    end

    L = zeros(n,n);

    for i = 1:n
        for j = 1:i
            s = 0;
            if j > 1
                s = L(i,1:j-1) * L(j,1:j-1).';
            end

            if i == j
                t = A(i,i) - s;
                if t <= 0
                    error('A is not SPD: encountered non-positive pivot at i=%d.', i);
                end
                L(i,i) = sqrt(t);
            else
                L(i,j) = (A(i,j) - s) / L(j,j);
            end
        end
    end
end