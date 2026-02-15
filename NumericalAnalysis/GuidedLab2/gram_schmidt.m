function Q = gram_schmidt( X )
    % Q = gram_schmidt( X )
    % X is an m by n matrix
    % Q is an m by r matrix with 1 <= r <= n that is linearly independent
    matrix_size = size(X);
    m = matrix_size(1,1);
    n = matrix_size(1,2);
    if X == zeros(m,n)
        error('There does not exist any type of basis for the zero vector space.');
    elseif n == 1
        Q = X(1:m,1)/norm(X(1:m,1));
    else
        if rank(X) ~= n
            X = basis_col(X);
        end
        matrix_size = size(X);
        m = matrix_size(1,1);
        n = matrix_size(1,2);
        Q = X(1:m,1)/norm(X(1:m,1));
        for i = 2:n
            u = X(1:m,i);
            v = zeros(m,1);
            for j = 1:(i - 1)
                v = v - dot(u,Q(1:m,j))*Q(1:m,j);
            end
            v_ = u + v;
            Q(1:m,i) = v_/norm(v_);
        end
    end
end