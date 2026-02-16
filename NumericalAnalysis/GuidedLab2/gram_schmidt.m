function Q = gram_schmidt( X )
    % Q = gram_schmidt( X )
    % X is an m by n matrix
    % Q is an m by r matrix with 1 <= r <= n that is linearly independent
    % computes the Gram-Schmidt left decomposition matrix Q
    % By Rob Tompkins, 20260215
    [d,n] = size(X);
    m = min(d,n);
    R = zeros(m,n);
    Q = zeros(d,m);
    for i = 1:m
        v = X(:,i);
        for j = 1:i-1
            R(j,i) = Q(:,j)'*v;
            v = v-R(j,i)*Q(:,j);
        end
        R(i,i) = norm(v);
        Q(:,i) = v/R(i,i);
    end
    R(:,m+1:n) = Q'*X(:,m+1:n);
end