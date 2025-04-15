%% Function computing the principal eigenvalue estimate

function p_eig = principal_eigenvalue(N, I, parameters)

    % definition of parameters
    n  = parameters(1);
    dz = parameters(2);
    D  = parameters(3);
    mu = parameters(4);
    K  = parameters(5);
    H  = parameters(6);
    m  = parameters(7);
    v  = parameters(8);

    % preallocation of matrices
    P = zeros(n, n);
    Q = zeros(n, n);

    % boundary conditions
    alpha = 2 - (8 * D) / (4 * D + v * dz);
    beta  = 2 - (8 * D) / (4 * D - v * dz);

    % matrix P
    for i = 1:n
        P(i, i) = 2 * D;
    end
    P(1, 1) = (alpha + 1) * D;
    P(n, n) = (beta + 1) * D;
    for i = 1:(n - 1)
        P(i + 1, i) = - D;
        P(i, i + 1) = - D;
    end
    P = P / dz^2;

    % matrix Q
    for i = 1:n
        Q(i, i) = m + (v^2 / (4 * D)) - growth_rate(N(i), I(i), mu, K, H);
    end

    % spectrum and principal eigenvalue
    M = P + Q;
    spectrum = - eig(M);
    p_eig = max(spectrum);

end
