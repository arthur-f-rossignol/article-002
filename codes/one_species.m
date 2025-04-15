%% Numerical scheme for the one-species model integration

function dU_dt = one_species(t, U, parameters)

    % definition of parameters
    n    = parameters(1);
    dz   = parameters(2);
    D    = parameters(3);
    N_0  = parameters(4);
    I_0  = parameters(5);
    a_bg = parameters(6);
    E    = parameters(7);
    r    = parameters(8);
    mu   = parameters(9);
    K    = parameters(10);
    H    = parameters(11);
    m    = parameters(12);
    v    = parameters(13);
    q    = parameters(14);
    a    = parameters(15);
    
    % preallocation of vectors
    dA_dz   = zeros(n, 1);
    dA_dz2  = zeros(n, 1);
    dA_dt   = zeros(n, 1);
    G       = zeros(n, 1);
    dN_dz2  = zeros(n, 1);
    dN_dt   = zeros(n, 1);
    I       = zeros(n, 1);

    % starting values of A and N
    A = U(1:n);
    N = U((n + 1):(2 * n));

    % light intensity
    I(1) = I_0 * exp(- a_bg * 0.5 * dz ...
                     - a * ((3 * A(1) - A(2)) / 8 + 3 * A(1) / 4) * dz);
    for i = 2:n
        S = - a * ((3 * A(1) - A(2)) / 8 + 3 * A(1) / 4) * dz;
        for k = 2:(i - 1)
            S = S - a * A(k) * dz;
        end
        I(i) = I_0 * exp(S - a_bg * (i - 0.5) * dz ...
                           - a * A(i) * 0.5 * dz);
    end

    % 1st spatial derivative of A (upwind 3rd-order scheme)
    dA_dz(3:(n - 1)) = (2 * A(4:n) ...
                     + 3 * A(3:(n - 1)) ...
                     - 6 * A(2:(n - 2)) ...
                     + A(1:(n - 3))) / (6 * dz); 
    dA_dz(1)         = (A(1) + A(2)) / (2 * dz);
    dA_dz(2)         = (- A(1) + 5 * A(2) + 2 * A(3)) / (6 * dz) ...
                     - (A(1) + A(2)) / (2 * dz);
    dA_dz(n)         = (A(n - 2) - 5 * A(n - 1) - 2 * A(n)) / (6 * dz);

    % 2nd spatial derivative of A (symmetrical 2nd-order scheme)
    dA_dz2(2:(n - 1)) = (A(3:n) ...
                       - A(2:(n - 1)) ...
                       - A(2:(n - 1)) ...
                       + A(1:(n - 2))) / dz^2;     
    dA_dz2(1)         = (A(2) - A(1)) / dz^2;
    dA_dz2(n)         = (A(n - 1) - A(n)) / dz^2;
    
    % 2nd spatial derivative of N (symmetrical 2nd-order scheme)
    dN_dz2(2:(n - 1)) = (N(3:n) ...
                      - N(2:(n - 1)) ...
                      - N(2:(n - 1)) ...
                      + N(1:(n - 2))) / dz^2; 
    dN_dz2(1)         = (N(2) - N(1)) / dz^2;
    dN_dz2(n)         = E * (N_0 - N(n)) / dz - (N(n) - N(n - 1)) / dz^2;

    % growth rate
    for i = 1:n
        G(i) = growth_rate(N(i), I(i), mu, K, H);
    end

    % time derivatives of A and N
    dA_dt(1:n) = (G(1:n) - m) .* A(1:n) ...
               - v * dA_dz(1:n) ...
               + D * dA_dz2(1:n);
    dN_dt(1:n) = q * (r * m - G(1:n)) .* A(1:n) ...
               + D * dN_dz2(1:n);
 
    % output
    dU_dt = [dA_dt; dN_dt];

end
