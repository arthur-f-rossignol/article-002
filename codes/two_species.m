%% Numerical scheme for the two-species model integration

function dU_dt = two_species(t, U, parameters)

    % definition of parameters
    n     = parameters(1);
    dz    = parameters(2);
    D     = parameters(3);
    N_0   = parameters(4);
    I_0   = parameters(5);
    a_bg  = parameters(6);
    E     = parameters(7);
    r     = parameters(8);
    mu_1  = parameters(9);
    K_1   = parameters(10);
    H_1   = parameters(11);
    m_1   = parameters(12);
    v_1   = parameters(13);
    q_1   = parameters(14);
    a_1   = parameters(15);
    mu_2  = parameters(16);
    K_2   = parameters(17);
    H_2   = parameters(18);
    m_2   = parameters(19);
    v_2   = parameters(20);
    q_2   = parameters(21);
    a_2   = parameters(22);

    % preallocation of vectors
    dA1_dz  = zeros(n, 1);
    dA1_dz2 = zeros(n, 1);
    dA1_dt  = zeros(n, 1);
    G1      = zeros(n, 1);
    dA2_dz  = zeros(n, 1);
    dA2_dz2 = zeros(n, 1);
    dA2_dt  = zeros(n, 1);
    G2      = zeros(n, 1);
    dN_dz2  = zeros(n, 1);
    dN_dt   = zeros(n, 1);
    I       = zeros(n, 1);

    % starting values of A1, A2, N
    A1 = U(1:n);
    A2 = U((n + 1):(2 * n));
    N  = U((2 * n + 1):(3 * n));

    % computation of light intensity
    I(1) = I_0 * exp(- a_bg * 0.5 * dz ...
                     - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
                     - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz);
    for i = 2:n
        S = - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz ...
            - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz;
        for k = 2:(i - 1)
            S = S - a_1 * A1(k) * dz ...
                  - a_2 * A2(k) * dz;
        end
        I(i) = I_0 * exp(S - a_bg * (i - 0.5) * dz ...
                           - a_1 * A1(i) * 0.5 * dz ...
                           - a_2 * A2(i) * 0.5 * dz);
    end

    % 1st spatial derivative of A1 (upwind 3rd-order scheme)
    dA1_dz(3:(n - 1)) = (2 * A1(4:n) ...
                      + 3 * A1(3:(n - 1)) ...
                      - 6 * A1(2:(n - 2)) ...
                      + A1(1:(n - 3))) / (6 * dz); 
    dA1_dz(1)         = (A1(1) + A1(2)) / (2 * dz);
    dA1_dz(2)         = (- A1(1) + 5 * A1(2) + 2 * A1(3)) / (6 * dz) ...
                      - (A1(1) + A1(2)) / (2 * dz);
    dA1_dz(n)         = (A1(n - 2) - 5 * A1(n - 1) - 2 * A1(n)) / (6 * dz);

    % 1st spatial derivative of A2 (upwind 3rd-order scheme)
    dA2_dz(3:(n - 1)) = (2 * A2(4:n) ...
                      + 3 * A2(3:(n - 1)) ...
                      - 6 * A2(2:(n - 2)) ...
                      + A2(1:(n - 3))) / (6 * dz); 
    dA2_dz(1)         = (A2(1) + A2(2)) / (2 * dz);
    dA2_dz(2)         = (- A2(1) + 5 * A2(2) + 2 * A2(3)) / (6 * dz)  ...
                      - (A2(1) + A2(2)) / (2 * dz);
    dA2_dz(n)         = (A2(n - 2) - 5 * A2(n - 1) - 2 * A2(n)) / (6 * dz);

    % 2nd spatial derivative of A1 (symmetrical 2nd-order scheme)
    dA1_dz2(2:(n - 1)) = (A1(3:n) ...
                       - A1(2:(n - 1)) ...
                       - A1(2:(n - 1)) ...
                       + A1(1:(n - 2))) / dz^2;     
    dA1_dz2(1)         = (A1(2) - A1(1)) / dz^2;
    dA1_dz2(n)         = (A1(n - 1) - A1(n)) / dz^2;

    % 2nd spatial derivative of A2 (symmetrical 2nd-order scheme)
    dA2_dz2(2:(n - 1)) = (A2(3:n) ...
                       - A2(2:(n - 1)) ...
                       - A2(2:(n - 1)) ...
                       + A2(1:(n - 2))) / dz^2;    
    dA2_dz2(1)         = (A2(2) - A2(1)) / dz^2;
    dA2_dz2(n)         = (A2(n - 1) - A2(n)) / dz^2;

    % 2nd spatial derivative of N (symmetrical 2nd-order scheme)
    dN_dz2(2:(n - 1)) = (N(3:n) ...
                      - N(2:(n - 1)) ...
                      - N(2:(n - 1)) ...
                      + N(1:(n - 2))) / dz^2; 
    dN_dz2(1)         = (N(2) - N(1)) / dz^2;
    dN_dz2(n)         = E * (N_0 - N(n)) / dz - (N(n) - N(n - 1)) / dz^2;

    % growth rate
    for i = 1:n
        G1(i) = growth_rate(N(i), I(i), mu_1, K_1, H_1);
        G2(i) = growth_rate(N(i), I(i), mu_2, K_2, H_2);
    end

    % time derivatives of A and N
    dA1_dt(1:n) = (G1(1:n) - m_1) .* A1(1:n) ...
                - v_1 * dA1_dz(1:n) ...
                + D * dA1_dz2(1:n);
    dA2_dt(1:n) = (G2(1:n) - m_2) .* A2(1:n) ...
                - v_2 * dA2_dz(1:n) ...
                + D * dA2_dz2(1:n);
    dN_dt(1:n)  = q_1 * (r * m_1 - G1(1:n)) .* A1(1:n) ...
                + q_2 * (r * m_2 - G2(1:n)) .* A2(1:n) ...
                + D * dN_dz2(1:n);
 
    % output
    dU_dt = [dA1_dt; dA2_dt; dN_dt];

end
