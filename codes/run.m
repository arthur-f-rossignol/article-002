%% Main script for the prediction procedure and the competition simulation

% abiotic parameters
l    = 100;         % maximum depth of the water column 
D    = 0.8;         % eddy diffusion coefficient
a_bg = 0.1;         % background turbidity
r    = 0.5;         % recycling rate
E    = 0.1;         % sediment interface permeability
N_0  = 10;          % nutrient concentration in the sediments
I_0  = 600;         % light intensity at the surface

% algal parameters of S1
mu_1 = 0.04;        % maximum growth rate
K_1  = 0.04;        % half-saturation constant for nutrient uptake
H_1  = 20;          % half-saturation constant for light uptake 
m_1  = 0.01;        % mortality rate
v_1  = 0.01;        % sinking velocity
q_1  = 2 * 1e-10;   % internal nutrient quota
a_1  = 6 * 1e-10;   % absorption coefficient of biomass

% algal parameters of S2
mu_2 = 0.04;        % maximum growth rate
K_2  = 0.12;        % half-saturation constant for nutrient uptake
H_2  = 5;           % half-saturation constant for light uptake 
m_2  = 0.01;        % mortality rate
v_2  = 0.03;        % sinking velocity
q_2  = 5 * 1e-10;   % internal nutrient quota
a_2  = 1 * 1e-10;   % absorption coefficient of biomass

% spatial discretization
dz = 0.25;
n  = floor(l / dz);
Z  = linspace(0, l, n);

% time
t_max = 1e8;

% solver options
options = odeset('nonnegative', 1, 'RelTol', 1e-8, 'AbsTol', 1e-10);

% resource profiles of trivial steady states
N = N_0 * ones(1, n);
I = zeros(n, 1);  
I(1) = I_0 * exp(- a_bg * 0.5 * dz);
for i = 2:n
    I(i) = I_0 * exp(- a_bg * (i - 0.5) * dz);
end

% eigenvalues λ₁* and λ₂*
lambda_1 = principal_eigenvalue(N, I, [n, dz, D, mu_1, K_1, H_1, m_1, v_1]);
lambda_2 = principal_eigenvalue(N, I, [n, dz, D, mu_2, K_2, H_2, m_2, v_2]);

% integration of S1's monoculture
parameters = [n, dz, D, N_0, I_0, a_bg, E, r, mu_1, K_1, H_1, m_1, v_1, q_1, a_1];
U0 = [1e3 * ones(n, 1);
      transpose(linspace(0, N_0, n))];
[t, monoculture_1] = ode15s(@(t, U) one_species(t, U, parameters), ...
                            [0, t_max], ...
                            U0, ...
                            options);  

% biomass and resource profiles from S1's monoculture
A1 = monoculture_1(end, 1:n);
N1 = monoculture_1(end, (n + 1):(2 * n));
I1 = zeros(n, 1);  
I1(1) = I_0 * exp(- a_bg * 0.5 * dz ...
                  - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz);
for i = 2:n
    S = - a_1 * ((3 * A1(1) - A1(2)) / 8 + 3 * A1(1) / 4) * dz;
    for k = 2:(i - 1)
        S = S - a_1 * A1(k) * dz;
    end
    I1(i) = I_0 * exp(S - a_bg * (i - 0.5) * dz ...
                        - a_1 * A1(i) * 0.5 * dz);
end

% eigenvalues Λ₂*
LAMBDA_2 = principal_eigenvalue(N1, I1, [n, dz, D, mu_2, K_2, H_2, m_2, v_2]);

% integration of S2's monoculture
parameters = [n, dz, D, N_0, I_0, a_bg, E, r, mu_2, K_2, H_2, m_2, v_2, q_2, a_2];
U0 = [1e3 * ones(n, 1);
      transpose(linspace(0, N_0, n))];
[t, monoculture_2] = ode15s(@(t, U) one_species(t, U, parameters), ...
                            [0, t_max], ...
                            U0, ...
                            options);  

% biomass and resource profiles from S2's monoculture
A2 = monoculture_2(end, 1:n);
N2 = monoculture_2(end, (n + 1):(2 * n));
I2 = zeros(n, 1);  
I2(1) = I_0 * exp(- a_bg * 0.5 * dz ...
                  - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz);
for i = 2:n
    S = - a_2 * ((3 * A2(1) - A2(2)) / 8 + 3 * A2(1) / 4) * dz;
    for k = 2:(i - 1)
        S = S - a_2 * A2(k) * dz;
    end
    I2(i) = I_0 * exp(S - a_bg * (i - 0.5) * dz ...
                        - a_2 * A2(i) * 0.5 * dz);
end

% eigenvalues Λ₁*
LAMBDA_1 = principal_eigenvalue(N2, I2, [n, dz, D, mu_1, K_1, H_1, m_1, v_1]);

% competition simulation with S1 as resident 
parameters = [n, dz, D, N_0, I_0, a_bg, E, r, ...
              mu_1, K_1, H_1, m_1, v_1, q_1, a_1 ...
              mu_2, K_2, H_2, m_2, v_2, q_2, a_2];
U0 = [1e3 * ones(n, 1);
      1e-3 * ones(n, 1);
      transpose(linspace(0, N_0, n))];
[t, competition_1] = ode15s(@(t, U) two_species(t, U, parameters), ...
                            [0, t_max], ...
                            U0, ...
                            options);

% biomass and resource profiles from competition with S1 as resident 
A1_comp_1 = competition_1(end, 1:n);
A2_comp_1 = competition_1(end, (n + 1):(2 * n));
N_comp_1  = competition_1(end, (2 * n + 1):(3 * n));

% competition simulation with S2 as resident
U0 = [1e-3 * ones(n, 1);
      1e3 * ones(n, 1);
      transpose(linspace(0, N_0, n))];
[t, competition_2] = ode15s(@(t, U) two_species(t, U, parameters), ...
                            [0, t_max], ...
                            U0, ...
                            options);

% biomass and resource profiles from competition with S2 as resident 
A1_comp_2 = competition_2(end, 1:n);
A2_comp_2 = competition_2(end, (n + 1):(2 * n));
N_comp_2  = competition_2(end, (2 * n + 1):(3 * n));
