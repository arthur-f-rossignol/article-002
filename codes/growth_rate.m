%% Function computing the growth rate

function g = growth_rate(N, I, mu, K, H)

    % Monod function for nutrient use efficiency
    N_monod = N / (K + N);

    % Monod function for light use efficiency
    I_monod = I / (H + I);

    % Liebig's law of the minimum
    g = mu * min(N_monod, I_monod);
  
end
