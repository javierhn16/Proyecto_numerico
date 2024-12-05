
function [t_vals, S_vals, V_vals] = black_scholes_explicit(K, T, sigma, r, M)
    %
    %   Función que calcula el precio de una opción europea mediante el
    %   método explícito en diferencias finitas.
    %
    %   Inputs
    %       K: Precio de ejercicio
    %       T: Tiempo hasta el vencimiento (en años)
    %       r: Tasa libre de riesgo
    %       sigma: Volatilidad
    %       M: Número de pasos en el precio
    %
    %   Outputs
    %       t_vals: Vector con los tiempos discretizados
    %       S_vals: Vector con los precios del subyacente
    %       V_vals: Matriz con los precios de la opción en cada nodo
    %

    Smax = 4 * K; % supuesto teórico

    % Discretización
    ds = Smax / M;
    dt = 0.9 / sigma ^ 2 / M ^ 2;
    n = floor(T / dt) + 1;
    dt = T / n;

    % Malla de puntos
    t_vals = linspace(0, T, n + 1);
    S_vals = linspace(0, Smax, M + 1)'; 
    V_vals = zeros(M + 1, n + 1);
     
    i = (0:M)';
    j = 0:n;

    % Condiciones de frontera
    V_vals(:, n + 1) = max(S_vals - K, 0); 
    V_vals(1, :) = 0; 
    V_vals(M + 1, :) = Smax - K * exp(-r * dt * (n - j)); 

    % Coeficientes del método explícito
    a = 0.5 * dt * ((sigma ^ 2) * i - r) .* i;
    b = 1 - dt * ((sigma ^ 2) * (i .^ 2) + r);
    c = 0.5 * dt * ((sigma ^ 2) * i + r) .* i;

    % Método explícito para calcular los precios hacia atrás en el tiempo
    for j = n:-1:1
        for i = 2:M
            V_vals(i, j) = a(i) * V_vals(i-1, j+1) + b(i) * V_vals(i, j+1) + c(i) * V_vals(i+1, j+1);
        end
    end
end
