
function [t_vals, S_vals, V_vals] = black_scholes_implicit(K, T, sigma, r, N, M)
    % 
    %   Función que calcula el precio de una opción europea tipo call
    %   utilizando el método implícito de diferencias finitas para la
    %   ecuación de Black-Scholes.
    %
    %   Inputs
    %       K: Precio de ejercicio de la opción.
    %       T: Tiempo hasta el vencimiento de la opción (en años).
    %       sigma: Volatilidad del activo subyacente.
    %       r: Tasa de interés libre de riesgo.
    %       N: Número de subdivisiones en el tiempo.
    %       M: Número de subdivisiones en el precio del subyacente.
    %
    %   Outputs
    %       t_vals: Vector de los valores de tiempo discretizados (0 a T).
    %       S_vals: Vector de precios del subyacente discretizados (0 a Smax).
    %       V_vals: Matriz que contiene el precio de la opción en cada punto
    %               de la malla tiempo-precio.
    %

    Smax = 4 * K; % supuesto teórico
    V_vals = zeros(N + 1, M + 1);
    
    % Discretización
    dt = T / N;
    dS = Smax / M;
    
    % Malla de puntos
    t_vals = 0:dt:T;
    S_vals = 0:dS:Smax;
    
    % Condiciones de frontera
    V_vals(:, 1) = 0;
    V_vals(:, end) = Smax - K * exp(-r * (T - t_vals));
    V_vals(end, :) = max(S_vals - K, 0);
    
    % Coeficientes del método implícito
    a = @(j) 0.5 * r * j * dt - 0.5 * sigma .^ 2 * j .^2 * dt;
    b = @(j) 1 + sigma .^ 2 * j .^ 2 * dt + r * dt;
    c = @(j) -0.5 * r * j * dt - 0.5 * sigma .^ 2 * j .^2 * dt;
    
    % Matriz tridiagonal
    A = diag(a(2:(M - 1)), -1) + diag(b(1:(M - 1)))+ diag(c(1:(M - 2)), 1);
    
    % Factorización LU de A
    [L, U] = lu(A);

    % Método implícito para calcular los valores hacia atrás en el tiempo
    for i = N:-1:1
        v = V_vals(i + 1, 2:M)';
        v(1) = v(1) - a(2) * V_vals(i, 1);
        v(end) = v(end) - c(M-1) * V_vals(i, M + 1); 

        y = L \ v;  
        V_vals(i, 2:M) = U \ y;

        % Condición de no arbitraje
        V_vals(i, 2:M) = max(V_vals(i, 2:M), S_vals(2:M) - K);
    end
end


