
function C = black_scholes_exact(S, K, T, r, sigma)
    % 
    % Función que calcula el precio exacto de una opción europea tipo call 
    % usando Black-Scholes en el tiempo t = 0.
    %
    % Inputs:
    %   S: Vector de precios del activo
    %   K: Precio de ejercicio
    %   T: Tiempo de vencimiento (años)
    %   r: Tasa libre de riesgo
    %   sigma: Volatilidad
    %
    % Output:
    %   C: Precio de la opción call en t - 0
    %

    d1 = (log(S / K) + (r + 0.5 * sigma^2) * T) ./ (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);

    C = S .* normcdf(d1) - K * exp(-r * T) .* normcdf(d2);
end