
%% Características de la opción
K = 13.5;
T = 1.5;
sigma = 0.25;
r = 0.0135;

Smax= 4 * K;
N = 100;
M = 100;

%% Método implícito de diferencias finitas 

[t_implicit, S_implicit, V_implicit] = black_scholes_implicit(K, T, sigma, r, N, M);

figure;
surf(t_implicit, S_implicit, V_implicit', 'EdgeColor', 'none'); 
xlabel('Tiempo (años)');
ylabel('Precio del activo subyacente');
zlabel('Precio de la opción');
title(['Valor de una opción tipo Call Europea usando el método ' ...
    'implícito de diferencias finitas']);
colorbar;
grid on;

% Valores en t = 0
C_implicit = V_implicit(1, :);

%% Método explícito de diferencias finitas

[t_explicit, S_explicit, V_explicit] = black_scholes_explicit(K, T, sigma, r, M);

figure;
surf(t_explicit, S_explicit, V_explicit, 'EdgeColor', 'none'); 
xlabel('Tiempo (años)');
ylabel('Precio del activo subyacente');
zlabel('Precio de la opción');
title(['Valor de una opción tipo Call Europea usando el método ' ...
    'explícito de diferencias finitas']);
colorbar;
grid on;

% Valores en t = 0
C_explicit = V_explicit(:, 1)';

%% Solución exacta en t = 0

C_exact = black_scholes_exact(linspace(0, Smax, M+1), K, T, r, sigma);

%% Comparación soluciones numérica vs solución exacta

figure;
plot(S_implicit, C_exact, 'k-', 'DisplayName', 'Solución exacta', 'LineWidth', 3);
hold on;
plot(S_implicit, C_implicit, 'b-', 'DisplayName', 'Solución método implícito', 'LineWidth', 1);
plot(S_implicit, C_explicit, 'r--', 'DisplayName', 'Solución método explícito', 'LineWidth', 2);
xlabel('Precio del activo subyacente');
ylabel('Precio de la opción');
legend show;
title('Comparación soluciones numéricas vs exacta');

%% Error

error_implicito = abs(C_exact - C_implicit);
error_explicito = abs(C_exact - C_explicit);

figure;
plot(S_implicit, error_implicito, 'r.-', 'DisplayName', 'Error método implícito', 'LineWidth', 1.5, 'MarkerSize', 10);
hold on;
plot(S_implicit, error_explicito, 'b.-', 'DisplayName', 'Error método explícito', 'LineWidth', 1.5, 'MarkerSize', 10);
xlabel('Precio del activo subyacente');
ylabel('Error');
legend show;
title('Error en las aproximaciones numéricas según el método utilizado');
grid on;

%% Malla de puntos

[t_grid, S_grid] = meshgrid(t_implicit, S_implicit);

figure;
plot(t_grid, S_grid, 'k.'); 
xlabel('Tiempo (años)');
ylabel('Precio del activo subyacente');
title('Malla de puntos');
grid on;