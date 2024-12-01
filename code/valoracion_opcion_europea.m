%% Datos de la opción

S_max = 100;  
E = 50;
T = 1;
sigma = 0.2;
r = 0.05;

%% Discretización

M = 100;
N = 100;

dS = S_max / M;  
dt = T / N;

%% Malla de puntos 
S = linspace(0, S_max, M+1)';
tau = linspace(0, T, N+1);

[S_grid, t_grid] = meshgrid(S, tau);

figure;
scatter(t_grid(:), S_grid(:), 10, 'filled');
xlabel('Tiempo (t)');
ylabel('Precio del subyacente (S)');
title('Malla de puntos');
grid on;

%% Condiciones de frontera

f = zeros(M+1, N+1);

f(:, end) = max(S - E, 0); % t = T, payoff de una opción call
f(1, :) = 0; % S = 0, el valor de la call es 0
f(end, :) = S_max - E * exp(-r * (T - (0:dt:T))); 

%% Coeficientes a, b, c

a = -0.5 * dt * (sigma^2 * S(2:end-1).^2 / dS^2 - r * S(2:end-1) / dS);
b = 1 + dt * (sigma^2 * S(2:end-1).^2 / dS^2 + r);
c = -0.5 * dt * (sigma^2 * S(2:end-1).^2 / dS^2 + r * S(2:end-1) / dS);

%% Matriz tridiagonal

A = diag(b) + diag(a(2:end), -1) + diag(c(1:end-1), 1);

for i = N:-1:1
    f(2:end-1, i) = A \ f(2:end-1, i+1);
end

%% Resultados

disp('Valores de la opción en t=0:');
disp(table(S, f(:, 1)));

figure;
plot(S, f(:, 1), 'LineWidth', 2);
xlabel('Precio del subyacente (S)');
ylabel('Valor de la opción (V)');
title('Valoración de una Opción Call Europea usando Diferencias Finitas');
grid on;