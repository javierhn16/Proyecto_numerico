%% Datos de la opción

S_max = 100;  
E = 50;
T = 1;
sigma = 0.2;
r = 0.05;

%% Malla de diferencias

M = 10;
N = 10;

dS = S_max / M;  
dt = T / N;

S = linspace(0, S_max, M+1)';
tau = linspace(0, T, N+1);

[S_grid, t_grid] = meshgrid(S, tau);

figure;
scatter(t_grid(:), S_grid(:), 10, 'filled');
xlabel('Tiempo (t)');
ylabel('Precio del subyacente (S)');
title('Malla de diferencias');
grid on;

%% Condiciones de frontera

V = zeros(M+1, N+1);

V(:, end) = max(S - E, 0); % t = T, payoff de una opción call
V(1, :) = 0; % S = 0, el valor de la call es 0
V(end, :) = S_max; % debo preguntarle a Álvaro
% V(:, 1) = % acá también