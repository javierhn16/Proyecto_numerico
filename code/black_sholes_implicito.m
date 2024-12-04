
function [t_vals, S_vals, V_vals] = black_scholes_implicito(K, T, sigma, r, N, M)

Smax = 4 * K;
V_vals = zeros(N + 1, M + 1);

% Discretización
dt = T / N;
dS = Smax / M;

% Malla de puntos
t_vals = 0:dt:T;
S_vals = 0:dS:Smax;

% Condiciones de frontera
V_vals(:,1) = 0;
V_vals(:,end) = Smax - K;
V_vals(end,:) = max(S_vals - K, 0);

% Valores a, b, c
a = @(j) 0.5 * r * j * dt - 0.5 * sigma .^ 2 * j .^2 * dt;
b = @(j) 1 + sigma .^ 2 * j .^ 2 * dt + r * dt;
c = @(j) -0.5 * r * j * dt - 0.5 * sigma .^ 2 * j .^2 * dt;

% Matriz tridiagonal
for i = N:-1:1  
    A = diag(a(2:M - 1), -1) + diag(b(2:M))+ diag(c(1:M - 2), 1);

    v = V_vals(i + 1, 2:M)';
    v(1) = v(1) - a(1) * V_vals(i,1);
    v(end) = v(end) - c(M + 1) * V_vals(i, M + 1);
    
    V_vals(i, 2:M) = A \ v;
    
    % Condición de no arbitraje
    V_vals(i, 2:M) = max(V_vals(i, 2:M), max(S_vals(2:M) - K, 0));
end
end


%% Aplicación

[t_vals, S_vals, V_vals] = black_scholes_implicito(13.5, 1.5, 0.25, 0.0135, 100, 100);

% Gráfica 3d
[S, T] = meshgrid(S_vals, t_vals)

figure;
surf(S, T, V_vals); 
xlabel('Precio del subyacente (S)');
ylabel('Tiempo (t)');
zlabel('Precio de la opción');
title('Superficie de precios de la opción Call Europea');
colorbar;

% Proyección 2d
N=100;
t_indices = [1, round(N/4), round(N/2), round(3*N/4), N]; 

figure;
hold on;
for i = 1:length(t_indices)
    plot(S_vals, V_vals(t_indices(i), :), 'DisplayName', ['t = ', num2str(t_vals(t_indices(i)))]);
end
hold off;
xlabel('Precio del subyacente (S)');
ylabel('Precio de la opción');
title('Evolución del precio de la opción Call Europea');
legend show;
grid on;
