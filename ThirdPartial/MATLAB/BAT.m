function [best, fmin, N_iter] = bat(n, A, r)
    % Parametetros
    n = 40;   % Tamaño de poblacion (10 to 25)
    A = 0.45; % Sonoridad (constante o decreciente)
    r = 0.5; % Tasa de pulso (constante o decreciente)
    
    % El rango de frecuencia determina la escala
    fmin = 0; % Frecuencia mínima
    fmax = 2; % Frecuencia máxima
    
    % Dimensión de las variables de búsqueda
    d = 2;
    
    % Inicializando matrices
    f = zeros(n, 1); % Frecuencia
    v = zeros(n, d); % Velocidades
    
    % Inicializando la población/soluciones
    sol = zeros(n, d);
    Fitness = zeros(n, 1);
    
    for i = 1:n
        sol(i, :) = randn(1, d);
        Fitness(i) = Fun(sol(i, :));
    end
    
    % Encontrando la mejor solución actual
    [fmin, idx] = min(Fitness);
    best = sol(idx, :);
    disp(['Aptitud:  ', num2str(fmin)]);
    
    % Iniciar iteraciones - Algoritmo Bat
    N_iter = 0;
    
    % Definir el rango para las variables x e y
    x = linspace(-2*pi(), 2*pi(), 100);
    y = linspace(-2*pi(), 2*pi(), 100);
    
    % Crear una cuadrícula de puntos
    [X, Y] = meshgrid(x, y);
    
    % Evaluar la función en cada punto de la cuadrícula
    Z = Fun([X(:), Y(:)]);
    Z = reshape(Z, size(X));
    
    % Figura para el gráfico 2D
    figure(1);
    hold on;
    
    % Graficar el contorno
    contour(X, Y, Z, 20, 'LineWidth', 1.5);
    title('Proceso de Optimización - Contorno');
    
    % Figura para el gráfico 3D
    figure(2);
    hold on;
    
    % Graficar la superficie 3D
    surf(X, Y, Z, 'FaceAlpha', 1);
    title('Proceso de Optimización - 3D');
    
    for k = 1:50
        for i = 1:n
            f(i) = fmin + (fmin - fmax) * rand;
            v(i, :) = v(i, :) + (sol(i, :) - best) * f(i);
            S(i, :) = sol(i, :) + v(i, :);
    
            % Tasa de pulso
            if rand > r
                S(i, :) = best + 0.01 * randn(1, d);
            end
    
            % Verificar si los puntos están dentro del rango deseado
            if (S(i, 1) >= min(x)) && (S(i, 1) <= max(x)) && (S(i, 2) >= min(y)) && (S(i, 2) <= max(y))
                % Gráfico de dispersión
                scatter3(S(i, 1), S(i, 2), Fun(S(i, :)), 'r');
            end
    
            % Evaluar las nuevas soluciones
            Fnew = Fun(S(i, :));
    
            % Si la solución mejora o no es demasiado ruidosa
            if (Fnew <= Fitness(i)) && (rand < A)
                sol(i, :) = S(i, :);
                Fitness(i) = Fnew;
            end
    
            % Actualizar la mejor solución actual
            if Fnew <= fmin
                best = S(i, :);
                fmin = Fnew;
            end
    
            pause(0.01);
    
            N_iter = N_iter + 1;
        end
    end
    
    scatter3(best(1), best(2), fmin, 'g', 'filled');
    hold off;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Proceso de Optimización');
    
    % Salida/mostrar
    disp(['Número de evaluaciones: ', num2str(N_iter)]);
    disp(['Mejor = ', num2str(best), ', fmin = ', num2str(fmin)]);
end

% Función objetivo: z = x.^2 + y.^2 + 25*(sin(x).^2 + sin(y).^2)
function z = Fun(u)
    x = u(:, 1);
    y = u(:, 2);
    z = x.^2 + y.^2 + 25*(sin(x).^2 + sin(y).^2);
end