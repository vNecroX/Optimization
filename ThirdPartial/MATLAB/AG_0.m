% Definición de la función objetivo
funcionObjetivo = @(x, y) (4-2.1*x.^2+x.^4/3)+x.^2 + x.*y + 4+(y.^2-1).*y.^2;

% Parámetros del algoritmo genético
numVariables = 2; % Número de variables de diseño
numIndividuos = 50; % Tamaño de la población
numGeneraciones = 100; % Número de generaciones
probabilidadCruce = 0.8; % Probabilidad de cruce
probabilidadMutacion = 0.01; % Probabilidad de mutación

% Rango de las variables de diseño
rangoX = [-5, 5];
rangoY = [-5, 5];

% Inicialización de la población
poblacion = rand(numIndividuos, numVariables);
poblacion(:, 1) = poblacion(:, 1) * (rangoX(2) - rangoX(1)) + rangoX(1); % Escalamiento de x
poblacion(:, 2) = poblacion(:, 2) * (rangoY(2) - rangoY(1)) + rangoY(1); % Escalamiento de y

% Ciclo de generaciones
for generacion = 1:numGeneraciones
    % Evaluación de la aptitud de cada individuo en la población
    aptitud = funcionObjetivo(poblacion(:, 1), poblacion(:, 2));
    
    % Selección de padres mediante torneo binario
    indicesPadres = seleccionTorneoBinario(aptitud, numIndividuos);
    
    % Operador de cruce
    descendientes = cruzarPoblacion(poblacion(indicesPadres, :), probabilidadCruce);
    
    % Operador de mutación
    descendientesMutados = mutarPoblacion(descendientes, probabilidadMutacion, rangoX, rangoY);
    
    % Reemplazo de la población anterior con los descendientes mutados
    poblacion = descendientesMutados;
    
    % Imprimir el mejor individuo de cada generación
    [~, indiceMejor] = min(aptitud);
    mejorIndividuo = poblacion(indiceMejor, :);
    mejorAptitud = aptitud(indiceMejor);
    fprintf('Generación %d - Mejor Individuo: %s - Aptitud: %.4f\n', generacion, num2str(mejorIndividuo), mejorAptitud);
end
% Generar una malla de puntos en el dominio
[x, y] = meshgrid(linspace(rangoX(1), rangoX(2), 100), linspace(rangoY(1), rangoY(2), 100));

% Evaluar la función objetivo en la malla de puntos
z = funcionObjetivo(x, y);

% Graficar la función objetivo en 3D
figure;
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('f(x, y)');
title('Función objetivo');

hold on;

% Punto óptimo encontrado por el algoritmo genético
mejorPuntoX = mejorIndividuo(1);
mejorPuntoY = mejorIndividuo(2);
mejorAptitud = funcionObjetivo(mejorPuntoX, mejorPuntoY);

% Graficar el mínimo global encontrado
scatter3(mejorPuntoX, mejorPuntoY, mejorAptitud, 'ro', 'LineWidth', 2);

hold off;

% Generar una malla de puntos en el dominio
[x, y] = meshgrid(linspace(rangoX(1), rangoX(2), 100), linspace(rangoY(1), rangoY(2), 100));

% Evaluar la función objetivo en la malla de puntos
z = funcionObjetivo(x, y);

% Graficar el contorno de la función objetivo
figure;
contour(x, y, z, 50);
xlabel('x');
ylabel('y');
title('Función objetivo (Contorno)');

hold on;

% Punto óptimo encontrado por el algoritmo genético
mejorPuntoX = mejorIndividuo(1);
mejorPuntoY = mejorIndividuo(2);
mejorAptitud = funcionObjetivo(mejorPuntoX, mejorPuntoY);

% Graficar el mínimo global encontrado
scatter(mejorPuntoX, mejorPuntoY, 'ro', 'LineWidth', 2);

hold off;
function indicesPadres = seleccionTorneoBinario(aptitud, numPadres)
    indicesPadres = zeros(numPadres, 1);
    for i = 1:numPadres
        indice1 = randi(length(aptitud));
        indice2 = randi(length(aptitud));
        if aptitud(indice1) < aptitud(indice2)
            indicesPadres(i) = indice1;
        else
            indicesPadres(i) = indice2;
        end
    end
end

% Operador de cruce
function descendientes = cruzarPoblacion(padres, probabilidadCruce)
    [numDescendientes, numVariables] = size(padres);
    descendientes = padres;
    for i = 1:numDescendientes-1
        if rand < probabilidadCruce
            puntoCruce = randi(numVariables-1) + 1;
            descendientes(i,:) = [padres(i, 1:puntoCruce), padres(i+1, puntoCruce+1:end)];
            descendientes(i+1,:) = [padres(i+1, 1:puntoCruce), padres(i, puntoCruce+1:end)];
        end
    end
end

% Operador de mutación
function descendientesMutados = mutarPoblacion(descendientes, probabilidadMutacion, rangoX, rangoY)
    [numDescendientes, numVariables] = size(descendientes);
    descendientesMutados = descendientes;
    for i = 1:numDescendientes
        for j = 1:numVariables
            if rand < probabilidadMutacion
                if j == 1 % Mutación para la variable x
                    descendientesMutados(i, j) = rand * (rangoX(2) - rangoX(1)) + rangoX(1);
                else % Mutación para la variable y
                    descendientesMutados(i, j) = rand * (rangoY(2) - rangoY(1)) + rangoY(1);
                end
            end
        end
    end
        subplot(2,1,1); plot(descendientesMutados); title ('1');
    subplot(2,1,2); plot(descendientes); title ('2');
end