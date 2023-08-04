% Genetic Algorithm Parameters
populationSize = 40; % Number of individuals in the population
chromosomeLength = 5; % Length of each individual's chromosome (number of variables)
crossoverRate = 0.9; % Probability of crossover
mutationRate = 0.002; % Probability of mutation
generations = 40; % Number of generations

% Domain Limits
lowerBound = -1;
upperBound = 1;

% Fitness Function
fitnessFunction = @(individual) sum(abs(individual).^((1:length(individual))+1));

% Initialization
population = lowerBound + (upperBound - lowerBound) * rand(populationSize, chromosomeLength);

% Main Loop
for generation = 1:generations
    % Evaluation
    fitness = zeros(populationSize, 1);
    for i = 1:populationSize
        fitness(i) = fitnessFunction(population(i, :));
    end
    
    % Selection
    selectionProbability = fitness / sum(fitness);
    cumulativeProbability = cumsum(selectionProbability);
    parents = zeros(size(population));
    for i = 1:populationSize
        randomNum = rand;
        idx = find(cumulativeProbability >= randomNum, 1);
        if isempty(idx)
            idx = populationSize; % Select the last individual if no valid individual is found
        end
        parents(i, :) = population(idx, :);
    end
    
    % Crossover
    crossoverPoints = randi(chromosomeLength-1, populationSize/2, 1);
    crossoverMask = rand(populationSize/2, 1) < crossoverRate;
    offspring = zeros(size(population));
    for i = 1:2:populationSize
        parent1 = parents(i, :);
        parent2 = parents(i+1, :);
        if crossoverMask((i+1)/2)
            crossoverPoint = crossoverPoints((i+1)/2);
            offspring(i, :) = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
            offspring(i+1, :) = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
        else
            offspring(i, :) = parent1;
            offspring(i+1, :) = parent2;
        end
    end
    
    % Mutation
    mutationMask = rand(size(offspring)) < mutationRate;
    mutatedOffspring = offspring;
    mutatedOffspring(mutationMask) = mutatedOffspring(mutationMask) + randn(sum(mutationMask(:)), 1);
    
    % Replacement
    population = mutatedOffspring;
    
    % Display the best individual in the current generation
    [~, bestIdx] = min(fitness);
    bestIndividual = population(bestIdx, :);

    % Plotting the best individual in the current generation
    clf; % Clear the previous plot
    plot(population(:, 1), population(:, 2), 'ko'); % Plot all individuals
    hold on;
    plot(bestIndividual(1), bestIndividual(2), 'ro', 'MarkerSize', 10); % Plot the best individual
    hold off;
    
    % Set plot properties
    xlim([lowerBound, upperBound]);
    ylim([lowerBound, upperBound]);
    xlabel('X');
    ylabel('Y');
    title(['Best Individual in Generation ', num2str(generation)]);
    legend('Individuals', 'Best Individual');
    grid on;
    
    pause(0.1); % Pause for visualization

    disp(['Generation: ', num2str(generation), ', Best Individual: ', num2str(bestIndividual)]);
end

% Find the best individual and its fitness
[bestFitness, bestIdx] = min(fitness);
bestIndividual = population(bestIdx, :);

% Set plot properties
figure(1);
xlim([lowerBound, upperBound]);
ylim([lowerBound, upperBound]);
xlabel('X');
ylabel('Y');
title('Best Individual and Contour Plot');
legend('Fitness Contour', 'Best Individual');
grid on;

% Display the best individual and its fitness
disp('Best Individual:');
disp(bestIndividual);
disp(['Best Fitness: ', num2str(bestFitness)]);

% Plotting the best individual and its fitness
figure(2);
x = linspace(lowerBound, upperBound, 100); % Generate x values
y = linspace(lowerBound, upperBound, 100); % Generate y values
[X, Y] = meshgrid(x, y); % Create a grid of points
Z = zeros(size(X)); % Initialize the fitness values for each point
for i = 1:numel(X)
    Z(i) = fitnessFunction([X(i), Y(i)]); % Evaluate fitness for each point
end
contour(X, Y, Z, 'LineWidth', 1.5); % Plot contour
hold on;
plot(bestIndividual(1), bestIndividual(2), 'ro', 'MarkerSize', 10); % Plot the best individual
hold off;

% Plotting the 3D surface
figure(3);
surf(X, Y, Z);
hold on;
plot3(bestIndividual(1), bestIndividual(2), bestFitness, 'ro', 'MarkerSize', 10); % Plot the best individual
hold off;
xlabel('X');
ylabel('Y');
zlabel('Fitness');
title('Fitness Landscape');