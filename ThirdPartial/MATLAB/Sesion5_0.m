% Genetic Algorithm Parameters
populationSize = 40; % Number of individuals in the population
chromosomeLength = 2; % Length of each individual's chromosome (number of variables)
crossoverRate = 0.9; % Probability of crossover
mutationRate = 0.002; % Probability of mutation
generations = 40; % Number of generations

% Domain Limits
lowerBound = [0, -2];
upperBound = [3, 2];

% Define the fitness function
fitnessFunction = @(x, y) (4 - 2.1*x^2 + x^4/3)*x^4 + x*y + 4*(y^2 - 1)*y^2;

% Initialization
population = repmat(lowerBound, populationSize, 1) + rand(populationSize, chromosomeLength) .* ...
             repmat((upperBound - lowerBound), populationSize, 1);

% Main Loop
for generation = 1:generations
    % Evaluation
    fitness = zeros(populationSize, 1);

    for i = 1:populationSize
        individual = population(i, :);

        % Extract x and y values from the chromosome
        x = individual(1);
        y = individual(2);

        % Check if individual is within the domain
        if any(individual < lowerBound) || any(individual > upperBound)
            fitness(i) = inf; % Assign high fitness to individuals outside the domain
        else
            fitness(i) = fitnessFunction(x, y);
        end
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

    % Find the best individual and its fitness in the current generation
    [~, bestIdx] = min(fitness);
    bestIndividual = population(bestIdx, :);

    % Plotting the best individual in the current generation
    clf; % Clear the previous plot
    hold on;
    plot(population(:, 1), population(:, 2), 'ko'); % Plot all individuals
    plot(bestIndividual(1), bestIndividual(2), 'ro', 'MarkerSize', 10); % Plot the best individual
    hold off;

    % Set plot properties
    xlim([lowerBound(1), upperBound(1)]);
    ylim([lowerBound(2), upperBound(2)]);
    xlabel('X');
    ylabel('Y');
    title(['Best Individual in Generation ', num2str(generation)]);
    legend('Individuals', 'Best Individual');
    grid on;
    
    pause(0.1); % Pause for visualization
    
    disp(['Generation: ', num2str(generation), ', Best Individual: ', mat2str(bestIndividual)]);
end

% Find the best solution within the domain
fitness = zeros(populationSize, 1);

for i = 1:populationSize
    individual = population(i, :);

    % Extract x and y values from the chromosome
    x = individual(1);
    y = individual(2);

    % Check if individual is within the domain
    if any(individual < lowerBound) || any(individual > upperBound)
        fitness(i) = inf; % Assign high fitness to individuals outside the domain
    else
        fitness(i) = fitnessFunction(x, y);
    end
end

% Find the best individual and fitness
[bestFitness, bestIdx] = min(fitness);
bestIndividual = population(bestIdx, :);

% Plotting the best individual and its fitness
figure;
hold on;
plot(population(:, 1), population(:, 2), 'ko'); % Plot all individuals
plot(bestIndividual(1), bestIndividual(2), 'ro', 'MarkerSize', 10); % Plot the best individual
hold off;

% Set plot properties
xlim([lowerBound(1), upperBound(1)]);
ylim([lowerBound(2), upperBound(2)]);
xlabel('X');
ylabel('Y');
title('Best Individual');
legend('Individuals', 'Best Individual');
grid on;

disp('Best Individual:');
disp(bestIndividual);
disp(['Best Fitness: ', num2str(bestFitness)]);