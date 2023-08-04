% Genetic Algorithm Parameters
populationSize = 40; % Number of individuals in the population
chromosomeLength = 5; % Length of each individual's chromosome (number of variables)
crossoverRate = 0.9; % Probability of crossover
mutationRate = 0.002; % Probability of mutation
generations = 40; % Number of generations

% Domain Limits
lowerBound = -5;
upperBound = 5;

% Fitness Function
fitnessFunction = @(individual) stochasticFitness(individual);

% Initialization
population = lowerBound + (upperBound - lowerBound) * rand(populationSize, chromosomeLength);

% Main Loop
for generation = 1:generations
    % Evaluation
    fitness = arrayfun(fitnessFunction, population);
    
    % Selection
    fitnessSum = sum(fitness);
    selectionProbability = fitness / fitnessSum;
    cumulativeProbability = cumsum(selectionProbability);
    parentsIdx = zeros(populationSize, 1);
    
    for i = 1:populationSize
        r = rand();
        idx = find(cumulativeProbability >= r, 1);
        parentsIdx(i) = idx;
    end
    
    % Crossover
    offspring = zeros(size(population));
    
    for i = 1:2:populationSize
        parent1 = population(parentsIdx(i), :);
        parent2 = population(parentsIdx(i+1), :);
        
        % Crossover Point
        crossoverPoint = randi([1, chromosomeLength-1]);
        
        % Create Offspring
        offspring(i, :) = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
        offspring(i+1, :) = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
    end
    
    % Mutation
    for i = 1:populationSize
        if rand() < mutationRate
            % Select a random gene
            geneIdx = randi(chromosomeLength);
            
            % Apply mutation
            offspring(i, geneIdx) = lowerBound + (upperBound - lowerBound) * rand();
        end
    end
    
    % Replacement
    population = offspring;

    % Find the best individual and its fitness in the current generation
    [bestFitness, bestIdx] = min(fitness);
    bestIndividual = population(bestIdx, :);

    % Plotting the best individual in the current generation
    plot(bestIndividual, 'ro', 'MarkerSize', 10);
    xlabel('Variable Index');
    ylabel('Variable Value');
    title(['Best Individual - Generation ', num2str(generation)]);
    % Update the plot
    drawnow;

    disp(['Generation: ', num2str(generation), ', Best Individual: ', mat2str(bestIndividual)]);
end

hold off;

% Set the axis labels and title
xlabel('Variable Index');
ylabel('Variable Value');
title('Best Individual');

% Display the final best individual and its fitness
disp('Best Individual:');
disp(bestIndividual);
disp(['Best Fitness: ', num2str(bestFitness)]);

% Stochastic Fitness Function
function fitness = stochasticFitness(individual)
    n = length(individual);
    fitness = sum(abs(individual).^((1:n)+1)) + randn(1) * 0.1;
end