% Particle Swarm Optimization (PSO)
clear all
close all

% Problem to optimize (minimize), the definition of the objective function
funObj=@(xi) 3*(1-xi(1))^2*exp(-(xi(1)^2)-(xi(2)+1)^2)-10*(xi(1)/5-xi(1)^3-xi(2)^5)* ...
             exp(-xi(1)^2-xi(2)^2)-1/3*exp(-(xi(1)+1)^2-xi(2)^2);

N = 10; % Particle number
d = 2; % Dimensions
lb = [-3 -3]; % Lower limit to search the space
ub = [3 3]; % Upper limit to search the space
k = 0; % Current iteration
kmax = 100; % Maximum numbr of iterations
c1 = 2; % Cognitive constant
c2 = 2; % Social constant

% Initialization of particles and velocity
for i=1:N
    x(i, :) = rand(1, d).*(ub-lb)+lb; % Initialization of particles
    v(i, :) = zeros(1, d); % Initialization of velocities.
end

% Evaluation of the initial particles in the objective function
for i=1:N
    xi = x(i, :); % Extraction of the particle xi
    fx(i, :) = funObj(xi); % Evaluation of the particle xi
end

% Record of the best global particle and the best local particles
[gfit, ind] = min(fx); % Fitness of the best global particle
g = x(ind, :); % Position of the best global particle
fp = fx; % Fitness of the best local particles
p = x; % Position of the best local particles

% Calculation of the search space for graphical purposes
axisx = linspace(min(lb), max(ub), 50); % Solutions for vector d = 1
axisy = axisx; % Solutions for vector d = 2
axisz = []; % Fitness matrix

for i=1:length(axisx)
    for j=1:length(axisy)
        axisz(i, j) = funObj([axisx(i) axisy(j)]);
    end
end

[axisy, axisx] = meshgrid(axisx, axisy); % Computation of the meshgrid

% Itrative process
while k<kmax
    k = k + 1;

    % Draw the search space surface
    figure(1);
    surf(axisx, axisy, axisz);
    hold on;

    % Draw the particles
    % Draw the particles in red color
    plot3(x(:, 1), x(:, 2), fx + 0.1, 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 7);
    % Draw best particles in green color
    plot3(p(:, 1), p(:, 2), fp + 0.1, 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    % Pause to allow viewing the graph
    pause(0.1);
    hold off;
    % Draw the function contour
    figure(2)
    contour(axisx, axisy, axisz, 20);
    hold on;
    % Draw the particles in the contour
    % Draw the particles in red color
    plot(x(:, 1), x(:, 2), 'o', 'MarkerFaceColor', 'm');
    % Draw the particles in gren color
    plot(p(:, 1), p(:, 2), 'o', 'MarkerFaceColor', 'g');
    % Pause to allow viewing the graph
    pause(0.3);
    hold off;

    % Computation of the new velocity for each particle
    for i=1:N
        % Extraction of the particle xi
        xi = x(i, :);
        % Extraction of the local particle pi
        pi = p(i, :);
        % Determination of the new velocity for each particle vi
        v(i, :) = v(i, :) + c1*rand(1, d).*(pi - xi) + c2*rand(1, d).*(g - xi);
    end

    % Determination of the new position of each particle
    x = x + v;
    % Verify that particles do not leave the limits lb and ub
    for i=1:N
        for j=1:d
            if x(i, j) < lb(j)
                x(i, j) = lb(j);
            elseif x(i, j) > ub(j)
                x(i, j) = ub(j);
            end
        end
    end

    % Evaluation of the new particles with the objective function
    for i=1:N
        xi = x(i, :); % Extraction of particle xi
        fx(i, :) = funObj(xi); % Evaluation of particle xi
    end

    % Record of the best global particle and the best local particles
    [gfitkplus1, ind] = min(fx);
    % If a better solution was found, update the global particle
    if gfitkplus1 < gfit
        % Update the fitness of the best global particle
        gfit = gfitkplus1;
        % Update the position of the best global particle
        g = x(ind, :);
    end

    for i=1:N
        % If ay particle is a better solution than the previous one
        % update your best local particle

        if fx(i, :) < fp(i, :)
            % Update the fitnees of the best local particles 
            fp(i, :) = fx(i, :);
            % Update the position of the best local particles
            p(i, :) = x(i, :);
        end
    end

    % Historical record of the best solutions found in each generation
    Evolution(k) = gfit;
end

% End of the iterative process, display of results
figure
% Graoh of the evolutionary process of the PSO
plot(Evolution);

disp(['The best solution: ', num2str(g)]);
disp(['The best fitness: ', num2str(gfit)]);