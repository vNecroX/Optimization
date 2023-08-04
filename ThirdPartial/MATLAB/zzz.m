% ALGORITMO GENETICO "SIMPLE"

function [bestSol, bestFun, count] = gasimple(funstr)
    global solNew sol pop popNew fitness fitOld f range;

    if nargin < 1
        funstr = '-cos(x)*exp(-(x-3.1415926)^2)';
    end

    range = [-10, 10];

    f = vectorize(inline(funstr));
    rand('state', 0'); % reset the random generator
    popsize = 20;
    maxGen = 100;
    count = 0;
    nsite = 2; % number of mutation
    pc = 0.95;
    pm = 0.05;
    nsbit = 16;
    
    popNew = init_gen(popsize, nsbit);
    fitness = zeros(1, popsize);

    x = range(1):0.1:range(2);
    plot(x, f(x));

    for i=1:popsize
        solNew(i) = bintodec(popNew(i, :));
    end

    for i=1:maxGen
        fitOld = fitness;
        pop = popNew;
        sol = solNew;
    

        for j=1:popsize
            ii = floor(popsize*rand) + 1;
            jj = floor(popsize*rand) + 1;
            
            if pc > rand
                [popNew(ii, :), popNew(jj, :)] = crossover(pop(ii, :), pop(jj, :));
    
                count = count + 2;
                evolve(ii);
                evolve(jj);
            end
    
            if pm > rand
                kk = floor(popsize*rand) + 1;
                count = count + 1;
                popNew(kk, :) = mutate(pop(kk, :), nsite);
                evolve(kk);
            end
        end

        bestFun(i) = max(fitness);
        bestSol(i) = mean(sol(bestFun(i) == fitness));
    end

    subplot(2, 1, 1); plot(bestSol); title('Best estimates');
    subplot(2, 1, 2); plot(bestFun); title('Fitness');
end

function pop = init_gen(np, nsbit)
    pop = rand(np, nsbit + 1) > 0.5;
end
    
function evolve(j)
    global solNew popNew fitness fitOld pop sol f;

    solNew(j) = bintodec(popNew(j, :));
    fitness(j) = f(solNew(j));

    if fitness(j) > fitOld(j)
        pop(j, :) = popNew(j, :);
        sol(j) = solNew(j);
    end
end

function [dec] = bintodec(bin)
    global range;

    nn = length(bin) - 1;
    num = bin(2:end);

    sign = 1 - 2*bin(1);
    dec = 0;
    dp = floor(log2(max(abs(range))));

    for i=1:nn
        dec = dec + num(i)*2^(dp - i);
    end

    dec = dec*sign;
end

function [c, d] = crossover(a, b)
    nn = length(a) - 1;

    cpoint = floor(nn*rand) + 1;
    c = [a(1:cpoint) b(cpoint + 1:end)];
    d = [b(1:cpoint) a(cpoint + 1:end)];
end

function anew = mutate(a, nsite)
    nn = length(a);
    anew = a;

    for i=1:nsite
        j = floor(rand*nn) + 1;
        anew(j) = mod(a(j) + 1, 2);
    end
end