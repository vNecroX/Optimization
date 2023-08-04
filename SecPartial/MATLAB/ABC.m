% ABC Example
function ABC(C)
    % Function Eggholder
    func = '((-1)*(y + 47)*(sin(sqrt(abs(y*(x/2) + 47))))) - (x*sin(sqrt(abs(x - (y + 47)))))';
    f = vectorize(inline(func));
    range = [-512 512 -512 512];
    % Initial parameter
    d = 2; Np = 100;
    food_source = round(Np/2);
    gmax = 150; % Maximum number generations
    Limit = 15; % Abandondment criterion
    Range = range(2)-range(1);
    Pop = (rand(food_source, d)*Range) + range(1);
    Ndiv = 100;
    dx = Range/Ndiv;
    dy = dx;
    [x, y] = meshgrid(range(1):dx:range(2), range(3):dy:range(4));
    z = f(x, y);
    figure, surfc(x, y, z);
    figure, contour(x, y, z, 15);
    
    for ii=1:food_source
        valfit(ii) = f(Pop(ii, 1), Pop(ii, 2));
        fitness(ii) = calculateFitness(valfit(ii));
    end
    
    figure, contour(x, y, z, 15); hold on;
    plot(Pop(:, 1), Pop(:, 2), 'b.', 'MarkerSize', 15);
    drawnow; hold on;
    test = zeros(1, food_source);
    BestInd = find(valfit == min(valfit));
    BestInd = BestInd(end);
    GlobalMin = valfit(BestInd);
    GlobalParams = Pop(BestInd, :);
    g = 0;
    
    while(g < gmax)
        for i=1:food_source
            Param2Change = fix(rand*d) + 1;
            neighbor = fix(rand*(food_source)) + 1;
    
            while(neighbor == 1)
                neighbor = fix(rand*(food_source)) + 1;
            end
    
            solutions = Pop(i, :);
            
            solutions(Param2Change) = Pop(i, Param2Change) + (Pop(i, Param2Change) - ...
                                      Pop(neighbor, Param2Change))*(rand - 0.5)*2;
            ind = find(solutions < range(1));
            solutions(ind) = range(1);
            ind = find(solutions > range(2));
            solutions(ind) = range(2);
            valfitSol = f(solutions(1), solutions(2));
            fitnessSol = calculateFitness(valfitSol);
    
            if(fitnessSol > fitness(1))
                Pop(i, :) = solutions;
                fitness(i) = fitnessSol;
                valfit(i) = valfitSol;
                test(i) = 0;
            else
                test(i) = test(i) + 1;
            end
        end
    
        probab = (0.9.*fitness./max(fitness)) + 0.1;
        i = 1;
        T = 0;
    
        while(T < food_source)
            if(rand < probab(i))
                T = T + 1;
                Param2Change = fix(rand*d) + 1;
                neighbor = fix(rand*(food_source)) + 1;
        
                while(neighbor == 1)
                    neighbor = fix(rand*(food_source)) + 1;
                end
    
                solutions = Pop(i, :);
    
                solutions(Param2Change) = Pop(i, Param2Change) + (Pop(i, Param2Change) - ...
                                      Pop(neighbor, Param2Change))*(rand - 0.5)*2;
                ind = find(solutions < range(1));
                solutions(ind) = range(1);
                ind = find(solutions > range(2));
                solutions(ind) = range(2);
                valfitSol = f(solutions(1), solutions(2));
                fitnessSol = calculateFitness(valfitSol);
    
                if(fitnessSol > fitness(1))
                    Pop(i, :) = solutions;
                    fitness(i) = fitnessSol;
                    valfit(i) = valfitSol;
                    test(i) = 0;
                else
                    test(i) = test(i) + 1;
                end
            end
            
            i = i + 1;
    
            if(i == (food_source) + 1)
                i = 1;
            end
        end
    
        % The best food source is stored
        ind = find(valfit == min(valfit));
        ind = ind(end);
        
        if(valfit(ind) < GlobalMin)
            GlobalMin = valfit(ind);
            GlobalParams = Pop(ind, :);
        end
        
        ind = find(test == max(test));
        ind = ind(end);
        
        if(test(ind) > Limit)
            test(ind) = 0;
            solutions = (Range).*rand(1, d) + range(1);
            valfitSol = f(solutions(1), solutions(2));
            fitnessSol = calculateFitness(valfitSol);
            Pop(ind, :) = solutions;
            fitness(ind) = fitnessSol;
            valfit(ind) = valfitSol;
        end
        
        g = g + 1;
        clc;
        
        disp(GlobalMin);
        disp(GlobalParams);
    end
end

function fFitness = calculateFitness(fobjv)
    fFitness = zeros(size(fobjv));
    ind = find(fobjv >= 0);
    fFitness(ind) = 1./(fobjv(ind) + 1);
    ind = find(fobjv < 0);
    fFitness(ind) = 1 + abs(fobjv(ind));
end