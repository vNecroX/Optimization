clear all
close all

D = 5000;
C = 5;
S = 49;
P = 0.2;
M = P*C;

funObj = @(xi) D*C + D/xi*S + xi/2*M;

N = 5;
d = 1;
lb = 400;
ub = 1200;

k = 0;
kmax = 150;
c1 = 2; 
c2 = 2;

for i=1:N
    x(i, :) = rand(1, d).*(ub - lb) + lb;
    v(i, :) = zeros(1, d);
end

for i=1:N
    xi = x(i, :);
    fx(i, :) = funObj(xi);
end

[gfit, ind] = min(fx);
g = x(ind, :);
fp = fx;
p = x;
axisx = lb:ub;
axisy = [];

for i=1:length(axisx)
    axisy(i) = funObj(axisx(i));
end

while k<kmax
    k = k + 1;

    figure(1);

    plot(x, fx, 'o', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
    pause(0.3);
    plot(x, fx, 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    pause(0.3);
    hold off;

    for i=1:N
        xi = x(i, :);
        pi = p(i, :);
        v(i, :) = v(i, :) + c1*rand(1, d).*(pi - xi) + c2*rand(1, d).*(g - xi);
    end

    x = x + v;
    
    for i=1:N
        for j=1:d
            if x(i, j) < lb(j)
                x(i, j) = lb(j);
            elseif x(i, j) > ub(j)
                x(i, j) = ub(j);
            end
        end
    end

    for i=1:N
        xi = x(i, :);
        fx(i, :) = funObj(xi);
    end

    [gfitkplus1, ind] = min(fx);

    if gfitkplus1 < gfit
        gfit = gfitkplus1;
        g = x(ind, :);
    end

    for i=1:N
        if fx(i, :) < fp(i, :)
            fp(i, :) = fx(i, :);
            p(i, :) = x(i, :);
        end
    end

    Evolution(k) = gfit;
end

figure;
plot(Evolution);
disp(['Best solution: ', num2str(g)]);
disp(['Best fitness: ', num2str(gfit)]);