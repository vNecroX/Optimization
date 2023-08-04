%"PSO 2D" Algorithm PSO Particle Swarm Optimization
%Quintero Sotomayor Esteban

function [best] = B3_pso(n, num_iterations)
    clear all; clc; 
    close all;
    if nargin<2, num_iterations=20; end
    if nargin<1, n=40; end
    fstr = '-sin(x)*(sin(x^2/pi))^20 - sin(y)*(sin(2*y^2/pi))^20';
    f = vectorize(inline(fstr));
    range = [1 3 1 3];
    alpha = 0.2; beta = 0.5;
    Ndiv = 100;
    dx = (range(2)-range(1))/Ndiv; 
    dy = (range(4)-range(3))/Ndiv;
    xgrid = range(1):dx:range(2); 
    ygrid = range(3):dy:range(4);
    [x, y] = meshgrid(xgrid, ygrid);
    z = f(x, y);
    figure(1); surfc(x, y, z);
    best = zeros(num_iterations, 3);
    [xn, yn] = init_pso(n, range);
    figure(2);
    
    for i=1:num_iterations
        contour(x, y, z, 15); hold on;
        zn = f(xn, yn);
        zn_min = min(zn);
        x0 = min(xn(zn==zn_min));
        y0 = min(yn(zn==zn_min));
        z0 = min(zn(zn==zn_min));
        plot(xn, yn, '.', x0, y0, '*'); axis(range);
        [xn, yn] = pso_move(xn, yn, x0, y0, alpha, beta, range);
        drawnow;
        hold off;
        best(i, 1) = x0; best(i, 2) = y0; best(i, 3) = z0;
    end

function [xn, yn] = init_pso(n, range)
    xrange = range(2) - range(1); yrange = range(4) - range(3);
    xn = rand(1, n)*xrange + range(1);
    yn = rand(1, n)*yrange + range(3);

function [xn, yn] = pso_move(xn, yn, x0, y0, a, b, range)
    nn = size(yn, 2);
    xn = xn.*(1 - b) + x0.*b + a.*(rand(1, nn) - 0.5);
    yn = yn.*(1 - b) + y0.*b + a.*(rand(1, nn) - 0.5);
    [xn, yn] = findrange(xn, yn, range);

function [xn, yn] = findrange(xn, yn, range)
    nn = length(yn);

    for i=1:nn
        if xn(i) <= range(1), xn(i) = range(1); end
        if xn(i) >= range(2), xn(i) = range(2); end
        if yn(i) <= range(3), yn(i) = range(3); end
        if yn(i) >= range(4), yn(i) = range(4); end
    end