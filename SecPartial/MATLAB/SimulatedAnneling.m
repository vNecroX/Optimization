% SIMULATED ANNELING
clear all;

% Definition of objective function
funstr = '3*(1 - x).^2.*exp(-(x.^2) - (y + 1).^2) - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2 - y.^2) - 1/3*exp(-(x + 1).^2 -y.^2)';
f = vectorize(inline(funstr));
range = [-3 3 -3 3];

% Draw the objetive function 
Ndiv = 50;
dx = (range(2) - range(1))/Ndiv;
dy = (range(4) - range(3))/Ndiv;
[x, y] = meshgrid(range(1):dx:range(2), range(3):dy:range(4));
z = f(x, y);

figure(1), surfc(x, y, z);

k = 1; valid = 0; niter = 150; temp = [];
T_ini = 10;
T_fin = 1e-10;
beta = 0.95;

% Initialization of the candidate solution
xrange = range(2) - range(1);
yrange = range(4) - range(3);
xn = rand*xrange + range(1);
yn = rand*yrange + range(3);

T = T_ini;

while(k < niter)
    figure(2);
    contour(x, y, z); hold on;
    plot(xn, yn, '.', 'MarkerSize', 10, 'MarkerFaceColor','g');
    drawnow; hold off;
    E_old = f(xn, yn);
    % disp(E_old);

    while(valid == 0)
        xnc = xn*rand*2.5*T;
        ync = yn*rand*2.5*T;

        if((xnc >= range(1)) && (xnc <= range(2)) && (ync >= range(3)) && (ync <= range(4)))
            valid = 0;
        end
    end

    valid = 0;
    data(k);
    E_new = f(xnc, ync);

    DeltaE = E_new - E_old;

    if(DeltaE < 0)
        xn = xnc;
        yn = ync;
    end

    if((DeltaE < 0) && (exp(DeltaE/T) > rand))
        xn = xnc;
        yn = ync;
    end

    T = beta*T;

    if(T < T_fin)
        T = T_fin;
    end

    k = k + 1;
end