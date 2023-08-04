% Gradient descent example

clear all

funstr = '10 - (exp(-1*(x^2 + 3*y^2)))';
f = vectorize(inline(funstr));
range = [-1 1 -1 1];

% Draw the function
Ndiv = 50;
dx = (range(2) - range(1))/Ndiv;
dy = (range(4) - range(3))/Ndiv;
[x, y] = meshgrid(range(1):dx:range(2), range(3):dy:range(4));
z = (f(x, y));
figure(1);
surfc(x, y, z);

% Define the number of iterations
k = 0;
niter = 200;

% Gradient step size h definition
hstep = 0.001;

% Step size of the Gradient descent method
alfa = 0.05;

% Initial point selection
xrange = range(2) - range(1);
yrange = range(4) - range(3);
x1 = rand*xrange + range(1);
x2 = rand*yrange + range(3);

% Optimization process
while(k < niter)
    % Function evaluation
    zn = f(x1, x2);

    % Computation of gradients gx1 & gx2
    vx1 = x1 + hstep;
    vx2 = x2 + hstep;
    gx1 = (f(vx1, x2) - zn)/hstep;
    gx2 = (f(x1, vx2) - zn)/hstep;

    % Draw the current position
    figure(2)
    contour(x, y, z, 15);
    hold on;
    plot(x1,x2,'.','MarkerSize','10','MarkerFaceColor','g');
    hold on;

    % Computation of the new solution
    x1 = x1 - alfa*gx1;
    x2 = x2 - alfa*gx2;
    k = k+1;
end