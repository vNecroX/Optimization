clear all; clc;
syms x y
F = 10*x^2 + 5*x*y + 10*(y-3)^2;
lambda = 0.01;
VF = inline([diff(F,x);diff(F,y)])
x = [8;-8];

for i=1:100
    x = x - lambda * VF(x(1),x(2));
    A(i,1) = x(1);
    A(i,2) = x(2);
    A(i,3) = 10*x(1)^2 + 5*x(1)*x(2) + 10*(x(2)-3)^2;

    x = linspace(-8,8,30);
    y = x';
    z = 10.*x.^2 + 5.*x.*y + 10.*(y-3).^2;
    figure(1);
    surf(x,y,z);
    hold on;
    plot3(x(1),x(2),A(i,3),'o','Color','r','MarkerFaceColor','red','MarkerSize',5)
    pause(0.85)
end