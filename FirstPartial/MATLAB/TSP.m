clear all; close all; clc;

distancias = [0 120 220 150 21;
              120 0 100 110 130;
              220 100 0 160 185;
              150 110 160 0 190;
              21 130 185 190 0]; % Aprox 56

disp('Problema 1: agente viajero');

res = input('Ingrese la ruta a seguir entre corchetes: \n'); % [1 2 3 4 5 1]
inicial = res;
ruta = [1 5 2 3 4 1];
dOptima = 0;
dPropuesta = 0;

for i=1:5
    dOptima = dOptima + distancias(ruta(i), ruta(i + 1)); % Dist ruta
    dPropuesta = dPropuesta + distancias(res(i), res(i + 1)); % Ruta propuesta
end

% Ingredientes 
alfa = 0.55;
T = 100;
Tmin = 50; % Puede cambiar / 56/ 57
k = 1;
Lk = 3; % Puede cambiar, numero de repeticiones
Sp = res; % Auxiliar de solucion propuesta

while T>Tmin
    for L=1:Lk
        posicionVec = sort(randperm(5, 2));

        while posicionVec(1)==1 || posicionVec(2)==1 
            % Para que no se intercambien los elementos de inicio a fin
            posicionVec = sort(randperm(5, 2));
            cambio1 = posicionVec(1);
            cambio2 = posicionVec(2);
        end

        if k>=2
            while posicionVec(1)==cambio1 || posicionVec(2)==cambio2
                posicionVec = sort(randperm(5, 2));
            end
        end

        Sp(posicionVec(1):posicionVec(2)) = flip(res(posicionVec(1):posicionVec(2)));
        dFinal = 0;

        for i=1:5
            dFinal = dFinal + distancias(Sp(i), Sp(i + 1));
            disp(Sp(i));
            disp('\n');
            disp(Sp(i + 1));
        end

        % Evolucion diferencial
        deltaE = dFinal-dPropuesta;

        if deltaE<=0
            resp = Sp;
            dPropuesta = dFinal;
        elseif rand(1, 1)<exp(deltaE/T)
            resp = Sp;
            dPropuesta = dFinal;
        end
    end

    Ti = T;
    k = k+1;
    T = alfa*Ti;
end

% Se obtiene el fitness
Eficiencia = 1-((dPropuesta-dOptima)/dOptima);
fprintf('\n Ruta ingresada por el usuario \n');
disp(inicial);
fprintf('Ruta de acuerdo con el algoritmo AS');
disp(res);
fprintf('Distancia de la ruta: %.3f \n', dPropuesta);
fprintf('El fitness: %.3f \n', Eficiencia);