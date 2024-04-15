%---------------------------------------------------------
% Programa de optimizacion utilizando fmincon
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

% Vector para inicializar los puntos iniciales de las variables x1, x2, x3, x4 y x5
x0 = ones(1,5);

% Limites de los valores minimos que pueden tomar (low bound)
lb = ones(1,5);

% Matriz de coeficientes de las funciones de restriccion
% A = [-8, -2, -15, -4, -30;
%      -9, -3, -3, -1, -9;
%      -35, -3, -17, -1, -16;
%      -100, -90, -350, -200, -410;
%      -10, -20, -40, -25, -40;
%      -1, 0, 0, 0, 0;
%      0, -1, 0, 0, 0;
%      0, 0, -1, 0, 0;
%      0, 0, 0, -1, 0;
%      0, 0, 0, 0, -1];

A = [-8, -2, -15, -4, -30;
     -9, -3, -3, -1, -9;
     -35, -3, -17, -1, -16;
     -100, -90, -350, -200, -410;
     -10, -20, -40, -25, -40];

% Desigualdades
%b = [-60; -100; -120; -2100; -400; -1; -1; -1; -1; -1];
b = [-60; -100; -120; -2100; -400];

% Llamada a la funcion fmincon
% Recibe como parametros la funcion objetivo, los puntos iniciales,
% la matriz de restricciones y sus coeficientes, y los límites mínimos
% que pueden tomar cada variable, en este caso 0
[x,fval] = fmincon(@f, x0, A, b, [], [], lb, []);

fprintf('Porciones por nutrientes\n');
fprintf('Leche (x1)= %f\n', x(1));
fprintf('Jugo (x2) = %f\n', x(2));
fprintf('Pescado (x3) = %f\n', x(3));
fprintf('Papas fritas (x4) = %f\n', x(4));
fprintf('Pollo (x5) = %f\n', x(5));

fprintf('El mínimo valor de gasto es: %f\n', fval);

% Función objetivo
function val_f = f(x)
    val_f = (1.1*x(1)) + (1.2*x(2)) + (2*x(3)) + (1.3*x(4)) + (3*x(5));
end