%---------------------------------------------------------
% Programa para generar numeros aleatorios
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

% Para generar un solo numero aleatorio entre 0 y 1
random_num = rand;

% Para generar un array de numeros aleatorios
N = 5 % Filas
random_array = rand(1, N);

% Para generar una matriz de números aleatorios
M = 2 % Columnas
random_matrix = rand(M, N);

% Para cambiar limites de numeros aleatorios
a = 7 %lb
b = 15 %ub

% Numero
random_num2 = a + (b - a) * rand;
% Array
random_array2 = a + (b - a) * rand(1, N);
% Matriz
random_matrix2 = a + (b - a) * rand(M, N);