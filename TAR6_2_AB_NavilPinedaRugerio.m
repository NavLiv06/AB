%---------------------------------------------------------
% Programa para generar una población aleatoria
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear

% Crear una matriz de individuos con los siguientes datos
% | lower bound  |   upper bound   |    press (filas)
% Numero de variables (columnas)
% individuos = [0,10,3;
%              -1,1,10;
%              -1, 1, 2];

individuos = [0,10,3;
             -100,100,10;
             -1, 1, 2];

% Semilla para generar numeros aleatorios
seed = 3456; % Semilla

% Llamar a la función del generador
generator(individuos, 5, seed);

% Generador de numeros de 32 bits basado en XORShift
function random_numbers = xorshift(seed, num)
    random_numbers = zeros(1, num); % Inicializamos un arreglo para guardar los numeros aleatorios
    state = uint32(seed); % El primer "estado" es la semilla, a partir de ahi iteramos
    for i = 1:num
        state = bitxor(bitshift(state, -13), state); % Desplazar 13 bits a la izquierda (shift) y hacer operacion xor bit a bit
        state = bitxor(bitshift(state, 17), state); % Desplazar 17 bits a la derecha (shift) y hacer operacion xor bit a bit
        state = bitxor(bitshift(state, -5), state); % Desplazar 5 bits a la izquierda (shift) y hacer operacion xor bit a bit
        random_numbers(i) = double(bitand(state, intmax('uint32'))) / double(intmax('uint32')); % Asegurar rango de 32 bits y convertir a rango de 0 a 1
    end
end

% Funcion para calcular numero de bits
function nb = nbits(ub, lb, press)
    nb = log((ub-lb)*(10^press)) + 0.9;
end

% Funcion para convertir de binario a decimal
function decimal = binary_decimal(number)
    string = num2str(number);
    decimal = bin2dec(string);
end

% Función para convertir de numero binario a real
function real = binary_real(binary, total, lb, ub)
    real = lb + ((binary/total)*(ub-lb));
end

% Generador de una población aleatoria
function generator(individuos, num_individuos, seed)
    % Calcular el numero de bits por variable
    num_filas = size(individuos, 1); % Numero de variables
    nb_total = zeros(1, num_filas); % Matriz para guardar el numero de bits por variable

    for i=1:num_filas
        % Determinamos los límites
        lb = individuos(i, 1);
        ub = individuos(i, 2);
        % Determinamos la presicion
        press = individuos(i, 3);

        % Calculamos el numero de bits
        nb = nbits(ub, lb, press);
        % Guardamos en la matriz el numero de bits ya redondeado
        nb_total(1, i) = ceil(nb);
    end

    % Sumamos el numero total de bits por todas las variables
    sum = 0;
    for i=1:size(nb_total, 2)
        sum = sum + nb_total(1, i);
    end
    
    %Generar matriz de binarios
    b_matriz = zeros(num_individuos, sum);
    % Rellenar con numeros aleatorios
    for i=1:num_individuos
        a = 0; b = 1; % Numeros de 0 a 1
        random_numbers = a + (b - a) * xorshift(seed, sum); % Llamar al generador de aleatorios
        random_numbers = round(random_numbers); % Redondear el numero
        b_matriz(i, :) = random_numbers; % Guardar en la matriz
        seed = seed + 1; % Cambiar la semilla
    end

    for i=1:size(nb_total,2)
        fprintf('Numero de bits de la variable %d: %d\n', i,nb_total(1,i));
    end

    fprintf('\nMatriz binaria\n');
    disp(b_matriz);

    % Separar variables segun el numero de bits y guardar en una matriz
    matriz_final = zeros(num_individuos,size(nb_total, 2));
    for i=1:num_individuos
        first = 1;
        last = nb_total(1);
        for j=1:size(nb_total, 2)
           % Convertir a decimal, convertir primero a string
           decimal = binary_decimal(b_matriz(i,first:last));

           if j~=size(nb_total,2)
               first = last+1;
               last = last + nb_total(j+1);
           end

           % Convertir a numero real
           real_number = binary_real(decimal, (2^nb_total(j)), individuos(j, 1), individuos(j, 2));
           % Agregar a la matriz final
           matriz_final(i,j) = real_number;
        end 
    end

    matriz_nueva = zeros(num_individuos,size(nb_total, 2)+1);
    for i=1:num_individuos
        matriz_nueva(i,1) = matriz_final(i,1);
        matriz_nueva(i,2) = matriz_final(i,2);
        matriz_nueva(i,3) = matriz_final(i,3);
        matriz_nueva(i,4) = ((matriz_final(i,1)^2) + (matriz_final(i,2)^2) + (matriz_final(i,3)^2));
    end
    

    % Imprimir la matriz
    fprintf('Población con %d individuos y %d variables\n', num_individuos, size(nb_total, 2));
    disp(matriz_final);
    disp(matriz_nueva);
end