%---------------------------------------------------------
% Programa para generar una población
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear

% Crear una matriz de individuos con los siguientes datos
% |lower bound|upper bound|press
% En este caso son tres individuos
individuos = [-5, 5, 4;
              10,30, 2];

% Semilla de numeros aleatorios
seed = 12345; % Semilla

generator(individuos, 5, seed);

% Generador de numeros de 32 bits
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

function real = binary_real(binary, total, lb, ub)
    real = lb + ((binary/total)*(ub-lb));
end

function generator(individuos, num_individuos, seed)
    % Calcular el numero de bits por variable
    num_filas = size(individuos, 1);

    nb_total = zeros(1, num_filas);
    %display(nb_total);

    for i=1:num_filas
        % Determinamos los límites
        lb = individuos(i, 1);
        ub = individuos(i, 2);
        press = individuos(i, 3);

        % Calculamos el numero de bits
        nb = log((ub-lb)*(10^press)) + 0.9;
        %display(nb);
        nb_total(1, i) = ceil(nb);
    end
    %display(nb_total);

    sum = 0;
    for i=1:size(nb_total, 2)
        sum = sum + nb_total(1, i);
    end
    %display(sum);
    
    %Generar matriz
    b_matriz = zeros(num_individuos, sum);
    %display(b_matriz);

    % Rellenar con numeros aleatorios
    for i=1:num_individuos
        a = 0; b = 1;
        random_numbers = a + (b - a) * xorshift(seed, sum);
        random_numbers = round(random_numbers);
        b_matriz(i, :) = random_numbers;
        seed = seed + 1;
    end

    %display(b_matriz);
    %display(b_matriz(1,1:8));

    % Separar variables segun el numero de bits y guardar en una matriz
    matriz_final = zeros(num_individuos,size(nb_total, 2));
    for i=1:num_individuos
        first = 1;
        last = nb_total(1);
        for j=1:size(nb_total, 2)
           % Convertir a decimal, convertir primero a string
           binary = num2str(b_matriz(i,first:last));
           %display(binary);
           decimal = bin2dec(binary);
           %display(decimal);

           if j~=size(nb_total,2)
               first = last+1;
               last = last + nb_total(j+1);
           end

         real_number = binary_real(decimal, (2^nb_total(j)), individuos(j, 1), individuos(j, 2));
         %display(real_number);
         matriz_final(i,j) = real_number;
        end 
    end

    fprintf('Población con %d individuos y %d variables\n', num_individuos, size(nb_total, 2));
    disp(matriz_final);
end
% Determinar numero de individuos
% Determinar numero de variables que se quieren
% Determinar numero de bits
% Determinar presición
% Llenar de 0 y 1 random
% Mostrar matriz de 1's y 0's
% Convertir a real
% Mostrar matriz con numeros reales