 %---------------------------------------------------------
% Programa de PSO
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear
rng(0,'twister');
% Crear una matriz de individuos con los siguientes datos
% | lower bound  |   upper bound   |    press (filas)
% 3 variables que cumplen -5 <= x <= 5
individuos = [-5,5,3;
             -5,5,3;
             -5,5,3];

%%%%%%%%%%%%%%%%%% FUNCIONES PARA POBLACIÓN ALEATORIA %%%%%%%%%%%%%%%%%%%%
% Semilla para generar numeros aleatorios
seed = 3456; % Semilla

% Llamar a la función del generador
for i=1:30
    pso(individuos, 50, seed, 0.3, 1.7, 1.5, 100);
end

% Funcion objetivo
function f = f_x1x2x3(x1,x2,x3)
    f = (x1^2)+(x2^2)+(x3^2);
end

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

%%%%%%%%%%%%%%%%% GENERADOR DE PARTÍCULAS ALEATORIAS %%%%%%%%%%%%%%%%%%%%
function [par, nb_total] = particles_genertor(individuos, num_individuos, seed)
    % Calcular el numero de bits por variable
    num_variables = size(individuos, 1); % Numero de variables
    nb_total = zeros(1, num_variables); % Matriz para guardar el numero de bits por variable

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % POBLACIÓN INICIAL ALEATORIA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:num_variables
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
        random_numbers = rand(1,sum); % Llamar al generador de aleatorios
        random_numbers = round(random_numbers); % Redondear el numero
        b_matriz(i, :) = random_numbers; % Guardar en la matriz
        seed = seed + 1; % Cambiar la semilla
    end

    for i=1:size(nb_total,2)
        %fprintf('Numero de bits de la variable %d: %d\n', i,nb_total(1,i));
    end

    %fprintf('\nMatriz binaria\n');
    %disp(b_matriz);

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

    %fprintf('%d  %d\n', size(matriz_final,1), size(matriz_final,2));

    par = matriz_final;
    %fprintf('%d  %d\n', size(par,1), size(par,2));
end

function pso(individuos, num_individuos, seed, w, c1, c2, generaciones)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERADOR DE PARTICULAS ALEATORIA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [matriz_final, nb_total] = particles_genertor(individuos, num_individuos, seed);

    %matriz_final = results(1);
    %nb_total = results(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUACIÓN DE LA FUNCIÓN OBJETIVO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_variables = size(individuos, 1);
    matriz_nueva = zeros(num_individuos, size(matriz_final,2)+1);

    %fprintf('Matriz final: %d %d\n', size(matriz_final,1), size(matriz_final,2))
    %fprintf('Matriz final: %d %d\n', size(matriz_nueva,1), size(matriz_nueva,2))

    for i=1:num_individuos
        for j=1:num_variables
            matriz_nueva(i,j) = matriz_final(i,j);
        end
        matriz_nueva(i,j+1) = f_x1x2x3(matriz_final(i,1), matriz_final(i,2), matriz_final(i,3));
    end

    % Imprimir la matriz
    %fprintf('Población con %d partículas y %d variables\n', num_individuos, size(nb_total, 2));
    %disp(matriz_final);
    %fprintf('Evaluación de partículas en la función objetivo\n');
    %disp(matriz_nueva);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % G BEST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    total_gbest = zeros(num_individuos, num_variables+1);
    g_best_val = min(matriz_nueva(:,4));
    %fprintf('G Best: '); disp(g_best_val);
    g_best_i = find(matriz_nueva(:,4)==g_best_val);
    %fprintf('G Best: '); disp(g_best_i);
    g_best = matriz_nueva(g_best_i,:);
    %fprintf('G Best: '); disp(g_best);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VELOCIDADES INICIALES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocidades inicializadas en 0
    velocidades = zeros(num_individuos, num_variables);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P BEST INICIAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_best = matriz_nueva;
    %fprintf('P best:\n');
    %disp(p_best);

    for gen=1:generaciones
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ACTUALIZAR VELOCIDADES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:num_individuos
            r1 = rand();
            %fprintf('Randi: %d  ', r1);
            r2 = rand();
            %fprintf('Randi: %d\n', r2);
            velocidades(i,:) = (w * velocidades(i,:)) + (c1 * r1 * (p_best(i,1:end-1) - matriz_nueva(i,1:end-1))) + (c2 * r2 * (g_best(1:end-1) - matriz_nueva(i,1:end-1)));
        end
        %fprintf('Nuevas Velocidades.\n');
        %disp(velocidades);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ACTUALIZAR POSICIÓN DE PARTÍCULAS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:num_individuos
            matriz_nueva(i,1:end-1) = matriz_nueva(i,1:end-1) + velocidades(i,:);
            matriz_nueva(i,end) = f_x1x2x3(matriz_nueva(i,1), matriz_nueva(i,2), matriz_nueva(i,3));
        end
        %fprintf('Nuevas partículas.\n');
        %disp(matriz_nueva);
    
        for i=1:num_individuos
            %fprintf('Matriz\n%d\n', p_best(i,end));
            if matriz_nueva(i,end) < p_best(i,end)
                p_best(i,:) = matriz_nueva(i,:);
                %fprintf('Es menor\n%d\n',matriz_nueva(i,:));
            else
                %fprintf('No es menor\n'); disp(matriz_nueva(i,:));
            end
        end
        %fprintf('Nuevo p best.\n');
        %disp(p_best);
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ACTUALIZAR G BEST
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       g_best_val = min(matriz_nueva(:,4));
       %fprintf('G Best: '); disp(g_best_val);
       g_best_i = find(matriz_nueva(:,4)==g_best_val);
       %fprintf('G Best: '); disp(g_best_i);
       if matriz_nueva(g_best_i,:) < g_best(4)
           g_best = matriz_nueva(g_best_i,:);
       end

       total_gbest(gen,:) = g_best;
       if g_best(4) <= 0.00001
           fprintf('Llego en la iteracion: %d\n', gen);
           break;
       end
    end
    fprintf('Global Best');
    disp(g_best);
    
    %fprintf('Mínimos en cada iteración: \n'); disp(sortrows(total_gbest,4));
    fprintf('Costo computacional: %d\n', num_individuos * generaciones);

end
