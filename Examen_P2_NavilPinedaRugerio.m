%---------------------------------------------------------
% Programa de Evolución Diferencial
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre

% El algoritmo que coloque en el examen comprobé que no me da los mejores
% resultados por lo que lo cambie a Evolucion Diferencial y utilizando las
% reglas de DEB comprobé aquellos valores de variables que cumplen con la
% restricción establecida en mi planteamiento de problema.
% 1. Se genera la poblacion inicial.
% 2. Se hace mutación y cruza.
% 3. Se comprueba con regla de DEB el mejor individuo ya sea el padre o el
% hijo.
% 4. Se imprimen los resultados-
%---------------------------------------------------------

clc
clear
rng(0,'twister');

% Crear una matriz de individuos con los siguientes datos
% | lower bound  |   upper bound   |    press (filas)
% Numero de variables (columnas)
individuos = [0,4,3;
              0,8,3];

% Semilla para generar numeros aleatorios
seed = 3456; % Semilla

% Llamar a la función del generador
% Llamar a la función del generador 30 veces
for i=1:30
rand_1_bin(individuos, 50, seed, 0.1, 0.5, 100, i);
end

% Comentario final
fprintf(['De este modo se comprueba que el minimo valor que puede tomar la función objetivo es 0, esto calculado mediante \n' ...
'Evolución Diferencias y que una buena sintonización es de 0.1 para factor de mutación, un factor de cruza de 0.5, 50 \n' ...
    'individuos, y 100 generaciones, lo que genera un costo computacional de 5000\n'])

% Funcion objetivo
function f = f_x1x2(x1,x2)
    %f= (3 * (1 - x1)^2 * exp(-x1^2 - (x2 + 1)^2)) + (10 * ((x1/5) - x1^3 - x2^5) * exp(-x1^2 - x2^2)) + ((1/3) * exp(-(x1 + 1)^2 - x2^2));
    f = 2*x1+3*(x2^2)-22;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo de Evolucion Diferencial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rand_1_bin(individuos, num_individuos, seed, mut_factor, cr_factor, num_generaciones, iter_ep)
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
        random_numbers = a + (b - a) * xorshift(seed, sum); % Llamar al generador de aleatorios
        random_numbers = round(random_numbers); % Redondear el numero
        b_matriz(i, :) = random_numbers; % Guardar en la matriz
        seed = seed + 1; % Cambiar la semilla
    end

    %for i=1:size(nb_total,2)
        %fprintf('Numero de bits de la variable %d: %d\n', i,nb_total(1,i));
    %end

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUACIÓN DE LA FUNCIÓN OBJETIVO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matriz_nueva = zeros(num_individuos,size(nb_total, 2)+1);
    for i=1:num_individuos
        for j=1:num_variables
            matriz_nueva(i,j) = matriz_final(i,j);
        end
        matriz_nueva(i,end) = f_x1x2(matriz_final(i,1), matriz_final(i,2));
    end
    % for i=1:num_individuos
    %     matriz_nueva(i,1) = matriz_final(i,1);
    %     matriz_nueva(i,2) = matriz_final(i,2);
    %     % Evaluar en la función objetivo
    %     matriz_nueva(i,3) = f_x1x2(matriz_final(i,1), matriz_final(i,2));
    % end

    % Imprimir la matriz
    %fprintf('Población con %d individuos y %d variables\n', num_individuos, size(nb_total, 2));
    %disp(matriz_final);
    %fprintf('Evaluacuón de la población en la función objetivo\n');
    %disp(matriz_nueva);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECORRER GENERACIONES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    may_child = -10000;
    min_child = 10000;
    val_x_min = 0;
    val_y_min = 0;

    for gen=1:num_generaciones
        % Matriz para los nuevos hijos
        hijos = zeros(num_individuos, num_variables+1);
    
        n = 1;
        while n <= num_individuos
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INDIVIDUOS AUXILIARES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            actual = matriz_nueva(n,1:end-1);
    
            indv = zeros(1,3);
    
            for i=1:3
                rand_indv = randi(num_individuos);
                % Verifica si el número generado ya existe en la lista
                while any(rand_indv == indv) || rand_indv == n
                    rand_indv = randi(num_individuos);
                end
                indv(i) = rand_indv;
            end
    
            % Individuos seleccionados
            %disp(indv);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DIFERENCIA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Se eligen los dos padres al azar y se calcula su diferencia
            dif_indv = zeros(1,3);
    
            for i=1:3
                rand_indv = randi(3);
                % Verifica si el número generado ya existe en la lista
                while any(rand_indv == dif_indv)
                    rand_indv = randi(3);
                end
                dif_indv(i) = rand_indv;
            end
    
            parent1 = indv(dif_indv(1));
            parent2 = indv(dif_indv(2));
            suma = indv(dif_indv(3));
    
            %fprintf('Padre 1: %d\n', parent1);
            %fprintf('Padre 2: %d\n', parent2);
    
            % Diferencia entre variables
            new_dif = matriz_nueva(parent1,1:end-1) - matriz_nueva(parent2,1:end-1);
            %fprintf('Diferencia: ');
            %disp(new_dif);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MUTACIÓN
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mutación (multiplicar por un factor de mutación)
            new_mut = new_dif * mut_factor;
            %fprintf('Mutacion: ');
            %disp(new_mut);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SUMA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fprintf('Suma de: ');
            %disp(matriz_nueva(suma,1:end-1));
            %fprintf('con: ');
            %disp(new_mut);
            new_mut = new_mut + matriz_nueva(suma,1:end-1);
            %fprintf('Resultado: ');
            %disp(new_mut);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % VERIFICACIÓN DE COTAS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Verificar si está dentro de las cotas
            % Manejo de cotas
            indx1 = [1];
            indx2 = [1];
            
            while ~isempty(indx1) || ~isempty(indx2)
                for l=1:num_variables
                    if new_mut(l) < individuos(l,1)
                        new_mut(l) = 2 * individuos(l,1) - new_mut(l);
                    end
                    if new_mut(l) > individuos(l,2)
                        new_mut(l) = 2 * individuos(l,2) - new_mut(l);
                    end
                    indx1 = find(new_mut < individuos(l,1));
                    indx2 = find(new_mut > individuos(l,2));
                end    
            end
            %fprintf('\nResultado dentro de las cotas: ');
            %disp(new_mut);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CRUZA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Aplicar cruza
            cross_child = zeros(1,num_variables);
            jrand = randi(num_variables);
            for m=1:num_variables
                r = rand;
                r = round(r, 1);
                %fprintf('Rand: %d', rand);
                if r < cr_factor || m == jrand
                    cross_child(1,m) = new_mut(m);
                else
                    cross_child(1,m) = actual(m);
                end
            end
            
            %fprintf('\nCombinación por cruza:');
            %fprintf('\nPadre:');
            %disp(actual);
            %fprintf('Ruido:');
            %disp(new_mut);
            %fprintf('Hijo:');
            %disp(cross_child);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ELEGIR PADRE VS HIJO (REGLAS DE DEB)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Checar si son factibles
            % Verifica las restricciones para el padre y para el hijo
            % Restriccion: que cumpla con 2x+y-8=0
            ap_ambos = 0;
            ap_padre = 0;
            ap_hijo = 0;
            % Cero gordo
            cero = 0.00001;
            if abs(2 * cross_child(1)+cross_child(2)-8.0) < cero
                % El hijo es apto
                ap_hijo = 1;
            end
            if abs(2 * actual(1)+actual(2)-8.0) < cero
                % El padre es apto
                ap_padre = 1;
            end

            if ap_hijo == 1 && ap_padre == 1
                % Ambos son aptos
                ap_ambos = 1;
            end

            func_padre = f_x1x2(actual(1),actual(2));
            func_hijo = f_x1x2(cross_child(1), cross_child(2));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ASIGNACIÓN DE HIJOS A LA MATRIZ
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ap_ambos == 1
                %fprintf('Ambos aptos\n');
                if func_padre < func_hijo
                    %fprintf('Pasa el padre\n');
                    for var=1:num_variables
                        hijos(n, var) = actual(var);
                    end
                elseif func_hijo < func_padre
                    %fprintf('Pasa el hijo\n');
                    for var=1:num_variables
                        hijos(n, var) = cross_child(var);
                    end
                end
                hijos(n, end) = f_x1x2(hijos(n,1), hijos(n,2));
                %disp(hijos(n,:));
            elseif ap_padre == 1
                %fprintf('Padre apto\n');
                for var=1:num_variables
                    hijos(n, var) = actual(var);
                end
            elseif ap_hijo == 1
                %fprintf('Hijo apto\n');
                for var=1:num_variables
                    hijos(n, var) = cross_child(var);
                end
            else
                %fprintf('Ninguno apto\n');
                for var=1:num_variables
                    hijos(n, var) = actual(var);
                end
                % for var=1:num_variables
                %     hijos(n, var) = cross_child(var);
                % end
                % hijos(n,:) = []
                % num_individuos = num_individuos - 1
            end
    
            n = n + 1;
        end

        % Hijos
        %fprintf('NUEVA GENERACIÓN:\n');
        %disp(hijos);

        % El mayor hijo y su iteración
        if min_child > min(hijos(:,end))
            min_child = min(hijos(:,end));
            val_var_i = find(hijos(:,end)==min_child);
            val_x_min = hijos(val_var_i,1);
            val_y_min = hijos(val_var_i,2);
        end

        if may_child < max(hijos(:,end))
            may_child = max(hijos(:,end));
        end

        matriz_nueva = hijos;
    end

    fprintf('Costo computacional de la iteración %d: %d\n', iter_ep ,num_individuos * num_generaciones);
    fprintf('Minimo (mejor valor) de la iteración %d: %d\n', iter_ep, min_child);
    %fprintf('Valores que generan el mínimo: x=%d y=%d\n', iter_ep, min_child);
    fprintf('Maximo (peor valor) de la iteración %d: %d\n', iter_ep,may_child);
    fprintf('Promedio de la iteracion %d: %d\n', iter_ep,(mean(matriz_nueva(:,end))));
    fprintf('Desviacion estandar de la iteración %d: %d\n\n', iter_ep, std(matriz_nueva(:,end)));

end