%---------------------------------------------------------
% Programa de Harmony Search
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear
rng('shuffle');
%rng(0,'twister');

individuos = [0,1200,5;
              0,1200,5;
             -0.55,0.55,5;
             -0.55,0.55, 5];

ejecuciones = 10;
finals = zeros(ejecuciones,5);
fprintf('Algoritmo Harmony Search, MOD 2...\n');
r_accept = 0.8;
r_pa = 0.6;
bw = 0.5;
r_ia = 0.3;
iteraciones = 10000;
for i=1:ejecuciones
    finals(i, :) = harmony_generator(individuos, 100, r_accept, r_pa, bw, r_ia, iteraciones);
end
fprintf('Mejores resultados:\n');
fprintf('    x1           x2        x3        x4      f(x)\n');
finals = sortrows(finals, size(finals, 2));
disp(finals);
%harmony_generator(individuos, 10, 0.9, 0.5, 0.7, 0.1, 100000);

% Función objetivo
function f = f_x(x)
    f = 3 * x(1) + 0.000001 * x(1)^3 + 2 * x(2) + (0.000002/3) * x(2)^3;
end

% Para restricciones
function [c, ceq] = constraints(x)
    c = [-x(4) + x(3) - 0.55;
         -x(3) + x(4) - 0.55];
    
    ceq = [1000*(sin(-(x(3)) - 0.25)) + 1000*(sin(-(x(4)) - 0.25)) + 894.8 - (x(1));
           1000*(sin((x(3)) - 0.25)) + 1000*(sin((x(3)) - (x(4)) - 0.25)) + 894.8 - (x(2));
           1000*(sin((x(4)) - 0.25)) + 1000*(sin((x(4)) - (x(3)) - 0.25)) + 1294.8];
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

%%%%%%%%%%%%%%%%% GENERADOR DE NOTAS ALEATORIA %%%%%%%%%%%%%%%%%%%%
function [par, nb_total, bounds] = first_hm(individuos, num_individuos)
    % Calcular el numero de bits por variable
    num_variables = size(individuos, 1); % Numero de variables
    nb_total = zeros(1, num_variables); % Matriz para guardar el numero de bits por variable
    bounds = zeros(num_variables, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % POBLACIÓN INICIAL ALEATORIA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:num_variables
        % Determinamos los límites
        lb = individuos(i, 1);
        ub = individuos(i, 2);
        bounds(i, :) = [lb,ub];
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
        random_numbers = rand(1,sum);
        random_numbers = round(random_numbers);
        b_matriz(i, :) = random_numbers; % Guardar en la matriz
    end

    for i=1:size(nb_total,2)
        fprintf('Numero de bits de la variable %d: %d\n', i,nb_total(1,i));
    end

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

    par = matriz_final;
end



function best_harmony = harmony_generator(individuos, num_individuos, r_accept, r_pa, bw, r_ia, iteraciones)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERADOR DE PARTICULAS ALEATORIA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[matriz_final, nb_total, bounds] = first_hm(individuos, num_individuos);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVALUACIÓN DE LA FUNCIÓN OBJETIVO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_variables = size(individuos, 1);
    %matriz_nueva = zeros(num_individuos, size(matriz_final,2)+1);
    %for i=1:num_individuos
        %for j=1:num_variables
            %matriz_nueva(i,j) = matriz_final(i,j);
        %end
        %matriz_nueva(i,j+1) = f_x(matriz_final(i,:));
    %end
    HM = zeros(num_individuos, num_variables+1);
    for i = 1:num_variables
        lower_bound = individuos(i, 1);
        upper_bound = individuos(i, 2);
        HM(:, i) = lower_bound + (upper_bound - lower_bound) * rand(num_individuos, 1);
    end

    for i = 1:num_individuos
        HM(i, end) = f_x(HM(i, :));
    end
    %fprintf('Evaluación de partículas en la función objetivo\n');
    %disp(HM);

    newGen = zeros(1,num_variables+1);
    Xworst = zeros(1,num_variables+1);
    n = 1;
    cero_gordo = 0.00001;
    X = zeros(iteraciones,num_variables+2);

    for gen=1:iteraciones
        for var=1:num_variables
            %bw = (individuos(var,2) - individuos(var,2)) / num_individuos;
            rand_1 = round(rand * 100) / 100;
            if rand_1 < r_accept
                index = randi(num_individuos);
                rand_2 = round(rand * 100) / 100;
                if rand_2 < r_pa
                    ajuste = HM(index, var) + bw * (-1+2*rand);
                    %fprintf('Cambia a ajuste en el indice %d: ', index); disp(ajuste);
                    newGen(var) = ajuste;
                else
                    % Modificacion
                    rand_3 = round(rand * 100) / 100;
                    if rand_3 < r_ia
                         best_val = min(HM(:,end));
                         best_i = find(HM(:,end)==best_val,1);
                         %fprintf('Cambia a best: '); disp(HM(best_i,:));
                         newGen(var) = HM(best_i,var);
                    else
                         %fprintf('Cambia a normal: '); disp(HM(ind,index));
                         newGen(var) = HM(index, var);
                    end
                    %newGen(var) = HM(index, var);
                end
            else
                % newRandHarm = rand;
                % rango = 5 - (-5);
                % newRandHarm = newRandHarm * rango + (-5);
                % %fprintf('Cambia a radom: '); disp(newRandHarm);
                % newGen(var) = newRandHarm;
                switch var
                    case {1, 2}
                        Li = 0;
                        Ui = 1200;
                    case {3, 4}
                        Li = -0.55;
                        Ui = 0.55;
                end
                newGen(var) = Li + rand * (Ui - Li);
            end
            newGen(end) = f_x(newGen(1:end-1));

            %disp(newGen);
            % Verifica cual es mejor (reglas de DEB y restricciones)
            [c, ceq] = constraints(newGen);
            % Cumple con ser factible
            sum = 0;
            %sum = max(0,c(1)) + max(0,c(2)) + abs(ceq(1)) + abs(ceq(2)) + abs(ceq(3));
            sum = max(0,c(1)) + max(0,c(2));
            % if c(1) > 0
            %     sum = sum + c(1);
            % end
            % if c(2) > 0
            %     sum = sum + c(2);
            % end
            if abs(ceq(1)) > cero_gordo
                sum = sum + abs(ceq(1));
            end
            if abs(ceq(2)) > cero_gordo
                sum = sum + abs(ceq(2));
            end
            if abs(ceq(3)) > cero_gordo
                sum = sum + abs(ceq(3));
            end


            worst = max(HM(:, end));
            Xworst = HM(find(HM(:, end) == worst, 1),:);
            
            [c2, ceq2] = constraints(Xworst);
            % Cumple con ser factible
            sum2 = 0;
            %sum2 = max(0,c2(1)) + max(0,c2(2)) + abs(ceq2(1)) + abs(ceq2(2)) + abs(ceq2(3));
            sum2 = max(0,c2(1)) + max(0,c2(2));
            % if c2(1) > 0
            %     sum = sum + c2(1);
            % end
            % if c2(2) > 0
            %     sum = sum + c2(2);
            % end
            if abs(ceq2(1)) > cero_gordo
                sum2 = sum2 + abs(ceq2(1));
            end
            if abs(ceq2(2)) > cero_gordo
                sum2 = sum2 + abs(ceq2(2));
            end
            if abs(ceq2(3)) > cero_gordo
                sum2 = sum2 + abs(ceq2(3));
            end
             
            % fprintf('Sum:\n'); disp(c); disp(ceq);
            %disp(newGen);
            %fprintf('Sum1: %d, Sum2: %d\n', sum, sum2);
            %fprintf('Worst: %d\n', worst);

            if sum < sum2
            %if abs(ceq(1)) < cero_gordo && abs(ceq(2)) < cero_gordo && abs(ceq(3)) < cero_gordo && c(1) <= 0 && c(2) <= 0
                % Busca la mejor F.O.

                %if newGen(end) < worst
                if newGen(end) < worst
                    %HM(find(HM(:, end) == worst, 1),:) = newGen;
                    Xworst = newGen;
                    worst = newGen(end);
                    %fprintf('New Worst: %d\n', worst);
                    %disp(Xworst);
                    X(gen, 1:end-1) = Xworst;
                    X(gen, end) = sum;
                end
            end
        end
        %disp(HM);
    end
    
    %disp(X);
    less_sum = 1000000;
    for i=1:iteraciones
        if X(i,end) < less_sum && X(i,end-1) > 4800.0
            %disp(X(i,end-1));
            Xworst = X(i,1:end-1);
        end
    end
    % fprintf('Mínimo valor de F.O: \n'); disp(final_min);
    % fprintf('Valores de variables:\n\nx1->%d\nx2->%d\nx3->%d\nx4->%d\n', final(1), final(2), final(3), final(4));
    % for i=1:num_individuos
    %     [c_f, ceq_f] = constraints(X(i,:));
    %     if c_f(1) > 0 && c_f(2) > 0 && abs(ceq_f(1)) > cero_gordo && abs(ceq_f(2)) > cero_gordo && abs(ceq_f(3)) > cero_gordo
    %         X(i,:) = [];
    %     end
    % end
    % %disp(HM);
    %final_min = min(HM(:, end));
    %final = HM(find(HM(:, end) == final_min, 1),:);
    % disp(final);
    best_harmony = Xworst;
    %best_harmony = final;
    %fprintf('Final Harmony: \n');
    %disp(HM);
    %fprintf('Worst: '); disp(Xworst);
end