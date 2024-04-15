%---------------------------------------------------------
% Programa de algoritmo ABC
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear
rng(0,'twister');

abejas = [0,1200,2;
          0,1200,2;
         -0.55,0.55,2;
         -0.55,0.55, 2];

reparticion = [10, 5;
               10, 2;
               30, 1];

% reparticion = [20, 3;
%                20, 1;
%                10, 2];

% reparticion = [10, 3;
%                30, 2;
%                10, 1];

% reparticion = [1, 5;
%                1, 2;
%                3, 1];

% reparticion = [100, 5;
%                100, 2;
%                300, 1];

ejecuciones = 10;
finals = zeros(ejecuciones,6);
fprintf('Algoritmo ABC...\n');

%limite = 4 * 100;
limite = 6;

iteraciones = 100;
for i=1:ejecuciones
    finals(i, :) = abc_generator(abejas, 100, reparticion, limite, iteraciones);
end
fprintf('Mejores resultados:\n');
fprintf('      x1        x2        x3       x4       f(x)        Limite\n');
finals = sortrows(finals, 5);
disp(finals);
%for i=1:size(finals,1)
%    fprintf('%.3f       %.3f        %.3f        %.3f        %.3f        %d\n', finals(i,1), finals(i,2), finals(i,3), finals(i,4), finals(i,5), finals(i,6));
%end

% Función objetivo
function f = f_x(x)
    f = 3 * x(1) + 0.000001 * x(1)^3 + 2 * x(2) + (0.000002/3) * x(2)^3;
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

% Para restricciones
function [g, h] = restricciones(x)
    g = [-x(4) + x(3) - 0.55;
         -x(3) + x(4) - 0.55];
    
    h = [1000*(sin(-(x(3)) - 0.25)) + 1000*(sin(-(x(4)) - 0.25)) + 894.8 - (x(1));
           1000*(sin((x(3)) - 0.25)) + 1000*(sin((x(3)) - (x(4)) - 0.25)) + 894.8 - (x(2));
           1000*(sin((x(4)) - 0.25)) + 1000*(sin((x(4)) - (x(3)) - 0.25)) + 1294.8];
end


function best_fa = abc_generator(abejas, num_abejas, reparticion, limlim, iteraciones)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERADOR DE PARTICULAS ALEATORIA
    % EVALUACIÓN DE LA FUNCIÓN OBJETIVO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_variables = size(abejas, 1);
    f_alimento = zeros(num_abejas, num_variables+2);
    for i = 1:num_variables
        lower_bound = abejas(i, 1);
        upper_bound = abejas(i, 2);
        f_alimento(:, i) = lower_bound + (upper_bound - lower_bound) * rand(num_abejas, 1);
    end

    for i = 1:num_abejas
        f_alimento(i, end-1) = f_x(f_alimento(i, :));
        f_alimento(i, end) = 0;
    end
    %fprintf('Evaluación de las fuentes de alimento en la función objetivo\n');
    %disp(f_alimento);

    %%%%%%%%%%%%%%%%%%% PRIMER PROCESO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cero_gordo = 0.000001;
    % Calcular nuevas velocidades
    for i=1:num_abejas
        % Buscar un k
        k = randi([1, num_abejas]);
        % Generar dos phi
        phi_1 = rand();
        phi_2 = rand();
        phi_3 = rand();
        phi_4 = rand();
        % Aplicar formula
        v1 = f_alimento(i, 1) + (phi_1 * (f_alimento(i, 1) - f_alimento(k, 1)));
        v2 = f_alimento(i, 2) + (phi_2 * (f_alimento(i, 2) - f_alimento(k, 2)));
        v3 = f_alimento(i, 3) + (phi_3 * (f_alimento(i, 3) - f_alimento(k, 3)));
        v4 = f_alimento(i, 4) + (phi_4 * (f_alimento(i, 4) - f_alimento(k, 4)));

        X = [v1, v2, v3, v4];

        indx1 = find(X < abejas(:, 1)');
        indx2 = find(X > abejas(:, 2)');

        while ~isempty(indx1) || ~isempty(indx2)
            for l=1:num_variables
                if X(l) < abejas(l,1)
                    X(l) = 2 * abejas(l,1) - X(l);
                end
                if X(l) > abejas(l,2)
                    X(l) = 2 * abejas(l,2) - X(l);
                end
                indx1 = find(X < abejas(:, 1)');
                indx2 = find(X > abejas(:, 2)');
            end
        end

        f_nueva = f_x(X);

        
        % APLICAR RESTRICCIONES A LA NUEVA F.A.
        [g, h] = restricciones(X);
        sum2 = 0;
        sum2 = max(0,g(1)) + max(0,g(2));
        if abs(h(1)) > cero_gordo
            sum2 = sum2 + abs(h(1));
        end
        if abs(h(2)) > cero_gordo
            sum2 = sum2 + abs(h(2));
        end
        if abs(h(3)) > cero_gordo
            sum2 = sum2 + abs(h(3));
        end

        % APLICAR RESTRICCIONES A LA ANTIGUA F.A.
        [g2, h2] = restricciones(f_alimento(i,1:4));
        sum = 0;
        sum = max(0,g2(1)) + max(0,g2(2));
        if abs(h2(1)) > cero_gordo
            sum = sum + abs(h2(1));
        end
        if abs(h2(2)) > cero_gordo
            sum = sum + abs(h2(2));
        end
        if abs(h2(3)) > cero_gordo
            sum = sum + abs(h2(3));
        end


        % Si es menor cambiarla
        if sum2 < sum
            %fprintf('Sum: %f < %f\n',sum2, sum);
            f_alimento(i,1) = v1;
            f_alimento(i,2) = v2;
            f_alimento(i,3) = v3;
            f_alimento(i,4) = v4;
            f_alimento(i,end-1) = f_nueva;
            f_alimento(i,end) = 0;
            %disp(f_alimento(i,:));
        else
            f_alimento(i,end) = f_alimento(i,end) + 1;
        end
    end

    for i=1:num_abejas
        X = f_alimento(i,1:4);

        indx1 = find(X < abejas(:, 1)');
        indx2 = find(X > abejas(:, 2)');

        while ~isempty(indx1) || ~isempty(indx2)
            for l=1:num_variables
                if X(l) < abejas(l,1)
                    X(l) = 2 * abejas(l,1) - X(l);
                end
                if X(l) > abejas(l,2)
                    X(l) = 2 * abejas(l,2) - X(l);
                end
                indx1 = find(X < abejas(:, 1)');
                indx2 = find(X > abejas(:, 2)');
            end
        end

        f_alimento(i,1:4) = X;
        f_alimento(i,end-1) = f_x(X);
    end
    %fprintf('Velocidades actualizadas\n');
    %disp(f_alimento);

    g_best = [1000000, 1000000, 1000000, 1000000, 1000000, 0];
    %fprintf('Segundo proceso\n');

    winners = floor(num_abejas/2);
    menores = zeros(winners, num_variables);
    f_alimento = sortrows(f_alimento, 5);
    for iter=1:iteraciones
        %%%%%%%%%%%%%%%%%%% SEGUNDO PROCESO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %disp(f_alimento);
        % Ordenar las cinco menores
        menores = f_alimento(1:winners,:);
        %disp(menores);
        %fprintf('Size: %d\n', size(menores, 1));
        for j=1:winners
            lim1 = 1;
            lim2 = reparticion(1,2);
            abeja = 1;
            %%%%%%%%%%%%%%%%%%%%% PRIMERA REPARTICION %%%%%%%%%%%%%%%%%%%%%%%%%
            %fprintf('Lim: %d, %d\n', abeja, (abeja + reparticion(1,1)));
            for abj=abeja:abeja + reparticion(1,1) - 2
                %fprintf('Lims: %d,%d\n', lim1, lim2);
                bees = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(bees,1)
                    %k = randi([1, size(abejas,1)]);
                    k = randi([1, num_abejas]);
                    phi_1 = rand();
                    phi_2 = rand();
                    phi_3 = rand();
                    phi_4 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    v3 = f_alimento(j, 3) + (phi_3 * (f_alimento(j, 3) - f_alimento(k, 3)));
                    v4 = f_alimento(j, 4) + (phi_4 * (f_alimento(j, 4) - f_alimento(k, 4)));
                    X = [v1, v2, v3, v4];
            
                    f_nueva = f_x(X);


                    % APLICAR RESTRICCIONES A LA NUEVA F.A.
                    [g, h] = restricciones(X);
                    sum2 = 0;
                    sum2 = max(0,g(1)) + max(0,g(2));
                    if abs(h(1)) > cero_gordo
                        sum2 = sum2 + abs(h(1));
                    end
                    if abs(h(2)) > cero_gordo
                        sum2 = sum2 + abs(h(2));
                    end
                    if abs(h(3)) > cero_gordo
                        sum2 = sum2 + abs(h(3));
                    end
            
                    % APLICAR RESTRICCIONES A LA ANTIGUA F.A.
                    [g2, h2] = restricciones(f_alimento(i,1:4));
                    sum = 0;
                    sum = max(0,g2(1)) + max(0,g2(2));
                    if abs(h2(1)) > cero_gordo
                        sum = sum + abs(h2(1));
                    end
                    if abs(h2(2)) > cero_gordo
                        sum = sum + abs(h2(2));
                    end
                    if abs(h2(3)) > cero_gordo
                        sum = sum + abs(h2(3));
                    end

                    % APLICAR RESTRICCIONES A GBEST
                    [g3, h3] = restricciones(g_best(1:4));
                    sum3 = 0;
                    sum3 = max(0,g3(1)) + max(0,g3(2));
                    if abs(h3(1)) > cero_gordo
                        sum3 = sum3 + abs(h3(1));
                    end
                    if abs(h3(2)) > cero_gordo
                        sum3 = sum3 + abs(h3(2));
                    end
                    if abs(h3(3)) > cero_gordo
                        sum3 = sum3 + abs(h3(3));
                    end


                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if sum2 < sum
                            %fprintf('Sum: %f\n',sum2);
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,3) = v3;
                            f_alimento(j,4) = v4;
                            f_alimento(j,end-1) = f_nueva;
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            %fprintf('No hay más SVR: %f < %f\n', sum2, sum);
                            %f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        %if g_best(end-1) > f_alimento(j,end-1)
                        if sum3 > sum
                            %fprintf('Sum: %f %f\n',sum2, sum3);
                            %disp(g_best);
                            %disp(f_alimento);
                            %fprintf('F: %.16f\n',f_alimento(j,end-1));
                            g_best = f_alimento(j,:);
                            %disp(g_best);
                        else
                            %fprintf('No se encontró una mejor que g_best');
                            %disp(g_best);
                            %fprintf('\n');
                        end
                        %f_alimento(j,:) = [0,0,0,0,0,0];
                        
                        %f_alimento(j,end) = 0;
                    end
                end
                lim1 = lim2 + 1;
                lim2 = lim2 + reparticion(1,2);
                %fprintf('Lims-a: %d, %d\n', lim1, lim2);
            end
    
            %%%%%%%%%%%%%%%%%%%%% SEGUNDA REPARTICION %%%%%%%%%%%%%%%%%%%%%%%%%
            abeja = abj;
            %fprintf('Lim: %d, %d\n', abeja, (abeja + reparticion(2,1)));
            lim1 = lim2 + 1;
            lim2 = lim2 + reparticion(2,2);
            %fprintf('Lims-bp: %d, %d\n', lim1, lim2);
            for abj=abeja:abeja + reparticion(2,1) - 2
                %fprintf('Lims: %d,%d\n', lim1, lim2);
                bees = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(bees,1)
                    %k = randi([1, size(abejas,1)]);
                    k = randi([1, num_abejas]);
                    phi_1 = rand();
                    phi_2 = rand();
                    phi_3 = rand();
                    phi_4 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    v3 = f_alimento(j, 3) + (phi_3 * (f_alimento(j, 3) - f_alimento(k, 3)));
                    v4 = f_alimento(j, 4) + (phi_4 * (f_alimento(j, 4) - f_alimento(k, 4)));
                    X = [v1, v2, v3, v4];
            
                    f_nueva = f_x(X);


                    % APLICAR RESTRICCIONES A LA NUEVA F.A.
                    [g, h] = restricciones(X);
                    sum2 = 0;
                    sum2 = max(0,g(1)) + max(0,g(2));
                    if abs(h(1)) > cero_gordo
                        sum2 = sum2 + abs(h(1));
                    end
                    if abs(h(2)) > cero_gordo
                        sum2 = sum2 + abs(h(2));
                    end
                    if abs(h(3)) > cero_gordo
                        sum2 = sum2 + abs(h(3));
                    end
            
                    % APLICAR RESTRICCIONES A LA ANTIGUA F.A.
                    [g2, h2] = restricciones(f_alimento(i,1:4));
                    sum = 0;
                    sum = max(0,g2(1)) + max(0,g2(2));
                    if abs(h2(1)) > cero_gordo
                        sum = sum + abs(h2(1));
                    end
                    if abs(h2(2)) > cero_gordo
                        sum = sum + abs(h2(2));
                    end
                    if abs(h2(3)) > cero_gordo
                        sum = sum + abs(h2(3));
                    end

                    % APLICAR RESTRICCIONES A GBEST
                    [g3, h3] = restricciones(g_best(1:4));
                    sum3 = 0;
                    sum3 = max(0,g3(1)) + max(0,g3(2));
                    if abs(h3(1)) > cero_gordo
                        sum3 = sum3 + abs(h3(1));
                    end
                    if abs(h3(2)) > cero_gordo
                        sum3 = sum3 + abs(h3(2));
                    end
                    if abs(h3(3)) > cero_gordo
                        sum3 = sum3 + abs(h3(3));
                    end


                    % Si alcanzo el limite de limites
                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if sum2 < sum
                            %fprintf('Sum: %f\n',sum2);
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,3) = v3;
                            f_alimento(j,4) = v4;
                            f_alimento(j,end-1) = f_nueva;
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        if sum3 > sum
                        %if g_best(end-1) > f_alimento(j,end-1)
                            g_best = f_alimento(j,:);
                        else
                            %fprintf('No se encontró una mejor que g_best');
                            %disp(g_best);
                            %fprintf('\n');
                        end
                        f_alimento(j,end) = 0;
                    end
                end
                lim1 = lim2 + 1;
                lim2 = lim2 + reparticion(2,2);
                %fprintf('Lims-b: %d, %d\n', lim1, lim2);
            end
    
            %%%%%%%%%%%%%%%%%%%%% TERCERA REPARTICION %%%%%%%%%%%%%%%%%%%%%%%%%
            abeja = abj + 1;
            lim1 = lim2 + 1;
            lim2 = lim2 + reparticion(3,2);
            %fprintf('Lim: %d, %d\n', abeja, abeja + reparticion(3,1));
            %fprintf('Lims-c: %d, %d\n', lim1, lim2);
            for abj=abeja:abeja + reparticion(3,1) - 1
                %fprintf('Lims: %d,%d\n', lim1, lim2);
                bees = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(bees,1)
                    %k = randi([1, size(abejas,1)]);
                    k = randi([1, num_abejas]);
                    phi_1 = rand();
                    phi_2 = rand();
                    phi_3 = rand();
                    phi_4 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    v3 = f_alimento(j, 3) + (phi_3 * (f_alimento(j, 3) - f_alimento(k, 3)));
                    v4 = f_alimento(j, 4) + (phi_4 * (f_alimento(j, 4) - f_alimento(k, 4)));
                    X = [v1, v2, v3, v4];
            
                    f_nueva = f_x(X);


                    % APLICAR RESTRICCIONES A LA NUEVA F.A.
                    [g, h] = restricciones(X);
                    %fprintf('Res:'); disp(g); disp(h);
                    sum2 = 0;
                    sum2 = max(0,g(1)) + max(0,g(2));
                    if abs(h(1)) > cero_gordo
                        sum2 = sum2 + abs(h(1));
                    end
                    if abs(h(2)) > cero_gordo
                        sum2 = sum2 + abs(h(2));
                    end
                    if abs(h(3)) > cero_gordo
                        sum2 = sum2 + abs(h(3));
                    end
            
                    % APLICAR RESTRICCIONES A LA ANTIGUA F.A.
                    [g2, h2] = restricciones(f_alimento(i,1:4));
                    sum = 0;
                    sum = max(0,g2(1)) + max(0,g2(2));
                    if abs(h2(1)) > cero_gordo
                        sum = sum + abs(h2(1));
                    end
                    if abs(h2(2)) > cero_gordo
                        sum = sum + abs(h2(2));
                    end
                    if abs(h2(3)) > cero_gordo
                        sum = sum + abs(h2(3));
                    end

                    % APLICAR RESTRICCIONES A GBEST
                    [g3, h3] = restricciones(g_best(1:4));
                    sum3 = 0;
                    sum3 = max(0,g3(1)) + max(0,g3(2));
                    if abs(h3(1)) > cero_gordo
                        sum3 = sum3 + abs(h3(1));
                    end
                    if abs(h3(2)) > cero_gordo
                        sum3 = sum3 + abs(h3(2));
                    end
                    if abs(h3(3)) > cero_gordo
                        sum3 = sum3 + abs(h3(3));
                    end


                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if sum2 < sum
                            %fprintf('Sum: %f\n',sum2);
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,3) = v3;
                            f_alimento(j,4) = v4;
                            f_alimento(j,end-1) = f_nueva;
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        if sum3 > sum
                        %if g_best(end-1) > f_alimento(j,end-1)
                            g_best = f_alimento(j,:);
                        else
                            %fprintf('No se encontró una mejor que g_best');
                            %disp(g_best);
                            %fprintf('\n');
                        end
                        f_alimento(j,end) = 0;
                    end
                end
                lim1 = lim2 + 1;
                lim2 = lim2 + reparticion(3,2);
            end
    
        end

        for i=1:num_abejas
            X = f_alimento(i,1:4);
    
            indx1 = find(X < abejas(:, 1)');
            indx2 = find(X > abejas(:, 2)');
    
            while ~isempty(indx1) || ~isempty(indx2)
                for l=1:num_variables
                    if X(l) < abejas(l,1)
                        X(l) = 2 * abejas(l,1) - X(l);
                    end
                    if X(l) > abejas(l,2)
                        X(l) = 2 * abejas(l,2) - X(l);
                    end
                    indx1 = find(X < abejas(:, 1)');
                    indx2 = find(X > abejas(:, 2)');
                end
            end
    
            f_alimento(i,1:4) = X;
            f_alimento(i,end-1) = f_x(X);
        end
    end

    %disp(f_alimento);
    %disp(g_best);
    
    best_fa = g_best;
end