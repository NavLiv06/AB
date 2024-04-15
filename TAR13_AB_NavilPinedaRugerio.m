%---------------------------------------------------------
% Programa de algoritmo ABC
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear
rng(0,'twister');

abejas = [-5,5,5;
          -5,5,5];

reparticion = [10, 5;
               10, 2;
               30, 1];

ejecuciones = 30;
finals = zeros(ejecuciones,4);
fprintf('Algoritmo ABC...\n');

iteraciones = 10;
for i=1:ejecuciones
    finals(i, :) = abc_generator(abejas, 100, reparticion, 10, iteraciones);
end
fprintf('Mejores resultados:\n');
fprintf('    x1           x2             f(x)        Limite\n');
finals = sortrows(finals, size(finals, 2));
for i=1:size(finals,1)
    fprintf('%.10f       %.10f        %.10f        %d\n', finals(i,1), finals(i,2), finals(i,3), finals(i,4))
end

% Función objetivo
function f = f_x(x)
    f = x(1)^2 + x(2)^2;
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

function best_fa = abc_generator(abejas, num_abejas, reparticion, limlim, iteraciones)
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
    % Calcular nuevas velocidades
    for i=1:num_abejas
        % Buscar un k
        k = randi([1, num_abejas]);
        % Generar dos phi
        phi_1 = rand();
        phi_2 = rand();
        % Aplicar formula
        v1 = f_alimento(i, 1) + (phi_1 * (f_alimento(i, 1) - f_alimento(k, 1)));
        v2 = f_alimento(i, 2) + (phi_2 * (f_alimento(i, 2) - f_alimento(k, 2)));
        f_nueva = f_x([v1, v2]);
        % Si es menor cambiarla
        if f_nueva < f_alimento(i,end-1)
            f_alimento(i,1) = v1;
            f_alimento(i,2) = v2;
            f_alimento(i,end-1) = f_nueva;
            f_alimento(i,end) = 0;
            %disp(f_alimento(i,:));
        else
            f_alimento(i,end) = f_alimento(i,end) + 1;
        end
    end
    %fprintf('Velocidades actualizadas\n');
    %disp(f_alimento);

    g_best = [0, 0, 100000, 0];

    winners = floor(num_abejas/2);
    menores = zeros(winners, num_variables);
    for iter=1:iteraciones
        %%%%%%%%%%%%%%%%%%% SEGUNDO PROCESO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_alimento = sortrows(f_alimento, 3);
        %disp(f_alimento2);
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
                abejas = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(abejas,1)
                    k = randi([1, size(abejas,1)]);
                    phi_1 = rand();
                    phi_2 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    f_nueva = f_x([v1, v2]);
                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if f_nueva < f_alimento(j,end-1)
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,end-1) = round(f_nueva,6);
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        if g_best(end-1) > f_alimento(j,end-1)
                            %fprintf('F: %.16f\n',f_alimento(j,end-1));
                            g_best = f_alimento(j,:);
                        end
                        f_alimento(j,end) = 0;
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
                abejas = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(abejas,1)
                    k = randi([1, size(abejas,1)]);
                    phi_1 = rand();
                    phi_2 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    f_nueva = f_x([v1, v2]);
                    % Si alcanzo el limite de limites
                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if f_nueva < f_alimento(j,end-1)
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,end-1) = f_nueva;
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        if g_best(end-1) > f_alimento(j,end-1)
                            g_best = f_alimento(j,:);
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
                abejas = f_alimento(lim1:lim2, :);
                %disp(abejas);
                for i=1:size(abejas,1)
                    k = randi([1, size(abejas,1)]);
                    phi_1 = rand();
                    phi_2 = rand();
                    % Aplicar formula
                    v1 = f_alimento(j, 1) + (phi_1 * (f_alimento(j, 1) - f_alimento(k, 1)));
                    v2 = f_alimento(j, 2) + (phi_2 * (f_alimento(j, 2) - f_alimento(k, 2)));
                    f_nueva = f_x([v1, v2]);
                    if f_alimento(j,end) <= limlim
                        % Si es menor cambiarla
                        if f_nueva < f_alimento(j,end-1)
                            f_alimento(j,1) = v1;
                            f_alimento(j,2) = v2;
                            f_alimento(j,end-1) = f_nueva;
                            f_alimento(j,end) = 0;
                            %disp(f_alimento(i,:));
                        else
                            f_alimento(j,end) = f_alimento(j,end) + 1;
                        end
                    else
                        if g_best(end-1) > f_alimento(j,end-1)
                            g_best = f_alimento(j,:);
                        end
                        f_alimento(j,end) = 0;
                    end
                end
                lim1 = lim2 + 1;
                lim2 = lim2 + reparticion(3,2);
            end
    
        end
    end

    %disp(f_alimento);
    %disp(g_best);
    
    best_fa = g_best;
end