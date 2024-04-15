clc;
clear;

rng(0,'twister');

abejas = [0,1200,2;
          0,1200,2;
         -0.55,0.55,2;
         -0.55,0.55,2];

reparticion = [10, 5;
               10, 2;
               30, 1];

ejecuciones = 10;
limlim = 10;
iteraciones = 1000;

for i = 1:ejecuciones
    % Inicialización de las fuentes de alimento
    f_alimento = initialize_food_sources(abejas, 100);
    %disp(f_alimento(10,:));
    
    % Fase de abejas empleadas y observadoras
    f_alimento = employed_bees_phase(f_alimento, abejas, reparticion, limlim, iteraciones);
    
    % Fase de abejas exploradoras
    f_alimento = scout_bees_phase(f_alimento, abejas, limlim, 4);

    f_alimento = verify_constraints(f_alimento, abejas);
    
    % Mostrar resultados o hacer lo que necesites con 'f_alimento'
    f_alimento_ordenado = sortrows(f_alimento, size(f_alimento, 2) - 1);
    umbral = 5126.49;
    f_alimento_ordenado = eliminar_soluciones_umbral(f_alimento_ordenado, umbral);
    %f_alimento_ordenado = f_alimeto_ordenado.*10;
    if i == ejecuciones
        formatSpec = '%.6f ';

        disp('Resultados:');
        for i = 1:size(f_alimento_ordenado, 1)
            fprintf(formatSpec, f_alimento_ordenado(i, 1:end-1));
            fprintf('\n');
        end
    end
end

function f_alimento = eliminar_soluciones_umbral(f_alimento, umbral)
    i = size(f_alimento, 1);
    while i > 0
        valor_objetivo = f_alimento(i, end - 1); % Obtener valor de la función objetivo
        
        if valor_objetivo < umbral
            f_alimento(i, :) = []; % Eliminar fila
        end
        
        i = i - 1; % Moverse hacia atrás
    end
end


function f = f_x(x)
    f = 3 * x(1) + 0.000001 * x(1)^3 + 2 * x(2) + (0.000002/3) * x(2)^3;
end

function v = apply_boundaries(v, bounds)
    for i = 1:numel(v)
        if v(i) < bounds(i, 1)
            v(i) = bounds(i, 1);
        elseif v(i) > bounds(i, 2)
            v(i) = bounds(i, 2);
        end
    end
end

function [g, h] = restricciones(x)
    g = [-x(4) + x(3) - 0.55;
         -x(3) + x(4) - 0.55];
    
    h = [1000 * (sin(-x(3) - 0.25)) + 1000 * (sin(-x(4) - 0.25)) + 894.8 - x(1);
         1000 * (sin(x(3) - 0.25)) + 1000 * (sin(x(3) - x(4) - 0.25)) + 894.8 - x(2);
         1000 * (sin(x(4) - 0.25)) + 1000 * (sin(x(4) - x(3) - 0.25)) + 1294.8];
end

function f_alimento = initialize_food_sources(abejas, num_abejas)
    num_variables = size(abejas, 1);
    f_alimento = zeros(num_abejas, num_variables + 2);
    
    for i = 1:num_variables
        lower_bound = abejas(i, 1);
        upper_bound = abejas(i, 2);
        f_alimento(:, i) = lower_bound + (upper_bound - lower_bound) * rand(num_abejas, 1);
    end

    for i = 1:num_abejas
        f_alimento(i, end-1) = f_x(f_alimento(i, :));
        f_alimento(i, end) = 0;
    end
end

function f_alimento = employed_bees_phase(f_alimento, abejas, reparticion, limlim, iteraciones)
    cero_gordo = 0.000001;
    num_variables = size(abejas, 1);
    num_abejas = size(f_alimento, 1);
    winners = floor(num_abejas / 2);

    for iter = 1:iteraciones
        for j = 1:num_abejas
            k = randi([1, num_abejas]);
            phi = rand(1, num_variables);
            v = f_alimento(j, 1:num_variables) + phi .* (f_alimento(j, 1:num_variables) - f_alimento(k, 1:num_variables));
            v = apply_boundaries(v, abejas); % Función para aplicar límites
            f_nueva = f_x(v);
            
            % Aplicar restricciones
            [g, h] = restricciones(v);
            sum2 = sum(max(0, g)) + sum(abs(h(h > cero_gordo)));
            
            [g2, h2] = restricciones(f_alimento(j, 1:num_variables));
            sum1 = sum(max(0, g2)) + sum(abs(h2(h2 > cero_gordo)));
            
            if sum2 < sum1
                f_alimento(j, 1:num_variables) = v;
                f_alimento(j, end-1) = f_nueva;
                f_alimento(j, end) = 0;
            else
                f_alimento(j, end) = f_alimento(j, end) + 1;
            end
        end
    end
end

function f_alimento = scout_bees_phase(f_alimento, abejas, limlim, num_variables)
    num_abejas = size(f_alimento, 1);
    cero_gordo = 0.000001;
    
    for i = 1:num_abejas
        if f_alimento(i, end) > limlim
            v = abejas(:, 1) + rand(size(abejas, 1), 1) .* (abejas(:, 2) - abejas(:, 1));
            f_nueva = f_x(v);
            
            [g, h] = restricciones(v);
            sum2 = sum(max(0, g)) + sum(abs(h(h > cero_gordo)));
            
            [g2, h2] = restricciones(f_alimento(i, 1:num_variables));
            sum1 = sum(max(0, g2)) + sum(abs(h2(h2 > cero_gordo)));
            
            if sum2 < sum1
                f_alimento(i, 1:num_variables) = v;
                f_alimento(i, end-1) = f_nueva;
                f_alimento(i, end) = 0;
            else
                f_alimento(i, end) = f_alimento(i, end) + 1;
            end
        end
    end
end

function [adjusted_solution, valid] = adjust_solution(solution)
    valid = true;
    adjusted_solution = solution;

    % Restricciones
    g_1 = -solution(4) + solution(3) - 0.55;
    g_2 = -solution(3) + solution(4) - 0.55;
    h_3 = 1000 * sin(-solution(3) - 0.25) + 1000 * sin(-solution(4) - 0.25) + 894.8 - solution(1);
    h_4 = 1000 * sin(solution(3) - 0.25) + 1000 * sin(solution(3) - solution(4) - 0.25) + 894.8 - solution(2);
    h_5 = 1000 * sin(solution(4) - 0.25) + 1000 * sin(solution(4) - solution(3) - 0.25) + 1294.8;

    % Verificar restricciones y ajustar si es necesario
    if g_1 > 0 || g_2 > 0 || abs(h_3) > 0 || abs(h_4) > 0 || abs(h_5) > 0
        valid = false;
        
        adjusted_solution = solution;
    end
end


% function f_alimento = verify_constraints(f_alimento)
%     cero_gordo = 0.000001;
%     num_variables = size(f_alimento, 2) - 2; % Obtener el número de variables
% 
%     for i = 1:size(f_alimento, 1)
%         % Verificar restricciones para cada solución
%         v = f_alimento(i, 1:num_variables);
%         [g, h] = restricciones(v);
%         sum_g = sum(max(0, g));
%         sum_h = sum(abs(h(h > cero_gordo)));
% 
%         % Si alguna restricción no se cumple, ajustar la solución
%         if sum_g > 0 || sum_h > 0
%             disp('¡Una solución no cumple con las restricciones!');
%             % Aquí puedes tomar acciones como ajustar la solución o marcarla como inválida
%         end
%     end
% end
% 

function f_alimento = verify_constraints(f_alimento, abejas)
    for i = 1:size(f_alimento, 1)
        v = f_alimento(i, 1:end - 2); % Obtener la solución
        [adjusted_solution, valid] = adjust_solution(v);
        
        if ~valid
            %disp('¡Se ha ajustado una solución!');
            % Reemplazar la solución no válida con la ajustada
            f_alimento(i, 1:end - 2) = adjusted_solution;
        end
    end
end
