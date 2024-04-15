step_size = 0.5; % Tamaño del paso


num_fuentes = 100;
num_bees = 20;
num_iter = 1000;
limites = 5;
ejecuciones = 30;

% Ejecucion para 30 veces
for e=1:ejecuciones
    source_positions = limites * (rand(num_fuentes, 2) - 0.5);
    %disp(source_positions);
    
    % Asignación de abejas a la mitad de fuentes de alimento
    bees_assigned = zeros(num_fuentes, 1);
    for i = 1:num_fuentes
        if i <= 10 % A 10/50 se les asignan 5
            bees_assigned(i) = 5;
        elseif i <= 20 % A 20/50 se les asignan 2
            bees_assigned(i) = 2;
        else % A 30/50 se les asignan 1
            bees_assigned(i) = 1;
        end
    end
    
    best_fitness = inf;
    best_position = zeros(1, 2);
    
    for iter = 1:num_iter
        % F.O.
        fitness = sum(source_positions.^2, 2); % f(x1, x2) = x1^2 + x2^2
        %fprintf('Fitness');
        %disp(fitness);

        %fprintf('%.6f %.6f  %.6f\n', source_positions(i,1), source_positions(i,2), fitness);
        
        
        [~, sorted_indices] = sort(fitness);
        selected_sources = sorted_indices(1:(num_fuentes/2)); % Seleccionar las mejores fuentes (la mitad)
        
        % Actualizar posicion
        for i = 1:length(selected_sources)
            num_bees_assigned = bees_assigned(selected_sources(i));
            for j = 1:num_bees_assigned
                source_position = source_positions(selected_sources(i), :);
                new_position = source_positions(j, :) + step_size * (source_position - source_positions(j, :)) * rand(); 
                % Verificar límites superior e inferior
                new_position(new_position > limites) = limites;
                new_position(new_position < -limites) = -limites;
                % Reemplazar la posición anterior si la nueva posición es mejor
                if sum((new_position - source_position).^2) < sum((source_positions(j, :) - source_position).^2)
                    source_positions(j, :) = new_position;
                    %fprintf('%d', sum((new_position).^2));
                    % Gbest
                    if sum((new_position).^2) < best_fitness
                        best_fitness = sum((new_position).^2);
                        %disp(best_fitness);
                        best_position = new_position;
                    end
                end
            end
        end
    end
    
    % Mostrar la mejor posición encontrada y su valor de función objetivo
    %disp('Mejor fuente encontrada (x1, x2): ');
    fprintf('%.6f %.6f  %.6f\n', best_position(1), best_position(2), best_fitness);
    %disp(best_position);
    %disp('Mejor F.O.');
    %disp(best_fitness);
end