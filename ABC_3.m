num_sources = 100; % Número de fuentes de alimento
num_bees = 20; % Número de abejas
num_iterations = 100; % Número de iteraciones
limit = -5; % Límite superior e inferior para las posiciones iniciales de las abejas
step_size = 0.5; % Tamaño del paso

% Inicialización aleatoria de las posiciones iniciales de las fuentes de alimento
source_positions = limit * (rand(num_sources, 2) - 0.5);

% Asignación de abejas a la mitad de fuentes de alimento
bees_assigned = zeros(num_sources, 1);
for i = 1:num_sources   
    if i <= 10 % A 10/50 se les asignan 5
        bees_assigned(i) = 5;
    elseif i <= 20 % A 20/50 se les asignan 2
        bees_assigned(i) = 2;
    else % A 30/50 se les asignan 1
        bees_assigned(i) = 1;
    end
end

for iter = 1:num_iterations
    % Evaluar la función objetivo en las posiciones actuales de las fuentes de alimento
    fitness = sum(source_positions.^2, 2); % f(x1, x2) = x1^2 + x2^2
    fprintf('Fitness');
    disp(fitness);
    
    % Encontrar las mejores fuentes de alimento ordenandolas
    [~, sorted_indices] = sort(fitness);
    selected_sources = sorted_indices(1:(num_sources/2)); % Seleccionar las mejores fuentes (la mitad)
    
    % Actualizar las posiciones de las abejas asignadas a las fuentes seleccionadas
    for i = 1:length(selected_sources)
        num_bees_assigned = bees_assigned(selected_sources(i));
        for j = 1:num_bees_assigned
            % Mover las abejas hacia las fuentes seleccionadas
            source_position = source_positions(selected_sources(i), :);
            
            % Generar una nueva posición para cada abeja basada en la posición de la fuente seleccionada
            new_position = source_positions(j, :) + step_size * (source_position - source_positions(j, :)) * rand();
            
            % Verificar límites superior e inferior
            new_position(new_position > limit) = limit;
            new_position(new_position < -limit) = -limit;
            
            % Reemplazar la posición anterior si la nueva posición es mejor
            if sum((new_position - source_position).^2) < sum((positions(j, :) - source_position).^2)
                positions(j, :) = new_position;
            end
        end
    end
end

% Mostrar la mejor posición encontrada y su valor de función objetivo
disp('Mejor posición encontrada:');
disp(best_position);
disp('Valor de la función objetivo en la mejor posición:');
disp(best_fitness);
