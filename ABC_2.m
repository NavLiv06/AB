num_bees = 20; % Número de abejas
num_iterations = 100; % Número de iteraciones
limit = 100; % Límite superior e inferior para las posiciones iniciales de las abejas
step_size = 0.5; % Tamaño del paso

% Inicialización aleatoria de las posiciones iniciales de las abejas
positions = limit * (rand(num_bees, 2) - 0.5);

for iter = 1:num_iterations
    % Evaluar la función objetivo en las posiciones actuales de las abejas
    fitness = sum(positions.^2, 2); % f(x1, x2) = x1^2 + x2^2
    
    % Encontrar la mejor posición (menor valor) de entre todas las abejas
    [best_fitness, best_idx] = min(fitness);
    best_position = positions(best_idx, :);
    
    % Actualizar las posiciones de las abejas
    for i = 1:num_bees
        % Generar una nueva posición para cada abeja basada en la mejor posición encontrada
        new_position = positions(i, :) + step_size * (randn(1, 2) * 2 - 1); % Movimiento aleatorio
        
        % Verificar límites superior e inferior
        new_position(new_position > limit) = limit;
        new_position(new_position < -limit) = -limit;
        
        % Reemplazar la posición anterior si la nueva posición es mejor
        if sum(new_position.^2) < fitness(i)
            positions(i, :) = new_position;
        end
    end
end

% Mostrar la mejor posición encontrada y su valor de la función objetivo
disp('Mejor posición encontrada:');
disp(best_position);
disp('Valor de la función objetivo en la mejor posición:');
disp(best_fitness);