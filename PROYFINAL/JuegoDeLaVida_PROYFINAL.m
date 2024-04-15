% Parámetros
nxC = 50;
nyC = 50;
width = 1000;
height = 1000;
dimCW = width / nxC;
dimCH = height / nyC;
num_iteraciones = 1000;

% Estado inicial de las celdas
gameState = zeros(nxC, nyC);

% Parpadeador (oscilador)
%gameState(6:8, 4) = 1;

% Automata movil
% gameState(22, 22:23) = 1;
% gameState(23, 23) = 1;
% gameState(22, 24) = 1;
% gameState(21, 24) = 1;
% gameState(20, 24) = 1;

% Glider (Planeador)
% gameState(25,25:27) = 1;
% gameState(24,27) = 1;
% gameState(23,26) = 1;

% Naves espaciales
% gameState(25:27,25) = 1;
% gameState(27,24) = 1;
% gameState(26,23) = 1;

% Gosper Glider Gun
gameState(6:7,2:3) = 1;
gameState(6:8,12) = 1;
gameState(6:8,12) = 1;
gameState(5,13) = 1;
gameState(9,13) = 1;
gameState(4,14:15) = 1;
gameState(10,14:15) = 1;
gameState(7,16) = 1;
gameState(5,17) = 1;
gameState(9,17) = 1;
gameState(6:8,18) = 1;
gameState(7,19) = 1;
gameState(4:6,22:23) = 1;
gameState(3,24) = 1;
gameState(7,24) = 1;
gameState(2:3,26) = 1;
gameState(7:8,26) = 1;
gameState(4:5,36:37) = 1;
gameState(12:13,25) = 1;
gameState(12:13,26) = 1;

% Segundo Glider Gun
gameState(15:16,46:47) = 1;
gameState(16,50) = 1;
gameState(15,51) = 1;
gameState(17,51) = 1;
gameState(14,52) = 1;
gameState(18,52) = 1;
gameState(14:18,53) = 1;
gameState(13:14,54) = 1;
gameState(18:19,54) = 1;
gameState(14:18,55) = 1;
gameState(15:17,56) = 1;
gameState(16,57) = 1;
gameState(20:21, 59) = 1;
gameState(19, 60) = 1;
gameState(21, 60) = 1;
gameState(21, 61) = 1;
gameState(17:19, 68) = 1;
gameState(16, 69) = 1;
gameState(20, 69) = 1;
gameState(15, 71:72) = 1;
gameState(21, 71:72) = 1;
gameState(16, 72) = 1;
gameState(20, 72) = 1;
gameState(18, 75) = 1;
gameState(17, 76:77) = 1;
gameState(19, 76:80) = 1;
gameState(16, 80) = 1;
gameState(17:18, 81) = 1;

% Estado inicial aleatorio
porcentaje = 100;
%gameState = randi([0, 1], nxC, nyC);
% total_celulas = (nxC * nyC) * porcentaje / 100;
% fprintf('Total de celulas:'); disp(nxC*nyC);
% fprintf('%d por ciento de celulas:', porcentaje); disp(total_celulas);
% celdas_si = rand(nxC, nyC) * 100;
% gameState(celdas_si <= porcentaje) = 1;


% Crear la ventana
figure;
title('Juego de la Vida.');
axis off;
axis equal;
xlim([0, width]);
ylim([0, height]);

imagesc(gameState, [0, 1]);
colormap([0.2, 0.2, 0.2; 1, 1, 1]);
axis off;

pause(2);

for iter=1:num_iteraciones
    % % Copiamos el estado actual
    % newGameState = gameState;
    % 
    % % Dibujamos el estado actual
    % for x = 1:nxC
    %     for y = 1:nyC
    %         % Calculamos el número de vecinos cercanos
    %         % Aplicamos el modulo para la forma toroidal
    %         n_neigh = gameState(mod(x-2, nxC)+1, mod(y-2, nyC)+1) + ...
    %                    gameState(mod(x-1, nxC)+1, mod(y-2, nyC)+1) + ...
    %                    gameState(mod(x  , nxC)+1, mod(y-2, nyC)+1) + ...
    %                    gameState(mod(x-2, nxC)+1, mod(y-1, nyC)+1) + ...
    %                    gameState(mod(x  , nxC)+1, mod(y-1, nyC)+1) + ...
    %                    gameState(mod(x-2, nxC)+1, mod(y  , nyC)+1) + ...
    %                    gameState(mod(x-1, nxC)+1, mod(y  , nyC)+1) + ...
    %                    gameState(mod(x  , nxC)+1, mod(y  , nyC)+1);
    % 
    %         % Reglas del juego
    %         % Regla 1: Una celula muerra con exactamente 3 vecinas vivas, "revive"
    %         if gameState(x, y) == 0 && n_neigh == 3
    %             newGameState(x, y) = 1;
    %         % Regla 2: Una celula viva con menos de 2 o mas de 3 vecinas vivas, "muere" (por abaondo o sobrepoblacion)
    %         elseif gameState(x, y) == 1 && (n_neigh < 2 || n_neigh > 3)
    %             newGameState(x, y) = 0;
    %         end
    % 
    %         % Dibujamos las celdas
    %         x_coords = [dimCW*(x-1), dimCW*x, dimCW*x, dimCW*(x-1)];
    %         y_coords = [dimCH*(y-1), dimCH*(y-1), dimCH*y, dimCH*y];
    % 
    %         % Coloreamos, si es celula muerta se clorea en gris
    %         if newGameState(x, y) == 0
    %             patch(x_coords, y_coords, [0.5, 0.5, 0.5]);
    %         % Coloreamos, si es celula viva se clorea en blanco
    %         else
    %             patch(x_coords, y_coords, [1, 1, 1]);
    %         end
    %     end
    % end
    % 
    % % Actualizar el estado
    % gameState = newGameState;
    % pause(0.05);
    % Calcular vecindario utilizando la convolución
    neighbors = conv2(double(gameState), ones(3), 'same') - double(gameState);

    % Aplicar reglas del juego de la vida
    newGameState = gameState;
    newGameState(gameState == 0 & neighbors == 3) = 1;
    newGameState(gameState == 1 & (neighbors < 2 | neighbors > 3)) = 0;

    % Mostrar el estado
    imagesc(newGameState, [0, 1]);
    colormap([0.1, 0.1, 0.1; 1, 1, 1]);
    axis off;

    % Actualizar el estado
    gameState = newGameState;
    
    % Pausa
    pause(0.05);
end