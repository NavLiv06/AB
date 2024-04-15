% Configuración inicial
width = 1000;
height = 1000;
bg = [25, 25, 25];
nxC = 25;
nyC = 25;
dimCW = width / nxC;
dimCH = height / nyC;
gameState = zeros(nxC, nyC);
iteraciones = 100;

% Automata palo
gameState(6, 4:6) = 1;

% Automata movil
gameState(22, 22:23) = 1;
gameState(21, [23, 21:23]) = 1;
gameState(20, 23) = 1;

% Control de la ejecución del juego
%pauseExec = false;

%while true
for i=1:iteraciones
    newGameState = gameState;

    % Dibujar celdas
    clf;
    hold on;
    
    for y = 1:nyC
        for x = 1:nxC
            % Calcular el número de vecinos cercanos
            n_neigh = sum([gameState(mod(x-2, nxC)+1, mod(y-2, nyC)+1), ...
                           gameState(mod(x-1, nxC)+1, mod(y-2, nyC)+1), ...
                           gameState(mod(x, nxC)+1, mod(y-2, nyC)+1), ...
                           gameState(mod(x-2, nxC)+1, mod(y-1, nyC)+1), ...
                           gameState(mod(x+1, nxC)+1, mod(y-1, nyC)+1), ...
                           gameState(mod(x-2, nxC)+1, mod(y, nyC)+1), ...
                           gameState(mod(x, nxC)+1, mod(y, nyC)+1), ...
                           gameState(mod(x+1, nxC)+1, mod(y, nyC)+1)]);
                
            % Regla 1: Una célula muere con exactamente 3 vecinas vivas, "revive"
            if gameState(x, y) == 0 && n_neigh == 3
                newGameState(x, y) = 1;
            % Regla 2: Una célula viva con menos de 2 o más de 3 vecinas vivas, "muere"
            elseif gameState(x, y) == 1 && (n_neigh < 2 || n_neigh > 3)
                newGameState(x, y) = 0;
            end
            
            % Creamos el polígono de cada celda a dibujar
            poly = [((x-1) * dimCW), (y-1) * dimCH; ...
                    (x * dimCW), (y-1) * dimCH; ...
                    (x * dimCW), y * dimCH; ...
                    (x-1) * dimCW, y * dimCH];
            
            % Y dibujamos la celda para cada par de x e y
            if newGameState(x, y) == 0
                fill(poly(:, 1), poly(:, 2), [128, 128, 128] / 255);
            else
                fill(poly(:, 1), poly(:, 2), [255, 255, 255] / 255);
            end
        end
    end
    
    hold off;
    axis off;
    %pause(0.1);
    
    % Actualizamos el estado del juego
    gameState = newGameState;
end
