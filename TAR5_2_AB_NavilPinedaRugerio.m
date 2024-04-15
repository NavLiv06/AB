%---------------------------------------------------------
% Programa para generar numeros aleatorios basado en XORShift
% Navil Pineda Rugerio
% Inteligencia Artificial
% 5to Semestre
%---------------------------------------------------------

clc
clear

seed = 45678; % Semilla
num = 20; % Numeros a generar

% Numeros pseudoaleatorios entre 0 y 1
%random_numbers = xorshift(seed, num);

% Para convertir a otro rango diferente de 0 y 1
a = 500;
b = 1000;
random_numbers = a + (b - a) * xorshift(seed, num);
random_numbers = round(random_numbers);
disp(random_numbers);

% Generador de numeros de 32 bits
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