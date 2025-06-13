%Ítem [4] Obtener el torque máximo que puede soportar el motor modelado mediante las Ecs. (1-5) 
%(1-6) y (1-7) cuando se lo alimenta con 12V, graficando para 5 segundos de tiempo la velocidad 
%angular y corriente ia para establecer su valor máximo como para dimensionar dispositivos 
%electrónicos. 
%Universidad Nacional de Catamarca. Facultad de Tecnologia y Ciencias Aplicadas
%Alumno: Ocampos Gustavo Ismael
%M.U.:983
clc;
clear all;

% Definimos las constantes del sistema cuyos valores fueron proporcionados
% por la catedra
L = 366e-6; 
J = 5e-9; 
R = 55.6; 
Bm = 0; 
Ki = 6.49e-3; 
Km = 6.53e-3;

% Tiempo de simulación
t_euler = 1e-7;            % Paso de integración
tiempo = 0:t_euler:0.1;   % Simular 5 segundos(reducimos el tiempo para acelerar la compilacion)

% Inicialización de variables
tita_c = zeros(1, length(tiempo));
w_c = zeros(1, length(tiempo));
I_a_c = zeros(1, length(tiempo));
Torque = zeros(1, length(tiempo));

% Estado inicial (corriente, posición, velocidad angular)
X = [0; 0; 0];

% Matrices del sistema
A = [(-R/L) 0 (-Km/L); 
      0     0   1; 
     (Ki/J) 0 (-Bm/J)];
Bva = [1/L; 0; 0];

% Entrada de tensión constante
u = 12;

% Integración por método de Euler
for i = 2:length(tiempo)
    X_P= A * X + Bva * u;
    X = X + t_euler * X_P;

    I_a_c(i) = X(1);
    tita_c(i) = X(2);
    w_c(i) = X(3);
    Torque(i) = Ki * X(1); % Torque electromagnético
end

% Gráficos
subplot(3,1,1);
plot(tiempo, w_c, 'g');
xlabel('Tiempo [s]');
ylabel('Velocidad angular [rad/s]');
title('Velocidad angular del motor');

subplot(3,1,2);
plot(tiempo, I_a_c, 'b');
xlabel('Tiempo [s]');
ylabel('Corriente de armadura [A]');
title('Corriente de armadura del motor');

subplot(3,1,3);
plot(tiempo, Torque, 'r');
xlabel('Tiempo [s]');
ylabel('Torque electromagnético [Nm]');
title('Torque desarrollado por el motor');

% Valores máximos
[w_max, idx_w] = max(w_c);
[ia_max, idx_ia] = max(I_a_c);
[Tmax, idx_T] = max(Torque);

fprintf('Velocidad máxima: %.2f rad/s a t=%.6f s\n', w_max, tiempo(idx_w));
fprintf('Corriente máxima: %.2f A a t=%.6f s\n', ia_max, tiempo(idx_ia));
fprintf('Torque máximo: %.6f Nm a t=%.6f s\n', Tmax, tiempo(idx_T));
