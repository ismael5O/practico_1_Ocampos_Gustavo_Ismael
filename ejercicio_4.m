%Ítem [4] Obtener el torque máximo que puede soportar el motor modelado mediante las Ecs. (1-5) 
%(1-6) y (1-7) cuando se lo alimenta con 12V, graficando para 5 segundos de tiempo la velocidad 
%angular y corriente ia para establecer su valor máximo como para dimensionar dispositivos 
%electrónicos. 
%Universidad Nacional de Catamarca. Facultad de Tecnologia y Ciencias Aplicadas
%Alumno: Ocampos Gustavo Ismael
%M.U.:983
clc;
clear all;

 datos = xlsread('Curvas_Medidas_Motor_2024_v.xls');% Definimos las constantes del sistema
% cuyos valores fueron proporcionados por la catedra
L = 366e-6; 
J = 5e-9; 
R = 55.6; 
Bm = 0; 
Ki = 6.49e-3; 
Km = 6.53e-3;

% Tiempo de simulación
t_euler = 1e-7;            % Paso de integración
tiempo = 0:t_euler: 0.2;   % Simular 5 segundos(reducimos el tiempo para acelerar la compilacion)

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
Btl=[0; 0; -1/J];
% Entrada de tensión constante
u = 12;
for j=1: length(tiempo)
    
end
% Integración por método de Euler
for i = 2:length(tiempo)
    
    if(tiempo(i)>0.04 && tiempo(i)<0.07)
       Torque(i)=0.0014; 
    else
        Torque(i)=0;
        end
    %definimos nuestro vector de variables derivadas
  X_Punto=A*X+Bva*u+Btl*Torque(i);
  %Integramos para obtener X y definir el comportamiento del sistema
  X=X+t_euler*X_Punto ;
  
  %obtenemos los valores de las variables para cada iteracion
  I_a_c(i)=X(1);
  tita_c(i)=X(2);
  w_c(i)=X(3);
  Torque(i) = Ki * X(1); % Torque electromagnético
  
 
end

% Gráficos
Tl_max=max(I_a_c)*Ki;
subplot(4,1,1);
plot(tiempo,w_c,'r'),title('Velocidad Angular');
subplot(4,1,2) ,plot(tiempo,I_a_c),title('Corrientede Armadura');
subplot(4,1,3) ,plot(tiempo,tita_c),title('Angulo');
subplot(4,1,4), plot(tiempo,Torque),title('Torque');
