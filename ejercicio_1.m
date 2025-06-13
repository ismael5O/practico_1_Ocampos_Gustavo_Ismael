%Ítem [1] Asignar valores a R=47?, L=1?Hy, y C=100nF. Obtener simulaciones que permitan 
%estudiar la dinámica del sistema, con una entrada de tensión escalón de 12V, que cada 1ms cambia 
%de signo.
%Universidad Nacional de Catamarca. Facultad de Tecnologia y Ciencias Aplicadas
%Alumno: Ocampos Gustavo Ismael
%M.U.:983
clear all;
clc;

% ingresamos los parametros del circuito brindados por la catedra
R = 47;           
L = 1e-6;         
C = 100e-9;       

% ingresamos las matrices que contienen los coeficientes del circuito
Mat_A = [-R/L, -1/L; 1/C, 0];
Mat_B = [1/L; 0];
Mat_C = [R, 0];  % la salida es la tensión en la resistencia

% Análisis de constantes de tiempo (autovalores de A)
Constantes_Tiempo = eig(Mat_A);
tau_min = min(abs(1 ./ Constantes_Tiempo));  % ajustamos el tiempo mas rapido
tau_max = max(abs(1 ./ Constantes_Tiempo));  % ajustamos el tiempo más lento

% Configuración de simulación
At = tau_min / 10; % Paso de integración pequeño
Tf = 5e-3; % Establecemos que la simulacion dure 5ms
N = round(Tf / At);% Número de pasos
t = 0:At:Tf;% Vector de tiempo

% procedemos a inicializar las variables 
x = [0; 0]; % Estado inicial
u = 12; % la entrada es una tension escalon de 12 V
Vo = zeros(1, N+1); % Tensión en R
I = zeros(1, N+1);  % Corriente
Vc = zeros(1, N+1); % Tensión en el capacitor
U = zeros(1, N+1);  % Entrada aplicada

for k = 1:N+1
    if mod(t(k), 1e-3) < At
        u = -u;
    end %en este if establecimos que la tension de entrada cambie de 12V a -12V cada 1ms

    % registro de señales
    Vo(k) = Mat_C * x;
    I(k) = x(1);
    Vc(k) = x(2);
    U(k) = u; %actualizamos los vectores con los resultados en cada instante de tiempo 

    xp = Mat_A * x + Mat_B * u;
    x = x + At * xp; % usamos euler para evaluar el cambio de las varibles del sistema
end

% Gráficos
figure;

subplot(4,1,1);
plot(t*1e3, U, 'r');
ylabel('u(t) [V]');
title('Entrada escalon alternante de 12 V');
grid on;

subplot(4,1,2);
plot(t*1e3, Vo, 'r');
ylabel('V_R(t) [V]');
title('Tensión en la resistencia');
grid on;

subplot(4,1,3);
plot(t*1e3, I, 'r');
ylabel('i(t) [A]');
title('Corriente en el circuito');
grid on;

subplot(4,1,4);
plot(t*1e3, Vc, 'r');
ylabel('v_C(t) [V]');
xlabel('Tiempo [ms]');
title('Tensión en el capacitor');
grid on;

