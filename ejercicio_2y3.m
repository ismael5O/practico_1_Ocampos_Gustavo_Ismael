%Ítem [3] Una vez determinados los parámetros R, L y C, emplear la serie de corriente desde 
%0.05seg en adelante para validar el resultado superponiendo las gráficas.
%Universidad Nacional de Catamarca. Facultad de Tecnologia y Ciencias Aplicadas
%Alumno: Ocampos Gustavo Ismael
%M.U.:983
clear all;
clc;
datos = xlsread('Curvas_Medidas_RLC_2024');
vector=transpose(datos);
filas=length(vector); %Obtenemos la cantidad de elementos

%Inicializamos los vectores, con la cantidad de filas que posee el archivo
%transpuesto, y solo 
u = zeros(1,filas); %entrada
I = zeros(1,filas);
V_resistencia = zeros(1,filas);
V_c = zeros(1,filas);
tiempo=zeros(1,filas);
% Guardar cada fila en una variable diferente
for i=1:1:4
    info{i} = vector(i, :);
end

tiempo=info{1};
I=info{2};
V_c=info{3};
u=info{4};

%Una vez que tenemos los datos cargados, debemos muestrear nuestra
%graficas, con el fin de obtener una funcion de transferencia del sistema
%y con ello luego obtener los parametros

h=0.0001; %viendo las tablas del tiempo, vemos que el paso es cada 0.0001 sg
StepAmplitude = 12;

%INICIO DE CALCULO DE CHEN%
 

t_inic=0.0105;  
t_retardo=0.0101; %Tiempo donde inicia la dinamica de la curva del Capacitor
t_diferencia=t_inic-t_retardo;

[val, lugar] =min(abs(t_diferencia+t_retardo-tiempo)); %Con esta operacion buscamos dentro del vector tiempo, le que se nos                                    
t_t1=tiempo(lugar);                                   %proporciona como dato, el tiempo mas cercano a la op t_dif+t_ret.
y_t1=V_c(lugar);                                    %Una vez obtenido se nos devuelve el indice de la posicion para igresar con el resto de datos

ii=1;
[val, lugar] =min(abs(2*t_diferencia+t_retardo-tiempo));
t_2t1=tiempo(lugar);
y_2t1=V_c(lugar);

[val, lugar] =min(abs(3*t_diferencia+ t_retardo-tiempo));
t_3t1=tiempo(lugar);
y_3t1=V_c(lugar);

K=12/StepAmplitude; %La funcion al final de su trayectoria tiene una salida de 12V, ya qe es la tension a la que se carga el cap
k1=(1/StepAmplitude)*(y_t1/K)-1;%Afecto el valor del Escalon
k2=(1/StepAmplitude)*(y_2t1/K)-1;
k3=(1/StepAmplitude)*(y_3t1/K)-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;

alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));

 beta=(k1+alfa2)/(alfa1-alfa2);

%beta=(2*k1^3+3*k1*k2+k3-sqrt(be))/(sqrt(be)); %Beta para dos polos distintos y un cero



T1_ang=-t_diferencia/log(alfa1);
T2_ang=-t_diferencia/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;

T1(ii)=T1_ang;
T2(ii)=T2_ang;
T3(ii)=T3_ang;

T3_ang=real(sum(T3/length(T3)));
T2_ang=real(sum(T2/length(T2)));
T1_ang=real(sum(T1/length(T1)));

FdT_chen=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))%Funcion que nos devuelve la funcion de Transferencia

%FIN DE CALCULO DE CHEN

h=log(0.95)/(2*max(roots(conv([T1_ang 1],[T2_ang 1])))); %Determinamos el tiempo de integracion
[y_chen, t_chen, entrada]=lsim(FdT_chen,u,tiempo,[0 0]); %simulamos la funcion de transferencia en base a parametros de entrada.
%posteriormente lo guardamos en vectores para su posterior ploteo

%y_chen la tension en el capacitor.
%Como Vc=1/c*integral(i(t)) por lo que si derivamos y despejamos encontraremos la
%corriente del circuito

Corriente=diff(y_chen)/h;

%Ahora debemos encontrar C para obtener I. 
%al plotear la corriente, esta tendra un valor muy grande; Por otro lado,
%en la hoja de datos proporcionada, tenemos la I medida, por lo
%que lo que separa ambas graficas sera nada mas que un factor de escala,
%osea C
%lo que significa que podremos obtener C al realizar un cociente entre la I
%medida de los datos y la I obtenida de corriente.

C=max(I)/max(Corriente);

I_id=C*Corriente;
tiempo_corr = tiempo(1:end-1);
I_identificada=transpose(I_id);

%La ecuacion caracteristica de un RLC es CL*S^2+S*RC+1, con lo que dado que
%tenemos C, podemos empezar a comparar terminos y asi sacar el resto de
%componentes: CL=1.02e-6, con lo que L sera:

L=1.02e-6/C;

%A su vez aplicamos lo mismo para R, en el cual RC=0.002699

R=0.002699/C;

%GRAFICAS

figure;

subplot(3,1,1);
plot(tiempo, I, 'b'); hold on;
plot(tiempo_corr, I_identificada, 'r');
legend('Curva real', 'Curva aproximada con Chen');
title('Gráficas de Corriente');

subplot(3,1,2);
plot(tiempo, V_c); hold on;
plot(t_chen, y_chen, 'r');
plot(t_t1, y_t1, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'y');
plot(t_2t1, y_2t1, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'y');
plot(t_3t1, y_3t1, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'y');
legend('Curva real','Curva con Chen','t1','2t1','3t1');
title('Tensión en el capacitor');
ylim([-15 15]);
line([tiempo(1), tiempo(end)], [12, 12], 'Color', 'g', 'LineStyle', '--');
line([tiempo(1), tiempo(end)], [-12, -12], 'Color', 'g', 'LineStyle', '--');

subplot(3,1,3);
plot(tiempo, u);
ylim([-15 15]);
line([tiempo(1), tiempo(end)], [12, 12], 'Color', 'r', 'LineStyle', '--');
line([tiempo(1), tiempo(end)], [-12, -12], 'Color', 'r', 'LineStyle', '--');
title('Tensión de entrada u');
