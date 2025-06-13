%Ítem [5] A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se requiere 
%obtener el modelo del sistema considerando como entrada un escalón de 12V, como salida a la 
%velocidad angular, y al torque de carga TL aplicado una perturbación. En el archivo 
%Curvas_Medidas_Motor.xls están las mediciones, en la primer hoja los valores y en la segunda 
%los nombres. Se requiere obtener el modelo dinámico, para establecer las constantes del modelo 
%(1-5) (1-6).  
%Universidad Nacional de Catamarca. Facultad de Tecnologia y Ciencias Aplicadas
%Alumno: Ocampos Gustavo Ismael
%M.U.:983
clear all;
clc;
datos = xlsread('Curvas_Medidas_motor_2024_v');
vector=transpose(datos);
filas=length(vector); %Obtenemos la cantidad de elementos
%Inicializamos los vectores, con la cantidad de filas que posee el archivo 
%transpuesto y solo 
V = zeros(1,filas); %entrada
w = zeros(1,filas);
I_a = zeros(1,filas);
tiempo = zeros(1,filas);
t_l = zeros(1,filas);
% Guardamos cada fila en una variable diferente
for i=1:1:5
    info{i} = vector(i, :);
end

tiempo=info{1};
w=info{2};
I_a=info{3};
V=info{4};
t_l=info{5};
StepAmplitude = 12;

t_inic=0.03506;  %obtenemos el tiempo de inicio y de retardo de la velocidad angular apartir del excel brindado por la catedra
t_retardo=0.03503;
t_diferencia=t_inic-t_retardo;

[val, lugar] =min(abs(t_diferencia+t_retardo-tiempo)); %buscamos dentro del vector tiempo, el valor mas cercano a la operacion t_dif+t_ret.
y_t1=w(lugar);%Una vez obtenido se nos devuelve el indice de la posicion para igresar con el resto de datos
t_t1=tiempo(lugar);
ii=1;
[val, lugar] =min(abs(2*t_diferencia+t_retardo-tiempo));
t_2t1=tiempo(lugar);
y_2t1=w(lugar);

[val, lugar] =min(abs(3*t_diferencia+ t_retardo-tiempo));
t_3t1=tiempo(lugar);
y_3t1=w(lugar);

K=198.2488/StepAmplitude; %tomando los datos del excl , sabemos que la funcion en su maximo tiene 198.2488 rad/s
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

FdT_chen_w_v=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))%Funcion que nos devuelve la funcion de Transferencia

%FIN DE CALCULO DE CHEN
h=log(0.95)/(2*max(roots(conv([T1_ang 1],[T2_ang 1])))); %Determinamos el tiempo de integracion
[w_va_chen, t_w_chen]=lsim(FdT_chen_w_v,V,tiempo,[0 0]); %con uste funcion, yo lo que hago es simular la funcion de transferencia en base a parametros de entrada, para luego guardarlo en vectores para su posterior ploteoo

t_inicio=0,18503 ;
t_end= 0,33502 ;
torque=t_l(18503:33502);
StepAmplitudeTL= mean(torque);
%Tomamos un valor promedio de la entrada, que en este caso queremos evaluar
%la velocidad con respecto a la entrada TL. como dicha entrada posee ruido,
%tomamos un promedio, que lo llamaremos StepAmplitudeTL


%INICIO DE CALCULO DE CHEN PARA W/TL%
t_inic=0.18502;  %Tomando este tiempo es que se logro la mejor aproximacion a la curva original
t_retardo=0.18501; %Tiempo donde inicia la dinamica debido al efecto de la carga TL
t_comparacion=t_inic-t_retardo;

[~, lugar] =min(abs(t_comparacion+t_retardo-tiempo)); %Con esta operacion buscamos dentro del vector tiempo, le que se nos                                    
t__t1=tiempo(lugar);                                    %proporciona como dato, el tiempo mas cercano a la op t_dif+t_ret.
y__t1=w(lugar);                                         %Una vez obtenido se nos devuelve el indice de la posicion para igresar con el resto de datos

ii=1;
[~, lugar] =min(abs(2*t_comparacion+t_retardo-tiempo));
t__2t1=tiempo(lugar);
y__2t1=w(lugar);

[val, lugar] =min(abs(3*t_comparacion+ t_retardo-tiempo));
t__3t1=tiempo(lugar);
y__3t1=w(lugar);


%

K=34.6806/StepAmplitudeTL; % 34.6806 es la diferencia entre el valor de tension maximos y el valor a cuando se produce el torque
k1=(1/StepAmplitudeTL)*(y__t1/K)-1;%Afecto el valor del Escalon
k2=(1/StepAmplitudeTL)*(y__2t1/K)-1;
k3=(1/StepAmplitudeTL)*(y__3t1/K)-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;

alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));

 beta2=(k1+alfa2)/(alfa1-alfa2);

T1_ang=-t_comparacion/log(alfa1);
T2_ang=-t_comparacion/log(alfa2);
T3_ang=beta2*(T1_ang-T2_ang)+T1_ang;

T1(ii)=T1_ang;
T2(ii)=T2_ang;
T3(ii)=T3_ang;
T3_ang=0.05*real(sum(T3/length(T3))); %Se considero que el valor de 0.05 permite una comparacion optima del ruido provocado por la carga
T2_ang=real(sum(T2/length(T2)));      % EN LA GRAFICA SIMULADA en comparacion con la curva real de la carga TL
T1_ang=real(sum(T1/length(T1)));

FdT_chen_w_TL=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))%Funcion que nos devuelve la funcion de Transferencia
%Esta representa la FdT de W(s)/Tension

%FIN DE CALCULO DE CHEN PARA W/TL
[w_tl_chen, t_chen, entrada]=lsim(FdT_chen_w_TL,t_l,tiempo,[0 0]); %simulo la funcion de transferencia en base a parametros de entrada, para luego guardarlo en vectores para su posterior ploteo
%
%Seccion de calculo de la respuesta completa de la Velocidid ángular W para
%las entradas de torque y tension.

simulada=w_va_chen-w_tl_chen;

%INICIO DE CALCULO DE CHEN MODIFICADO PARA I/V

%como la curva de la corriente no presenta un valor en estado estacionario, por lo que
%no podemos aplicar de manera indicada el metodo de Chen, ademas de que su funcion de transferencia
%posee un polo en el origen. Esto nos lleva a introducir una variacion de
%dicho metodo que nos fue brindada por la catedra
t_retardo=0.03502; %Tiempo donde inicia la dinamica de la curva de la velocidad angular con respecto a Va


 %Recorte de los datos
Y_aux = I_a(31:16860);
t_aux = tiempo(31:16860);
u_aux = V(31:16860);

y1i = max(Y_aux);
uMAX = max(V);

% t1i, y1i
[~, lugar] = min(abs(Y_aux - y1i));
t1i = t_aux(lugar);
y1i = y1i / uMAX;

% t2i, y2i
[~, lugar] = min(abs((t1i - 0.03502)*2 + 0.03502 - tiempo));
t2i = tiempo(lugar);
y2i = I_a(lugar) / uMAX;

% t3i, y3i
[~, lugar] = min(abs((t1i - 0.03502)*3 + 0.03502 - tiempo));
t3i = tiempo(lugar);
y3i = I_a(lugar) / uMAX;

ii=1;

% Cálculo de alfa1 (según fórmula equivalente en MATLAB)
ret = 0.03502;
den = 4*y1i*y3i - 3*y2i^2;
alfa1_expr = (-4*y1i*y3i * sqrt(1/den) + 3*y2i^2 * sqrt(1/den) + y2i) / (2*y1i);
alfa1 = abs(alfa1_expr);

% alfa2
alfa2 = y2i / y1i - alfa1;

% Constantes de tiempo T1 y T2
T1 = -(t1i - ret) / log(alfa1);
T2 = -(t1i - ret) / log(alfa2);

% Parámetros del sistema
beta3 = y1i / (alfa1 - alfa2);
beta4 = y2i / (alfa1^2 - alfa2^2);
beta5 = y3i / (alfa1^3 - alfa2^3); 

K_2_C = y1i * (T1 - T2) / (alfa1 - alfa2); % y1i ya normalizado

% Función de transferencia
s = tf('s');
FdT_I_Va = K_2_C * (s) / ((T1*s + 1)*(T2*s + 1))

% Simulación del sistema
[y__1, t__1, ~] = lsim(FdT_I_Va, V, tiempo, [0; 0]);


%FIN DE CALCULO DE CHEN MODIFICADO PARA I/va

%EXTRACCION DE PARAMETROS:

%compararemos los terminos de las funciones transferencia obtenidos con las
%teoricas para encontrar el valor de los parametros.
%Los terminos a obtener seran R, L, Km, J y Ki

% Su Funcion es (S/L)/(S^2 + SR/L + KmKi/JL)


%De la FdT de la CORRIENTE
[num_IA, den_IA] = tfdata(FdT_I_Va, 'v');

NORMALIZADO_IA=den_IA(1);  %Factor de normalizacion

num_n_IA=num_IA/NORMALIZADO_IA;
den_n_IA=den_IA/NORMALIZADO_IA;

FdT_I_Va_n=tf(num_n_IA,den_n_IA)

Ian=num_n_IA(2);
Iad1=den_n_IA(1);
Iad2=den_n_IA(2);
Iad3=den_n_IA(3);


%Despejand de la FdT podemos llegar a que:

L=1/Ian
R=Iad2*L

%De la FdT de la celocidad angular con respecto a la tension W/V:
% Su Funcion es (Ki/J)/(S^2 + SR/L + KmKi/JL)

[num_WVa, den_WVa] = tfdata(FdT_chen_w_v, 'v');

NORMALIZADO_WVa=den_WVa(1); %Factor de normalizacion

% *Normalizamos numerador y denominador*
num_n_WVa=num_WVa/NORMALIZADO_WVa;
den_n_WVa=den_WVa/NORMALIZADO_WVa;

FdT_WVa_n=tf(num_n_WVa,den_n_WVa)

Wva_n2=num_n_WVa(3);
Wvad1=den_n_WVa(1);
Wvad2=den_n_WVa(2);
Wvad3=den_n_WVa(3);

Km=(Wvad3/Wva_n2)*L

%De la FdT de la velocidad angular con respecto al torque w/tl:
% Su Funcion es (S/J + R/LJ)/(S^2 + SR/L + KmKi/JL)

[num_WTL, den_WTL] = tfdata(FdT_chen_w_TL, 'v');

NORMALIZADO_WTL=den_WTL(1);  %Factor de normalizacion

% *Normalizamos numerador y denominador*
num_n_WTL=num_WTL/NORMALIZADO_WTL;
den_n_WTL=den_WTL/NORMALIZADO_WTL;

FdT_WTL_n=tf(num_n_WTL,den_n_WTL)

W_WTL_n2=num_n_WTL(2);
W_WTL_n3=num_n_WTL(3);
W_WTL_d1=den_n_WTL(1);
W_WTL_d2=den_n_WTL(2);
W_WTL_d3=den_n_WTL(3);
%con el nimerador y denominador normalizados realizamos las ecuaciones con
%las que obtendremos J y Ki
J= (1/W_WTL_n2)
Ki=((W_WTL_d3*J*L)/Km)
%una vez obtenidos los datos , concluimos con las graficas de los valores
%obtenidos y comparandolas con las graficas de referencia 

 figure;

subplot(4,1,3);
plot(tiempo, V, 'b', 'LineWidth', 1.5);
ylabel('Tensión [V]');
title('Tensión / Tiempo');
grid on;

subplot(4,1,2);
plot(tiempo, I_a , 'r', 'LineWidth', 1.5);
hold on;
plot(t__1, y__1, 'b');
plot([t1i t2i t3i], [y1i y2i y3i]*uMAX, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
legend('modelo medido','modelo ajustado','puntos medidos');
ylabel('Corriente [A]');
title('Corriente de Armadura / Tiempo');
grid on;

subplot(4,1,1);
plot(tiempo,w,'r'),hold on, plot(t_w_chen,simulada,'b');
hold on;

plot(t_t1, y_t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
plot(t_2t1, y_2t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot(t_3t1, y_3t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
plot(t__t1, y__t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
plot(t__2t1, y__2t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot(t__3t1, y__3t1, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'c')
legend('Velocidad medida', 'Aprox. Chen');
ylabel('Velocidad [rad/s]');
title('Velocidad Angular / Tiempo');
grid on;

subplot(4,1,4);
plot(tiempo, t_l, 'k', 'LineWidth', 1.5);
ylabel('T_L [Nm]');
xlabel('Tiempo [s]');
title('Torque de Carga / Tiempo');
grid on;