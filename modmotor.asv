function [X]=modmotor(t_etapa, xant, accion)
%Valores identificados
La=4.990001e-04; Ra=1.997577e+01; Jm=2.513807e-09;
Km=6.052994e-02; Ki=9.970607e-03; Bm=0;

%Variables para simular el motor usando Euler (se reciben del codigo principal del PID)
ia=xant(1);
theta=xant(2);
wr=xant(3);

Va=accion(1);
TL=accion(2);

t_e=5e-6;%Definicion del tiempo de integracion
pasos=floor(t_etapa/t_e);
%Ciclo for para Forward Euler
for i=1:pasos
ip=-(Ra/La)*ia -(Km/La)*wr +(1/La)*Va;
wrp=(Ki/Jm)*ia -(Bm/Jm)*wr -(1/Jm)*TL;
ia=ia+ip*t_e;
wr=wr+wrp*t_e;
theta=theta+wr*t_e;
end

X=[ia,theta,wr];
end