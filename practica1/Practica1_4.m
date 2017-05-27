%% PRACTICA 1.4 FILTRADO
%% EJERCICIO 1.4.1
clc
clear all
close all
 %Generamos los �ndices de la se�ales x h
n_x = 0:14;
n_h = 0:2;
%Generamos los valores de la se�ales x h
x = [ones(1,5), zeros(1,5), ones(1,5)];
h = [1 3 1];%deltas solo valores
%Generamos los �ndices de la se�al y
ini_y = n_x(1) + n_h(1);
fin_y = n_x(end) + n_h(end);
n_y = ini_y:fin_y;
%Generamos los valores de la se�al
y = conv(x,h);
%Dibujamos las se�ales
figure;
subplot(2,2,1)
stem(n_x,x)
title('x[n]')
subplot(2,2,2)
stem(n_h,h)
title('h[n]')
subplot(2,2,[3,4])
stem(n_y,y);
title('y[n]');
%% EJERCICIO 1.4.2
clear all
clc
close all
%Generamos los indices de la se�al x
n_x = 0:14;
%Generamos los valores de la se�al x 
x = [ones(1,5) zeros(1,5) ones(1,5)];
%Generamos los coeficientes del filtro A= Y//B=X
b = [1 3 1];
a = 1;
%Generamos los valores de la se�al y los �ndices
y = filter(b,a,x);
n_y = n_x;
%Dibujamos las se�ales
figure; 
stem(n_y,y);
title('y[n]');

%% EJERCICIO 1.4.3
clear all
close all
clc
w_0 = 2*pi*5/32;
n = 0:50;
%Calculamos el coeficiente para que la respuesta en frecuencia
% se %anule en w_0
a_coef = -(1+exp(-2*j*w_0))/exp(-j*w_0);
%Generamos la sinusoide
x = cos(w_0*n);
%Coeficientes de la ecuacion en diferencias
b = [1 a_coef 1];
%Calculamos y como la salida del filtro que implementa la
% ecuacion en diferencias
y = filter(b, 1, x);
%Representamos las se�ales x e y para comprobar la salida
figure
plot(n,x)
hold on
stem(n,y,'r')
legend('x[n]','y[n]')

%% EJERCICIO 1.4.4
clc
clear all
close all
%Generamos la se�al, para ello tenemos en cuenta que:
N = 500;
fs =  10000;%Frecuencia de muestro 
n = (-N/2):(N/2-1); %Eje de tiempos 
x_ideal = 10*sinc(n./10); % (n/fs)/0.001
ruido = randn(1,500);
x = x_ideal + ruido;
%Calculamos y mostramos las potencias (en el intervalo)
% para hacernos una mejor idea del nivel de ruido
P_x_ideal = mean(abs(x_ideal).^2);
P_ruido = mean(abs(ruido).^2);
P_x = mean(abs(x).^2);
%Dise�amos los distintos filtros paso bajo y los aplicamos a la
% se�al ruidosa
f_corte = 1500;
Kf = 10; %Orden del filtro si aumento 10 o lo reduzco la sinc se me deformA
Ff= f_corte/(fs/2);
[b,a] = butter(Kf,Ff);
y_butter = filter(b, a, x);
[b,a] = cheby1(Kf,1, Ff); %R es la atenuacion maxima 
y_cheby1 = filter(b, a, x);
[b,a] = cheby2(Kf,20, Ff);% si R aumenta la sinc se deforma
y_cheby2 = filter(b, a, x);
[b,a] = fir1(Kf, Ff);
y_fir1 = filter(b, a, x);
%Calculamos y mostramos la potencia de las se�ales filtradas
P_y_butter = mean(abs(y_butter).^2);
P_y_cheby1 = mean(abs(y_cheby1).^2);
P_y_cheby2 = mean(abs(y_cheby2).^2);
P_y_fir1 = mean(abs(y_fir1).^2);
%Dibujamos las distintas se�ales
figure
subplot(3,2,1)
stem(n,x_ideal)
title('x ideal')
axis([-100 100 min(x_ideal) max(x_ideal)])
subplot(3,2,2)
stem(n,x)
title('x con ruido')
axis([-100 100 min(x) max(x)])
subplot(3,2,3)
stem(n,y_butter)
title('x filtrada butter')
axis([-100 100 min(x) max(x)])
subplot(3,2,4)
stem(n,y_cheby1)
title('x filtrada cheby1')
axis([-100 100 min(x) max(x)])
subplot(3,2,5)
stem(n,y_cheby2)
title('x filtrada cheby2')
axis([-100 100 min(x) max(x)])
subplot(3,2,6)
stem(n,y_fir1)
title('x filtrada fir1')
axis([-100 100 min(x) max(x)])

%% EJERCICIO 1.4.6
clear all
close all
clc

% Generamos la se�al, para ello tenemos en cuenta que:
N = 500;
fs = 10000; % Frecuencia de muestrEo 
n = (-N/2):(N/2-1); % Eje de tiempos 
x1 = 2*sinc((n/fs)./0.002); % BAJAS
x2 = sinc((n/fs)./0.002).*cos(2*pi*3000*n./fs); %ALTAS
x = x1+x2;

figure
subplot(2,2,1)
stem(n,x1)
title('x1[n]')
subplot(2,2,2)
stem(n,x2)
title('x2[n]')
subplot(2,2,[3,4])
stem(n,x)
title('x[n] = x1[n] + x2[n]')

load('lp_filter.mat');
load('hp_filter.mat');
y1 = filter(lp_filter, x);
y2 = filter(hp_filter, x);

figure
subplot(2,2,1)
stem(n,x1)
title('x1[n]')
subplot(2,2,2)
stem(n,y1)
title('y1[n]')
subplot(2,2,3)
stem(n,x2)
title('x2[n]')
subplot(2,2,4)
stem(n,y2)
title('y2[n]')