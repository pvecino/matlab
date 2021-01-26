%% PRACTICAS 1_5 MUESTREO
%% EJERCICIO 1.5.1 MUestreo de señales continuas
clear all
close all
clc
f0_1=1200;
fs_1=8000;
%Obtenemos la se�al 
[s,t] = sinus_c(50,f0_1,45/(2*pi),fs_1,0,7e-3);
%Calculamos la longitud y frecuencia equivalente
len_s_1 = length(s)
siz_s = size(s)
feq_1 = (f0_1/fs_1).*(2*pi)
%Representamos la sinusoide como una se�al continua y como una secuencia
figure
subplot(2,1,1)
plot(t,s)
xlabel('Tiempo [s]         t')
title('Señal continua')
subplot(2,1,2)
n=0:(len_s_1-1);
stem(n,s)
title('Secuencia discreta')
xlabel('Indice             n')

%Obtenemos la se�al 2
f0_2 = 6800;
fs_2 = 8000;
[s,t] = sinus_c(50,f0_2,45/(2*pi),fs_2,0,7e-3);
%Calculamos la longitud  
len_s_2 = length(s);
siz_s = size(s);
feq_2 = (f0_2/fs_2).*(2*pi)
%Representamos la sinusoide como una se�al continua y como una secuencia
figure
subplot(2,1,1)
plot(t,s)
xlabel('Tiempo [s]         t')
title('Señal continua')
subplot(2,1,2)
n=0:(len_s_2-1);
stem(n,s)
title('Secuencia discreta')
xlabel('Indice             n')

%% EJERCICIO 1.5.2
clear all
close all
clc
N_1= 2;
N_2= 3;
N_3 = 6;
% Genero la señal
f0_1=1200;
fs_1=8000;
[x,t] = sinus_c(50,f0_1,45/(2*pi),fs_1,0,7e-3);
len_x = length(x);
n = 0: (len_x-1);

% diezmado
[y1, n1] = diezma(x, N_1);
[y2, n2] = diezma (x, N_2);
[y3, n3] = diezma(x, N_3);


figure
subplot(2, 2, 1)
stem(n, x)
title('x[n]')

subplot(2, 2, 2)
stem(n1, y1)
title('x[n2]')

subplot(2, 2, 3)
stem(n2, y2)
title('x[n3]')

subplot(2, 2, 4)
stem(n3, y3)
title('x[n6]')

%% EJERCICIO 1.5.3

close all
clear all
clc
f0_1=1200;
fs_1=8000;
%Obtenemos la se�al 
[s,t] = sinus_c(50,f0_1,45/(2*pi),fs_1,0,7e-3);
x_original = s;
N = 2;
x_diezmada = diezma(x_original, N);
L = 2;
y1 = expande(x_diezmada, L);
h1 = filtro_interp(L);
x_recuperada_1 = conv(y1,h1);
L = 4;
y2 = expande(x_diezmada, L);
h2 = filtro_interp(L);
x_recuperada_2 = conv(y2, h2);
figure
subplot(3,1,1)
stem(x_original)
title('x_{1original}[n]')
subplot(3,1,2)
stem(x_diezmada)
title('x1[nL]')
subplot(3,1,3)
stem(x_recuperada_1)
title('x_{1recuperada}[n]')
figure
subplot(3,1,1)
stem(x_original)
title('x_{2original}[n]')
subplot(3,1,2)
stem(x_diezmada)
title('x_{2diezmada}[n]')
subplot(3,1,3)
stem(x_recuperada_2)
title('x_{2recuperada}[n]')


%% EJERCICIO 1.5.5
clc
clear all
close all
n = -999:1000;
x1 = 0.5*sin(2*pi*n*1200);
x2 = sin(2*pi*n*1000);
x3 = 2*sin(2*pi*n*1500);
x4 = sinc(n/100);
x5 = sinc(n/200).^2;

x = x1 + x2 + x3 + x4 + x5;
L = 4;
y = diezma(x,L);
m = -999:L:1000;

z1 = interp(y,L);
z2 = interp(y,L*2);
z3 = interp1(m,y,n,'nearest');
z4 = interp1(m,y,n,'linear');
z5 = interp1(m,y,n,'spline');

figure
subplot(3,2,1)
stem(n,x)
title('Señal original')
subplot(3,2,2)
stem(n,z1)
title('Señal interpolada L=4')
subplot(3,2,3)

n2=-1999:2000;
stem(n2,z2)
title('Señal interpolada L=8')
subplot(3,2,4)
stem(n,z3)
title('Interpolacion "nearest"')
subplot(3,2,5)
stem(n,z4)
title('Interpolacion "linear"')
subplot(3,2,6)
stem(n,z5)
title('Interpolacion "spline"')

e1 = x - z1;
e2 = x - z2(1:2:length(z2));
e3 = x - z3;
e3(isnan(e3))=0;
e4 = x - z4;
e4(isnan(e4))=0;
e5 = x - z5;

e1_suma = sum(abs(e1))
e2_suma = sum(abs(e2))
e3_suma = sum(abs(e3))
e4_suma = sum(abs(e4))
e5_suma = sum(abs(e5))