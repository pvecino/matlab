%% LSC Practica 1
%% Ejercicio 1.1.1
clc
clear all
close all
% x con longitud 1000
len = 1000;
xaux = rand(1,len);
% para convertir, todos los valores de rand, se pone la condicion que el
% 60% sean 0, y el 40% sean 1.
x = [(xaux > 0.6)];


%x e y tienen la misma longitud
ylen= length(x);
yaux = rand(1, ylen);
% x e y tienen el 90% de ser iguales
yaux = [yaux > 0.9];
%por tanto el 10% es de ser distintas, con xor
y = xor(x, yaux);

figure
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
hist(s(1), x)
hist(s(2), y)

res = [x~= y];
dif = sum(res); 
%% Ejercicio 1.1.2
clc
clear all
close all
% x con longitud 1000
len = 1000;
xaux = rand(1,len);
% para convertir, todos los valores de rand, se pone la condicion que el
% 60% sean 0, y el 40% sean 1.
x = [(xaux > 0.6)];


%x e y tienen la misma longitud
ylen= length(x);
yaux = rand(1, ylen);
% x e y tienen el 90% de ser iguales
yaux = [yaux > 0.995];
%por tanto el 10% es de ser distintas, con xor
y = xor(x, yaux);

figure
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
hist(s(1), x)
hist(s(2), y)

res = [x~= y];
dif = sum(res); 

%% Ejercicio 1.1.3
clc
clear all
close all
%inicializo/creo la matriz 
A = zeros(10, 20);


for n=1:10,
    for m=1:20,
        A(n,m) = 5 * n + (-1).^m;
    end
end

% Suma con el metodo 1
suma_for = 0;
pos = m*n;
for k=1:pos
    suma_for = suma_for + A(k);
end
  
% Suma con el metodo 2
% sum(A,2) devuelve una matriz de una columna con la suma de cada fila.
suma2 = sum(sum(A,2)); 

% Suma con el metodo
% sum(A,1) devuelve una matriz de una fila con la suma de cada columna.
suma3 = sum(sum(A,1)); 
%% Ejercicio 1.1.4
clc
clear all
close all
% ind =find(x==numero_buscado)
% z(ind)
%TODO DEPENDE DE LO QUE PASOS QUE DES
h = 0.001;% el paso debe ser pequeño, 
%para que la aproximacion sea buena
x = 1:h:1000;%debe ser positivo
y = log10(1+x);
z = diff(y/h);

figure
subplot(1,2,1)
plot(y)
title('log10(1+x)')
subplot(1,2,2)
plot(z)

a = [0:h:4];
f = (cos(a./2)).^2;
int = cumsum(f)*h; 
% lo que va valiendo la integral de 0 a 4 en pasos 
%de 1, si variamos el paso de 1 a 0,01
int2 = sum(f)*h 
% el sumatorio final, la cifra final del sumatorio, se debe 
%multipicar por el paso que se de.
%% Ejercicio 1.1.5
clc
close all
clear all

x1 = -1:0.04:1;
y1 =  (3.*x1) - 1 ;

x2 = 1:0.02:3;
y2 = ((x2.^3) + (-2*(x2.^2)) + 3);

plot(x1, y1)
hold on
plot(x2, y2)
title('Funcion y');
xlabel('x');
ylabel('y');

%% Ejercicio 1.1.6
clc
clear all
close all
x = 0:0.02:8;

y = 2 .* cos(x ./ 2);

figure; 
subplot(2,1,1); 
plot(x, y, 'r:', 2, 2, 'r+')
title('y=2cos(x/2)')


z = 5 - (x-4).^2;
subplot(2,1,2); 

[m, b] = max(z);
x_max= x(b);

plot(x, z, 'k')
hold on;
plot(x_max, m, 'go') 
title('Y=5-(x-4)^2')

%%
clear all
clc
x = [1 2 18 25 69 3 7 8 9 10];

y = mifuncion(x, 5)