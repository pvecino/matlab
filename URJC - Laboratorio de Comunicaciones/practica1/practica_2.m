%% Practica 1-2: Representacion de señales
%% 1. Secuencias basicas
%% 1.A Secuencias Reales
%% Ejercicio 1.2.1
clear all
close all
clc

n1 = 0:20;
x1 = (0.9.*delta(n1-5));
figure
subplot(3,1, 1)
stem(x1)
title('0.9.*delta(n1-5)');

n2 = -10:10;
x2 = (4.5.*delta(n2+7));
subplot(3,1, 2)
stem(x2)
title('4.5.*delta(n2+7)');

n3 = 0:20;
x3 = (4.5.*delta(n3+17));
subplot(3,1, 3)
stem(x3)
title('4.5.*delta(n3+17)');

%% Ejercicio 1.2.2
clear all
close all
clc
%METODO BUCLE FOR
tic
n= 1:500;
P = 5;
s = zeros(1, length(n));
for k=1:length(n)
    if mod(k-1,P) == 0
        s(k) = 1;
    end
end
toc

figure
subplot(2,1,1)
stem(n, s)

% METODO MOD
tic
n2=1:50;
s2=zeros(1,length(n2));
s2(mod(n2,P)==0) = 1;
toc
subplot(2,1,2)
stem(n2, s2)

%% Ejercicio 1.2.3
clear all
close all
clc
[x, n] = sinus(5, 10, 3, 0, 3);
plot(n, x)
%% Ejercicio 1.2.4
clear all
close all
clc
b = 0.9
n0 = 0
L = 21;
% A)
[x, n] = expon(b, n0, L);
figure
subplot(2, 1, 1)
plot(n, x)
title('x[n] = (0.9)^n');
xlabel('n');
% B)
sum1 = sum(x)
sum2 = ((1-(b.^L))./(1-b))
% C)

n = 0:20;
x2 = delta(n);

a = [1 -0.9];
b2 =[1];

y = filter(b2, a, x2);

subplot(2, 1, 2)
t = 1: length(x);
plot(t, y)
title('con filter');
xlabel('n');
%% 1.A Secuencias Reales
%% Ejercicio 1.2.5
clear all
close all
clc

n1 = 0:20;
x1 = (7*i*exp((-1*i*pi*n1)/7)+i*exp((i*pi*n1)/7))/2;

figure(2)
subplot(2,2,1)
stem(n1,real(x1))
title('Parte real')
subplot(2,2,2)
stem(n1,imag(x1))
title('Parte imaginaria')

subplot(2,2,3)
stem(n1,abs(x1))
title('Modulo')
subplot(2,2,4)
stem(n1,angle(x1))
title('Fase')

n2 = 0:50;
x2 = ((1.1).^n2).*cos(((pi*n2)./11)+(pi./4));

figure(3)
subplot(2,2,1)
stem(n2,real(x2))
title('Parte real')
subplot(2,2,2)
stem(n2,imag(x2))
title('Parte imaginaria')
subplot(2,2,3)
stem(n2,abs(x2))
title('Modulo')
subplot(2,2,4)
stem(n2,angle(x2))
title('Fase')

%% 2.Representacion de secuencias en el dominio de la frecuencia
%% Ejercicio 1.2.6
a = [1 1 1 1 1];
N = 5;
[x,n]=dsf_sintesis(a,N);
a2 = dsf_analisis_fft(x,N);

figure
subplot(2,1,1);
stem(a);
axis([1,5,0,1.5]);
title('Coeficientes iniciales');
subplot(2,1,2);
stem(abs(a2));
title('Coeficientes finales');

%% Ejercicio 1.2.7

%Generamos la se�al
x = [ones(1,8) zeros(1, 24)];
P_vec = [4,8,16,32];
%Calculamos los coeficientes
a = dsf_analisis_fft(x,32);
figure;
for count = 1:length(P_vec)
    P = P_vec(count);
    %A priori inicializamos todos los coeficientes a cero
    a_P = zeros(size(a));
    %Nos quedamos con los coeficientes que van desde 0 hasta P/2
    a_P(1:(P/2+1)) = a(1:(P/2+1));
    %Y nos quedamos con los coeficientes que van desde -P/2 hasta -1 (o
    %equivalentemente desde N-P/2 hasta N-1
    a_P(end-(P/2-1):end) = a(end-(P/2-1):end);
    %Calculamos la se�al con esos coeficientes
    x_P = dsf_sintesis(a_P,32);
    
    %Dibujamos la se�al
    subplot(4,1,count)
    titulo=sprintf('P = %d', P);
    title(titulo);
    hold on
    stem(abs(x))
    grid on
    %Dibujamos el valor absoluto para evitar problemas con la parte
    %imaginaria que aparece (del orden de 1e-10*i) debido al error 
    %numerico 
    stem(abs(x_P),'r')
    legend('x', 'x_P');
end

