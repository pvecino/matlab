%% Practica 1.3 Dominio transformado y enventanado
%% EJERCICIO 1.3.1
clc
clear all
close all
p = ones(1,7);
%Generamos la DFT original fft(señal, longitud de la dft)
P_dft = fft(p,16);
k = 0:15;
figure
subplot(2,1,1)
stem(k,abs(P_dft))
xlabel('k')
ylabel('M�dulo')
subplot(2,1,2)
stem(k, angle(P_dft))
xlabel('k')
ylabel('Fase')

%% EJERCICIO 1.3.2
clc
close all
clear all
L = 21;
%% A) N = L
N = L;
F = 1/N
n = 0:(L-1);
s = cos(2 * pi / N * n);
%Generamos la DFT original
S_dft = fft(s, L);
k = 0:(L-1);
figure
subplot(3,1,1)
stem(n,s)
xlabel('n')
ylabel('Señal')
subplot(3,1,2)
stem(k,abs(S_dft))
xlabel('k')
ylabel('Modulo DFT')
subplot(3,1,3)
stem(k, angle(S_dft))
xlabel('k')
ylabel('Fase DFT')
%% B) s = seno y N=L
N = L;
n = 0:(L-1);
s = sin(2 * pi / N * n);
%Generamos la DFT original
S_dft = fft(s, L);
k = 0:(L-1);
figure
subplot(3,1,1)
stem(n,s)
xlabel('n')
ylabel('Señal')
subplot(3,1,2)
stem(k, abs(S_dft))
xlabel('k')
ylabel('Modulo DFT')
subplot(3,1,3)
stem(k, angle(S_dft))
xlabel('k')
ylabel('Fase DFT')
%% c) N=L/3
N = L/3;
F2 = 1/N
n = 0:(L-1);
s = cos(2 * pi / N * n);
%Generamos la DFT original
S_dft = fft(s, L);
k = 0:(L-1);
figure
subplot(3,1,1)
stem(n,s)
xlabel('n')
ylabel('Señal')
subplot(3,1,2)
stem(k,abs(S_dft) )
xlabel('k')
ylabel('Modulo DFT')
subplot(3,1,3)
stem(k,angle(S_dft) )
xlabel('k')
ylabel('Fase DFT')
%% D) F = 3,1
N = L*3.1;
n = 0:(L-1);
s = cos(2 * pi / N * n);
%Generamos la DFT original
S_dft = fft(s, L);
k = 0:(L-1);
figure
subplot(3,1,1)
stem(n,s)
xlabel('n')
ylabel('Señal')
subplot(3,1,2)
stem(k, abs(S_dft))
xlabel('k')
ylabel('Modulo DFT')
subplot(3,1,3)
stem(k, angle(S_dft))
xlabel('k')
ylabel('Fase DFT')

%% EJERCICIO 1.3.3
clc
clear all
close all
N = 16;
n = 0:(N-1);
x = 0.9.^n;
% a) fft
X_dft = fft(x, N);
k = 0:(N-1);
figure
subplot(2,1,1)
stem(k,abs(X_dft))
xlabel('k')
title('Módulo')

% b) teorica
%linspace(inicio, fin, valores que tiene +1) si no se pone N+1 la
% diferencia no te sale conts.
Omega=linspace(0,2*pi,N+1);
% omega es la longitud de la dft, por tanto , (1:N)
Omega=Omega(1:N);
Xomega=1./(1-0.9*exp(-1j*Omega));

subplot(2,1,2)
stem(Omega, abs(Xomega));
title('|X(\Omega)|');
xlabel('\Omega');
% c) comparacion

figure (3)
stem(Omega,abs(X_dft), 'k')
xlabel('k')
ylabel('Módulo')
hold on
stem(Omega, abs(Xomega), 'b');
title('|X(\Omega)|');
xlabel('\Omega');

legend('|X(k)|', '|X(\Omega)|');

diferencia = abs(X_dft) ./ abs(Xomega)

%% EJERCICIO 1.3.4
clc
close all
clear all
%Generamos la se�al original muestreada
Ts = 1/10;
Fs = 1/Ts;
t = 0:Ts:2;% pongo Ts porque estoy muestreando
x = 1- (abs(t-1));
% A) Generamos la DFT original
L = 50;
X_dft = fft(x,L);
k = 0:(L-1);
% B) Generamos la DF de la secuencia, en este caso podemos dibujar la TF en
%cualquier intervalo de longitud 2*pi
X_TF_sec = fft(x, L);
Omega = (-pi:2*pi/(L-1):pi);
% C) Generamos la DF de la se�al original, en este caso la TF debe 
%representarse con un espectro sim�trico (debe utilizarse la funci�n fftshift) 
X_TF_sen = fftshift(fft(x, L));
omega = (-Fs/2: Fs/(L-1): Fs/2);
% Comparacion
figure
subplot(3,1,1)
stem (k,abs(X_dft)) %discreta
xlabel('k')
ylabel('DFT: X[k]')
subplot(3,1,2)
stem (Omega,abs(X_TF_sec)) %discreta con periodo 2pi
xlabel('\Omega')
ylabel('TF: X(\Omega)')
subplot(3,1,3)
 plot (omega,abs(X_TF_sen)) %continua 
xlabel('\omega')
ylabel('TF: X(\omega)')

%% EJERCICIO 1.3.5
clc
clear all
close all

L = 10;
N = 10;
x = [ones(1,L) zeros(1,N-L)];
k = 0:N-1;
dft_x = fft(x,N);
dft_y = (dft_x).^2;
idft_y = ifft(dft_y,N);

figure
subplot(3,1,1)
stem(k,idft_y)
title('N=10');
ylabel('y[n]');
xlabel('n');

N = 20;
x = [ones(1,L) zeros(1,N-L)];
k = 0:N-1;
dft_x = fft(x,N);
dft_y = (dft_x).^2;
idft_y = ifft(dft_y,N); 
subplot(3,1,2)
stem(k,idft_y)
title('N=20');
ylabel('y[n]');
xlabel('n');

N = 40;
x = [ones(1,L) zeros(1,N-L)];
k = 0:N-1;
dft_x = fft(x,N);
dft_y = (dft_x).^2;
idft_y = ifft(dft_y,N);
subplot(3,1,3)
stem(k,idft_y)
title('N=40');
ylabel('y[n]');
xlabel('n');

%% EJERCICIO 1.3.6
clc
clear all
close all

%Generamos el impulso (delta)
n = 0:9;
d = [1, zeros(1,9)];
%Dise�amos el sistema de ecuaciones en diferencias
%A = y(n)
%B = x(n)
A = [1 1/6 -1/6];
B = [1 3 -1];
%Generamos la salida cuando la entrada es una delta (respuesta al impulso)
h = filter(B,A,d);
%Aunque no nos lo piden la dibujamos
figure
stem(n,h)
xlabel('n')
title('Respuesta al impulso')
%Generamos el diagrama de polos y ceros
figure
zplane(B,A)

%% EJERCICIO 1.3.7
clc
clear all
close all
N = 512;
%Generamos las ventanas en el dominio del tiempo
w_rec_21 = [ones(1,21) zeros(1,N-21)];
w_rec_8 = [ones(1,8) zeros(1,N-8)];
w_rec_64 = [ones(1,64) zeros(1,N-64)];
n_i = 0:(N-1);
%Calculamos DTFs. Queremos 512 muestras en el intervalo
% [-pi,pi), as� que tomamos 513 en [-pi,pi] y tiramos la �ltima 
Omega = linspace(-pi,pi,N+1);
Omega = Omega(1:N);
W_rec_21 = fftshift(fft(w_rec_21,N));
W_rec_8 = fftshift(fft(w_rec_8,N));
W_rec_64 = fftshift(fft(w_rec_64,N));
%Normalizamos los espectros y los pasamos a dBs
W_rec_21_dB = 10*log10(W_rec_21./max(W_rec_21));
W_rec_8_dB = 10*log10(W_rec_8./max(W_rec_8));
W_rec_64_dB = 10*log10(W_rec_64./max(W_rec_64)); 
%Generamos las figuras correspondientes
figure
subplot(2,3,1)
stem(n_i,w_rec_21)
ylabel('w_2_1[n]')
axis([0 128 -.5 1.5])
subplot(2,3,2)
stem(n_i(1:256),w_rec_8(1:256))
ylabel('w_8[n]')
axis([0 128 -.5 1.5])
subplot(2,3,3)
stem(n_i(1:256),w_rec_64(1:256))
ylabel('w_6_4[n]')
axis([0 128 -.5 1.5])
subplot(2,3,4)
plot(Omega,W_rec_21_dB)
ylabel('W_2_1(\Omega)')
axis([-pi pi -30 1])
subplot(2,3,5)
plot(Omega,W_rec_8_dB)
ylabel('W_8(\Omega)')
axis([-pi pi -30 1])
subplot(2,3,6)
plot(Omega,W_rec_64_dB)
ylabel('W_6_4(\Omega)')
axis([-pi pi -30 1])
%% EJERCIO 1.3.8
close all
clear all
clc
N = 5;
% APARTADO A)
n = 0:10;
x1 = bartlett(length(n));
figure
stem(n, x1)
title('Funcion bartlett')

% APARTADO B)
clc

N1=31;
N2=61;

n1=0:30;
n2=0:60;

w1=triang(N1);
w2=triang(N2);

W1=fftshift(fft(w1,N*N1));
W2=fftshift(fft(w2,N*N2))

W1_db=10*log10(W1./max(W1));
W2_db=10*log10(W2./max(W2));

k1=-(N1*N-1)./2:(N1*N-1)./2;
k2=-(N2*N-1)./2:(N2*N-1)./2;

figure(2)
subplot(2,2,1);
stem(n1,w1);
title('x_3_1[n]')
subplot(2,2,2);
plot(k1,W1_db);
title('X_3_1(/Omega)')
subplot(2,2,3);
stem(n2,w2);
title('x_6_1[n]')
subplot(2,2,4);
plot(k2,W2_db);
title('X_6_1(/Omega)')

% APARTADO C)
clc

L=21;
n_i=0:L-1;

w_i=triang(L);
w_p=[ones(1,((L+1)/2)) zeros(1, (L+1)/2-1)];

W_i=abs(fftshift(fft(w_i,N*L)));
W_p=abs(fftshift(fft(w_p,N*L)));

Omega=linspace(-pi,pi,N*L+1);
Omega=Omega(1:N*L);

W_rec_dB = 10*log10(W_p./max(W_p));
W_tri_dB = 10*log10(W_i./max(W_i));

figure (3)
subplot(2,2,1)
stem(n_i,w_p)
ylabel('w_{rect}[n]')
subplot(2,2,2)
stem(n_i,w_i)
ylabel('w_{triang}[n]')
subplot(2,2,3)
plot(Omega,W_rec_dB)
ylabel('W_{rect}(\Omega)')
axis([-pi pi -30 1])
subplot(2,2,4)
plot(Omega,W_tri_dB)
ylabel('W_{triang}(\Omega)')
axis([-pi pi -30 1])

% APARTADO D)
L=41;
n_i=0:L-1;

w_i=triang(L);
w_p=[ones(1,((L+1)/2)) zeros(1, (L+1)/2-1)];

W_i=abs(fftshift(fft(w_i,N*L)));
W_p=abs(fftshift(fft(w_p,N*L)));

Omega=linspace(-pi,pi,N*L+1);
Omega=Omega(1:N*L);

W_rec_dB = 10*log10(W_p./max(W_p));
W_tri_dB = 10*log10(W_i./max(W_i));

figure(4)
subplot(2,2,1)
stem(n_i,w_p)
ylabel('w_{rect}[n]')
subplot(2,2,2)
stem(n_i,w_i)
ylabel('w_{triang}[n]')
subplot(2,2,3)
plot(Omega,W_rec_dB)
ylabel('W_{rect}(\Omega)')
axis([-pi pi -30 1])
subplot(2,2,4)
plot(Omega,W_tri_dB)
ylabel('W_{triang}(\Omega)')
axis([-pi pi -30 1])

L=61;
n_i=0:L-1;

w_i=triang(L);
w_p=[ones(1,((L+1)/2)) zeros(1, (L+1)/2-1)];

W_i=abs(fftshift(fft(w_i,N*L)));
W_p=abs(fftshift(fft(w_p,N*L)));

Omega=linspace(-pi,pi,N*L+1);
Omega=Omega(1:N*L);

W_rec_dB = 10*log10(W_p./max(W_p));
W_tri_dB = 10*log10(W_i./max(W_i));

figure(5)
subplot(2,2,1)
stem(n_i,w_p)
ylabel('w_{rect}[n]')
subplot(2,2,2)
stem(n_i,w_i)
ylabel('w_{triang}[n]')
subplot(2,2,3)
plot(Omega,W_rec_dB)
ylabel('W_{rect}(\Omega)')
axis([-pi pi -30 1])
subplot(2,2,4)
plot(Omega,W_tri_dB)
ylabel('W_{triang}(\Omega)')
axis([-pi pi -30 1])
%% EJERCICIO 1.3.9
clear all
close all
clc
%% APARTADO A) L = 64
clear all
close all
clc
%Generamos las se�ales
N = 1024; %Long de la se�al
L = 64; %Long. de la ventana
w1 = 2*pi/5;
%Vector que contiene los valores de DeltaOmega
vector_dw = 2*pi*[1/128 1/64 1/32 1/16];
for k = 1:length(vector_dw)
    n = 0:(N-1);
    dw = vector_dw(k);
    %phase1 = 2*pi*rand;
    %phase2 = 2*pi*rand;
    phase1 = 0; %A partir de la 1� simulacion, comentar
    phase2 = 0; %A partir de la 1� simulacion, comentar
    w2 = w1 + dw;
    y = exp(j*w1*n+phase1) + exp(j*w2*n+phase2);
    %Generamos la ventana y la se�al enventanada
    w_rec_64 = [ones(1,L), zeros(1,N-L)];
    y_rec_64 = y.*w_rec_64;
    %Calculamos las TF a partir de las DTFs.
    Omega = linspace(-pi,pi,L+1);
    Omega = Omega(1:L);
	 Y = fftshift(fft(y,L));
    W_rec_64 = fftshift(fft(w_rec_64,L));
    Y_W_rec_64 = fftshift(fft(y_rec_64,L)); 
    %Dibujamos los resultados
    figure
    subplot(2,3,1)
    stem(n(1:256),real(y(1:256)))
    ylabel('Re{y[n]}')
    subplot(2,3,2)
    stem(n(1:256),w_rec_64(1:256))
    ylabel(['w[n]: L=', num2str(L) ])
    subplot(2,3,3)
    stem(n(1:256),real(y_rec_64(1:256)))
    ylabel('Re{y_w[n]}')
    subplot(2,3,4)
    plot(Omega,abs(Y))
    ylabel('|Y(\Omega)|')
    axis([-pi pi min(Y) max(Y)])
    subplot(2,3,5)
    plot(Omega,abs(W_rec_64))
    ylabel('|W(\Omega)|')
    axis([-pi pi min(W_rec_64) max(W_rec_64)])
    subplot(2,3,6)
    plot(Omega,abs(Y_W_rec_64))
    ylabel('|Y_W(\Omega)|')
    axis([-pi pi min(Y_W_rec_64) max(Y_W_rec_64) ])
end
%% APARTADO B) L = 128
clear all
close all
clc
%Generamos las se�ales
N = 1024; %Long de la se�al
L = 128; %Long. de la ventana
w1 = 2*pi/5;
%Vector que contiene los valores de DeltaOmega
vector_dw = 2*pi*[1/128 1/64 1/32 1/16];
for k = 1:length(vector_dw)
    n = 0:(N-1);
    dw = vector_dw(k);
    phase1 = 2*pi*rand;
    phase2 = 2*pi*rand;
    %phase1 = 0; %A partir de la 1� simulacion, comentar
    %phase2 = 0; %A partir de la 1� simulacion, comentar
    w2 = w1 + dw;
    y = exp(j*w1*n+phase1) + exp(j*w2*n+phase2);
    %Generamos la ventana y la se�al enventanada
    w_rec_64 = [ones(1,L), zeros(1,N-L)];
    y_rec_64 = y.*w_rec_64;
    %Calculamos las TF a partir de las DTFs.
    Omega = linspace(-pi,pi,L+1);
    Omega = Omega(1:L);
	 Y = fftshift(fft(y,L));
    W_rec_64 = fftshift(fft(w_rec_64,L));
    Y_W_rec_64 = fftshift(fft(y_rec_64,L)); 
    %Dibujamos los resultados
    figure
    subplot(2,3,1)
    stem(n(1:256),real(y(1:256)))
    ylabel('Re{y[n]}')
    subplot(2,3,2)
    stem(n(1:256),w_rec_64(1:256))
    ylabel(['w[n]: L=', num2str(L) ])
    subplot(2,3,3)
    stem(n(1:256),real(y_rec_64(1:256)))
    ylabel('Re{y_w[n]}')
    subplot(2,3,4)
    plot(Omega,abs(Y))
    ylabel('|Y(\Omega)|')
    axis([-pi pi min(Y) max(Y)])
    subplot(2,3,5)
    plot(Omega,abs(W_rec_64))
    ylabel('|W(\Omega)|')
    axis([-pi pi min(W_rec_64) max(W_rec_64)])
    subplot(2,3,6)
    plot(Omega,abs(Y_W_rec_64))
    ylabel('|Y_W(\Omega)|')
    axis([-pi pi min(Y_W_rec_64) max(Y_W_rec_64) ])
end

%% APARTADO C) N = 256
clear all
close all
clc
%Generamos las se�ales
N = 1024; %Long de la se�al
L = 256; %Long. de la ventana
w1 = 2*pi/5;
%Vector que contiene los valores de DeltaOmega
vector_dw = 2*pi*[1/128 1/64 1/32 1/16];
for k = 1:length(vector_dw)
    n = 0:(N-1);
    dw = vector_dw(k);
    phase1 = 2*pi*rand;
    phase2 = 2*pi*rand;
    %phase1 = 0; %A partir de la 1� simulacion, comentar
    %phase2 = 0; %A partir de la 1� simulacion, comentar
    w2 = w1 + dw;
    y = exp(j*w1*n+phase1) + exp(j*w2*n+phase2);
    %Generamos la ventana y la se�al enventanada
    w_rec_64 = [ones(1,L), zeros(1,N-L)];
    y_rec_64 = y.*w_rec_64;
    %Calculamos las TF a partir de las DTFs.
    Omega = linspace(-pi,pi,L+1);
    Omega = Omega(1:L);
	 Y = fftshift(fft(y,L));
    W_rec_64 = fftshift(fft(w_rec_64,L));
    Y_W_rec_64 = fftshift(fft(y_rec_64,L)); 
    %Dibujamos los resultados
    figure
    subplot(2,3,1)
    stem(n(1:256),real(y(1:256)))
    ylabel('Re{y[n]}')
    subplot(2,3,2)
    stem(n(1:256),w_rec_64(1:256))
    ylabel(['w[n]: L=', num2str(L) ])
    subplot(2,3,3)
    stem(n(1:256),real(y_rec_64(1:256)))
    ylabel('Re{y_w[n]}')
    subplot(2,3,4)
    plot(Omega,abs(Y))
    ylabel('|Y(\Omega)|')
    axis([-pi pi min(Y) max(Y)])
    subplot(2,3,5)
    plot(Omega,abs(W_rec_64))
    ylabel('|W(\Omega)|')
    axis([-pi pi min(W_rec_64) max(W_rec_64)])
    subplot(2,3,6)
    plot(Omega,abs(Y_W_rec_64))
    ylabel('|Y_W(\Omega)|')
    axis([-pi pi min(Y_W_rec_64) max(Y_W_rec_64) ])
end