%% Practica 1-6
%% EJERCICIO 1.6.1
clear all
close all
clc
%Generamos un proceso gaussiano iid con randn. sqrt(var).*randn(1, l) +mean
var = 4;
k = 50;
x_iid = sqrt(var).*randn(k, 1);
n_xiid = 0:k-1;
%Estimamos la correlaci�n con xcorr. Utilizamos xcorr con las cuatro
%opciones posibles (ver doc xcorr). Las m�s comunes son unbiased y coeff 
[R1,k_i1] = xcorr(x_iid,'unbiased');
[R2,k_i2] = xcorr(x_iid,'biased');
[R3,k_i3] = xcorr(x_iid,'coeff');
[R4,k_i4] = xcorr(x_iid,'none');
%Dibujamos la se�al y las cuatro correlaciones
figure
subplot(2,3,1)
stem(n_xiid,x_iid)
ylabel('x[n] gauss iid')
subplot(2,3,2)
stem(k_i1,R1)
ylabel('R[k]: unbiased')
subplot(2,3,3)
stem(k_i2,R2)
ylabel('R[k]: biased')
subplot(2,3,5)
stem(k_i3,R3)
ylabel('R[k]: coeff')
subplot(2,3,6)
stem(k_i4,R4)
ylabel('R[k]: none')

%% EJERCICIO 1.6.2
clear all
close all
clc

%Generamos un proceso gaussiano iid con randn. 
var = 4;
k = 50;
x_iid = sqrt(var).*randn(k,1);
n_xiid = 0:k-1;
%Filtramos dicha secuencia (filtro MA s�lo coeff b y filtro AR un �nico
%coef a -adem�s del a(0), que siempre hay que ponerlo-) A = Y; B =X
a_MA = [1];
b_MA = [1 1 1 1 1 1 1]; %x[n] + x[n-1] + x[n-2] + x[n-3] + x[n-4] + x[n-5] + x[n-6];
x_MA = filter(b_MA,a_MA,x_iid);

a_AR = [1 0.8];
b_AR = [1];
x_AR = filter(b_AR,a_AR,x_iid);
%Respuestas al impulso correspondientes (obtenidas de forma te�rica)
h_MA = [ones(1,7), zeros(1,k-7)];
h_AR = (0.8).^n_xiid;
%Nota: las respuestas al impulso tambi�n podr�an haberse obtenido de forma 
%'pr�ctica' con las siguientes 3 l�neas
% impulso = [1, zeros(1,50)];
% h_MA = filter(ones(1,7),1,impulso);
% h_AR = filter(1,[1 -.8],impulso);
figure
subplot(2,3,1)
stem(n_xiid,x_iid)
ylabel('x[n] gauss iid')
subplot(2,3,2)
stem(n_xiid,x_MA)
ylabel('x[n] gauss MA')
subplot(2,3,3)
stem(n_xiid,x_AR)
ylabel('x[n] gauss AR')
subplot(2,3,5)
stem(n_xiid,h_MA)
ylabel('h[n] filtro MA')
subplot(2,3,6)
stem(n_xiid,h_AR)
ylabel('h[n] filtro AR')

%%Estimamos la AUTOcorrelaci�n con xcorr. 
[R_iid,k_i1] = xcorr(x_iid,'biased');
[R_MA,k_i2] = xcorr(x_MA,'biased');
[R_AR,k_i3] = xcorr(x_AR,'biased');
%%Para obtener la primera estimaci�n de la DEP, calculamos tambi�n su TF
Omega1 = linspace(-pi,pi, k);
S_iid1 = abs(fftshift(fft(R_iid,k))).^2;
S_MA1 = abs(fftshift(fft(R_MA,k))).^2;
S_AR1 = abs(fftshift(fft(R_AR,k))).^2;
%Para obtener la segunda estimaci�n de la DEP, utilizamos pwelch 
[S_iid2, Omega2] = pwelch(x_iid); 
[S_MA2, Omega2] = pwelch(x_MA);
[S_AR2,Omega2] = pwelch(x_AR);
%Podemos comparar estos resultado con un tercer m�todo para estimar
%la DEP que calcula el periodograma considerando la totalidad de la
%se�al como el �nico segmento (es decir, equivale a hacer 
%la DFT y elevarla alcuadrado)
[S_iid3, Omega3] = periodogram(x_iid); 
[S_MA3, Omega3] = periodogram(x_MA);
[S_AR3,Omega3] = periodogram(x_AR);
%Dibujamos la se�al y las cuatro correlaciones
figure
subplot(4,3,1)
stem(n_xiid,x_iid)
ylabel('x[n] gauss iid')
subplot(4,3,2)
stem(n_xiid,x_MA)
ylabel('x[n] gauss MA')
subplot(4,3,3)
stem(n_xiid,x_AR)
ylabel('x[n] gauss AR')
subplot(4,3,4)
plot(Omega1,S_iid1)
ylabel('S_x(\Omega): TF{R_i_i_d}')
subplot(4,3,5)
plot(Omega1,S_MA1)
ylabel('S_x(\Omega): TF{R_M_A}')
subplot(4,3,6)
plot(Omega1,S_AR1)
ylabel('S_x(\Omega): TF{R_A_R}')
subplot(4,3,7)
plot(Omega2,S_iid2)
ylabel('S_x(\Omega): Period')
subplot(4,3,8)
plot(Omega2,S_MA2)
ylabel('S_x(\Omega): Period')
subplot(4,3,9)
plot(Omega2,S_AR2)
ylabel('S_x(\Omega): Period')
subplot(4,3,10)
plot(Omega3,S_iid3)
ylabel('S_x(\Omega): Period 1 shot')
subplot(4,3,11)
plot(Omega3,S_MA3)
ylabel('S_x(\Omega): Period 1 shot')
subplot(4,3,12)
plot(Omega3,S_AR3)
ylabel('S_x(\Omega): Period 1 shot')

%%Calculamos te�ricamente la autocorrelaci�n
%Para el caso i.i.d. es trivial 
R_iid_analit = [zeros(1,k-1), 1, zeros(1,k-1)];
%Para los otros dos casos, convolucionamos h[n]con h[-n] (help fliplr)
R_MA_analit = conv(h_MA, fliplr(h_MA));
R_AR_analit = conv(h_AR, fliplr(h_AR));
%Generamos el �ndice temporal
k_ind_analit = -(k-1):(k-1);
%Estimamos las TF correspondientes usando fft
k2 = 256;
S_iid_analit  = abs(fftshift(fft(R_iid_analit,k2))).^2;
S_MA_analit  = abs(fftshift(fft(R_MA_analit,k2))).^2;
S_AR_analit  = abs(fftshift(fft(R_AR_analit,k2))).^2;
%Generamos el eje de frecuencias: 256 muestras en el intervalo [0,pi)
Omega_analit = linspace(0, pi, k2);
%Dibujamos las se�ales, las tres autocorrelaciones y las tres DEP
figure
subplot(3,3,1)
stem(n_xiid,x_iid)
ylabel('x[n] gauss iid')
subplot(3,3,2)
stem(n_xiid,x_MA)
ylabel('x[n] gauss MA')
subplot(3,3,3)
stem(n_xiid,x_AR)
ylabel('x[n] gauss AR')
subplot(3,3,4)
stem(k_ind_analit,R_iid_analit)
ylabel('R_i_i_d[k]')
subplot(3,3,5)
stem(k_ind_analit,R_MA_analit)
ylabel('R_M_A[k]')
subplot(3,3,6)
stem(k_ind_analit,R_AR_analit)
ylabel('R_A_R[k]')
subplot(3,3,7)
plot(Omega_analit,S_iid_analit)
ylabel('S_x(\Omega): Analit')
subplot(3,3,8)
plot(Omega_analit,S_MA_analit)
ylabel('S_x(\Omega): Analit')
subplot(3,3,9)
plot(Omega_analit,S_AR_analit)
ylabel('S_x(\Omega): Analit')

