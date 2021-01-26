%% Ejercicio 5
clc
clear all
close all
Mtx = 4;
N = 1e5;
snr=10^(0.6);

H = sqrt(snr/2)*(randn(Mtx,N) + 1i*randn(Mtx,N));
h_max=max(H,[],1);% te devuelve un vector fila con los valores maximos
                  % de los 4 posibles, en cada instante. 

%Representamos el histograma 4tx--1rx
figure
subplot(2,1,1)
hist(abs(h_max).^2,50)
title('Histograma 4Tx y 1Rx')
xlabel('SNR')

%Representamos histograma 1tx--1rx
h=H(1,:);
subplot(2,1,2)
hist(abs(h).^2,50)
title('Histograma 1Tx y 1Rx')
xlabel('SNR')

