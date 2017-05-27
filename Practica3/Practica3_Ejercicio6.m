%% Ejercicio 6
clc
clear all
close all
%% Ejercicio 6a
Mtx = 3;
Mrx = 5;
N = 1e5;
snr=10^(0.6);

H = sqrt(snr/2)*(randn(Mtx,Mrx,N) + 1i*randn(Mtx,Mrx,N));
%Buscamos el m�ximo canal de los 15 posibles en cada instante
h_max1 =max(H,[],1);
h_max2=max(h_max1);
h_max=squeeze(h_max2);

figure
% Histograma 3Tx y 5Rx
subplot(3,1,1)
hist(abs(h_max).^2,100)
title('Histograma SNR 3Tx y 5Rx')

%Histohrama 1Tx y 5Rx 
subplot(3,1,2)
h=H(1,:,:);
h = squeeze(h);
h_15 = max(h,[],1);
hist(abs(h_15).^2,100)
title('Histograma SNR 1Tx y 5Rx')


%Histograma 3Tx y 1Rx
subplot(3,1,3)
h=H(:,1,:);
h = squeeze(h);
h_31 = max(h,[],1);
hist(abs(h_31).^2,100)
title('Histograma SNR 3Tx y 1Rx')
xlabel('SNR')

%% Ejercicio 6b
Mtx = 3;
Mrx = 5;
N = 1e5;
snr=10^(0.6);

H = sqrt(snr/2)*(randn(Mrx,Mtx,N) + 1i*randn(Mrx,Mtx,N));
h_eq = zeros(1,N);

%Implementamos MIMO
for n=1:N
     %Hacemos descomposici�n SVD
     [U,S,V] = svd(H(:,:,n));
     
     v=V(:,1); %Vector por el que multiplicamos en tx (se lo tiene que enviar el Rx)
     u=U(:,1); %En el Rx multiplicamos por el conjugado de este
     %Canal equivalente
     h_eq(n) = (u')*H(:,:,n)*v; %En Matlab ' es herm�tico y .' transpuesto 
end

figure;
hist(abs(h_eq).^2,50);
xlabel('SNR MIMO');

%% Ejercicio 6c
Mtx = 3;
Mrx = 5;
N = 1e4;
snr=10^(0.6);

H = sqrt(snr/2)*(randn(Mtx,Mrx,N) + 1i*randn(Mtx,Mrx,N));
%Buscamos el m�ximo canal de los 15 posibles en cada instante
%caso a
    h_max1 = max(H,[],1);
    h_max2 = max(h_max1);
    h_max = squeeze(h_max2);
%caso b 
    h=H(1,:,:);
    h = squeeze(h);
    h_15 = max(h,[],1);
    h_15 = h_15.';
%caso c
    h=H(:,1,:);
    h = squeeze(h);
    h_31 = max(h,[],1);
    h_31 = h_31.';
    
    
    h_eq = zeros(1,N);
for n=1:N
     %Hacemos descomposici�n SVD
     [U,S,V] = svd(H(:,:,n));
     
     v=V(:,1); %Vector por el que multiplicamos en tx (se lo tiene que enviar el Rx)
     u=U(:,1); %En el Rx multiplicamos por el conjugado de este
     %Canal equivalente
     h_eq(n) = (u')*H(:,:,n)*v; %En Matlab ' es herm�tico y .' transpuesto
end
     h_eq = h_eq.';
         
     
N_bits = 3*N;
vector_SNR = [2 6 10 14 18];
ber = zeros(4,length(vector_SNR));
canales =[h_max h_15 h_31 h_eq];
Ax = modem.qammod('M',8,'SymbolOrder','gray', 'InputType', 'bit');

for m = 1:4
    h = canales(:,m);
    for k = 1:length(vector_SNR)
        SNR = vector_SNR(k);
        %Transmisor habitual
        b_v = round(rand(N_bits,1));
        x_v = modulate(Ax, b_v);
        E_Simb = mean(abs(x_v).^2);
        x_v = (1/sqrt(E_Simb)).*x_v;
        % Ruido
        snr = 10^(SNR/10);
        ruido = (randn(size(x_v)) + 1i*randn(size(x_v)))/sqrt(2*snr);
        % Canal
        y = x_v.*h + ruido;

        % Recepci�n
        y_r = y./h;

        % Demodulamos
        brx_v = demodulate(modem.qamdemod(Ax),y_r);
        brx_v = brx_v(:);
        ber(m,k) = mean(abs(brx_v-b_v));
    end
end
ber