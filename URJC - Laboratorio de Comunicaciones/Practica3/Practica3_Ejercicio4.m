%% Ejercicio 4
%% Enunciado
clc
clear all
close all
%OFDM con canal Multipath y Rayleigh


% TRANSMISOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Definiciones
N = 64;
L = 5;
N_bits = 2*N*100;
snr = 10;
 
%Transmisor habitual
Ax = modem.qammod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
b_v = round(rand(N_bits,1));
x_v = modulate(Ax, b_v);
 
%Procesamiento OFDM
x_m = reshape(x_v,N,ceil(length(x_v)/N));
X_m = sqrt(N)*ifft(x_m,N);
X_m_cp = [X_m((N-(L-1)):N,:); X_m];
X_v_cp = X_m_cp(:);

% CANAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perfil de taps
h = [1,1,1,1]; % En amplitud
h_pot = abs(h).^2; % En potencia;
h_pot_dB = 10*log10(h_pot); % En potencia y en dB
%Rayleigh (sin Multipath)
ch = rayleighchan(1e-5,10,(0:3)*1e-5,h_pot_dB);
ch.StorePathGains = 1;
HX_v_cp = filter(ch,X_v_cp);
H = ch.PathGains;  
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));

%RECEPTOR: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Procesamiento OFDM
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
H_m_cp = reshape(H,L+N,ceil(length(X_v_cp)/(L+N)),length(h_pot_dB));
Y_m = Y_m_cp(L+1:end,:);
H_m = H_m_cp(L+1:end,:,:); 
y_m = sqrt(1/N)*fft(Y_m,N);

%Eliminamos el efecto del canal
H_mean = mean(H_m,1);  
H_mean = squeeze(H_mean); 
h_m = fft(H_mean.',N);                      
y_m = y_m./h_m;
y_v = reshape(y_m,1,length(x_v));
 
%Receptor (detecci�n) habitual
brx_v = demodulate(modem.qamdemod(Ax),y_v);
brx_v = brx_v(:);

%Estimamos la BER
disp(['snr = ', num2str(snr),'(canal multitrayecto y Rayleigh)  BER = ', num2str(mean(abs(brx_v-b_v)))]);
%% Ejercicio 4a
%inserto pilotos cada 5 instantes y 4 portadoras
H_p = h_m(1:4:N,1:5:100);
[X,Y] = meshgrid((1:4:N),(1:5:100));
%interpolacion
X = X.';
Y = Y.';
Xq = 1:100;
Yq = 1:N;
V = interp2(Y,X,H_p,Xq,Yq.');   
for k =1:3;%portadoras
    V(61+k,:) = V(61,:);
end
for k = 1:4 %instantes
    V(:,96+k) = V(:,96);
end
% Representamos
figure
subplot(2,1,1)
imagesc(abs(h_m));
title('Canal original')
subplot(2,1,2)
imagesc(abs(V));
title('Canal estimado')
ECM_a = abs(mean(mean((h_m-V).^2)))
%% Ejercicio 4b 
H_p = h_m(1:16:64,1:5:100);%inserto pilotos cada 5 instantes y 16 portadoras
[X,Y] = meshgrid((1:16:64),(1:5:100));
X = X.';Y = Y.';
Xq = 1:100;
Yq = 1:64;
V = interp2(Y,X,H_p,Xq,Yq.');
% Representamos
for k =1:15;%portadoras
    V(49+k,:) = V(49,:);
end
for k = 1:4 %instantes
    V(:,96+k) = V(:,96);
end
figure
subplot(2,1,1)
imagesc(abs(h_m));
title('Canal original')
subplot(2,1,2)
imagesc(abs(V));
title('Canal estimado')
ECM_b = abs(mean(mean((h_m-V).^2)))
