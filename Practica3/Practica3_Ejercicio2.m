%% Ejercicio 2.
clear all
close all
clc

%% Ejercicio 2a.
%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
Fd = 250;
 
%Transmisor habitual
Ax = modem.pskmod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
b_v = round(rand(N_bits,1));
x_v = modulate(Ax, b_v);
 
%Procesamiento OFDM
%Bloques (conversor S/P)
x_m = reshape(x_v,N,ceil(length(x_v)/N));
%DFT (inversa)
X_m = sqrt(N)*ifft(x_m,N);
%Prefijo ciclico (a�adimos)
X_m_cp = [X_m((N-(L-1)):N,:); X_m];
%Bloques (conversor P/S)
X_v_cp = X_m_cp(:);


% CANAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rayleigh (sin Multipath)
ch = rayleighchan(1e-5,Fd);
ch.StorePathGains = 1;
HX_v_cp = filter(ch,X_v_cp);    % Equivale a channel.filter
H = ch.PathGains;   % Guardamos los path gains.
                    % En teor�a, me deber�an salir iguales para cada portadora (por estar dentro del BW de coherencia),
                    % pero la realidad es que var�an (aunque poco, porque lo hacen lentamente).
                    % Ya veremos c�mo lo arreglamos (media, escoger un representante...).
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));


%RECEPTOR: hay que eliminar el efecto del canal %%%%%%%%%%%%%%%%%%%%%%%%%%

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
H_m_cp = reshape(H,L+N,ceil(length(X_v_cp)/(L+N))); % Voy haciendo lo mismo con los coeficientes del canal
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
H_m = H_m_cp(L+1:end,:); % Ver H_m(:,1:4). No deber�an variar en cada columna, pero en realidad lo hacen (aunque despacio)
                         % Variar el desplazamiento Doppler del rayleighchan para comprobar como var�an m�s o menos.
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
%Eliminamos el efecto del canal
H_mean = mean(H_m); % Tenemos que quedarnos con "un representante" por columna. En este caso, optamos por la media. Ver H_mean(1:4)
h_m = fft([H_mean;zeros(size(H_mean))],N); % Ponemos fila de 0s para que haga fft por columnas. Ver diferencia con fft(H_mean,N)
                                           % Ver como queda replicada h_m(:,1:4)
                                           % Todas las frecuencias (portadoras) se ven afectadas por un coeficiente; pero ese coeficiente var�a en el tiempo.
y_m = y_m./h_m;
%Bloques (conversor P/S)
y_v = reshape(y_m,1,length(x_v));
 
%Receptor (detecci�n) habitual
brx_v = demodulate(modem.pskdemod(Ax), y_v);
brx_v = brx_v(:);

%Estimamos la BER
disp(['snr = ', num2str(snr),'(canal Rayleigh)  BER = ', num2str(mean(abs(brx_v-b_v)))])

%% Ejercicio 2b.
%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
pos_port = randi(64, 1, 4);
Fd = 250;
 
%Transmisor habitual
Ax = modem.pskmod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
b_v = round(rand(N_bits,1));
x_v = modulate(Ax, b_v);
 
%Procesamiento OFDM
%Bloques (conversor S/P)
x_m = reshape(x_v,N,ceil(length(x_v)/N));
%DFT (inversa)
X_m = sqrt(N)*ifft(x_m,N);
%Prefijo ciclico (a�adimos)
X_m_cp = [X_m((N-(L-1)):N,:); X_m];
%Bloques (conversor P/S)
X_v_cp = X_m_cp(:);


% CANAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rayleigh (sin Multipath)
ch = rayleighchan(1e-5,Fd);
ch.StorePathGains = 1;
HX_v_cp = filter(ch,X_v_cp);    % Equivale a channel.filter
H = ch.PathGains;   % Guardamos los path gains.
                    % En teor�a, me deber�an salir iguales para cada portadora (por estar dentro del BW de coherencia),
                    % pero la realidad es que var�an (aunque poco, porque lo hacen lentamente).
                    % Ya veremos c�mo lo arreglamos (media, escoger un representante...).
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));


%RECEPTOR: hay que eliminar el efecto del canal %%%%%%%%%%%%%%%%%%%%%%%%%%

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
H_m_cp = reshape(H,L+N,ceil(length(X_v_cp)/(L+N))); % Voy haciendo lo mismo con los coeficientes del canal
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
H_m = H_m_cp(L+1:end,:); % Ver H_m(:,1:4). No deber�an variar en cada columna, pero en realidad lo hacen (aunque despacio)
                         % Variar el desplazamiento Doppler del rayleighchan para comprobar como var�an m�s o menos.
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
y_m_2 = y_m;

k=0;
figure
for i = pos_port
	k=k+1;
	snr_p = abs(y_m_2(k,:)).*10;
	subplot(2,2,k)
	hist(snr_p, 40)
	title(['Portadora numero ' num2str(i)])
end

%% Ejercicio 2c.
%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
pos_port = randi(64, 1, 4);
Fd = 250;
 
%Transmisor habitual
Ax = modem.pskmod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
b_v = round(rand(N_bits,1));
x_v = modulate(Ax, b_v);
 
%Procesamiento OFDM
%Bloques (conversor S/P)
x_m = reshape(x_v,N,ceil(length(x_v)/N));
%DFT (inversa)
X_m = sqrt(N)*ifft(x_m,N);
%Prefijo ciclico (añadimos)
X_m_cp = [X_m((N-(L-1)):N,:); X_m];
%Bloques (conversor P/S)
X_v_cp = X_m_cp(:);


% CANAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rayleigh (sin Multipath)
ch = rayleighchan(1e-5,Fd);
ch.StorePathGains = 1;
HX_v_cp = filter(ch,X_v_cp);    % Equivale a channel.filter
H = ch.PathGains;   % Guardamos los path gains.
                    % En teor�a, me deber�an salir iguales para cada portadora (por estar dentro del BW de coherencia),
                    % pero la realidad es que var�an (aunque poco, porque lo hacen lentamente).
                    % Ya veremos c�mo lo arreglamos (media, escoger un representante...).
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));


%RECEPTOR: hay que eliminar el efecto del canal %%%%%%%%%%%%%%%%%%%%%%%%%%

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
H_m_cp = reshape(H,L+N,ceil(length(X_v_cp)/(L+N))); % Voy haciendo lo mismo con los coeficientes del canal
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
H_m = H_m_cp(L+1:end,:); % Ver H_m(:,1:4). No deber�an variar en cada columna, pero en realidad lo hacen (aunque despacio)
                         % Variar el desplazamiento Doppler del rayleighchan para comprobar como var�an m�s o menos.
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
%Eliminamos el efecto del canal
H_mean = mean(H_m); % Tenemos que quedarnos con "un representante" por columna. En este caso, optamos por la media. Ver H_mean(1:4)
h_m = fft([H_mean;zeros(size(H_mean))],N); % Ponemos fila de 0s para que haga fft por columnas. Ver diferencia con fft(H_mean,N)
                                           % Ver como queda replicada h_m(:,1:4)
                                           % Todas las frecuencias (portadoras) se ven afectadas por un coeficiente; pero ese coeficiente var�a en el tiempo.
y_m = y_m./h_m;
 
%Receptor (detecci�n) habitual


%Estimamos la BER
num = 1;
for pos = pos_port
	%Bits originales.
	x_mp = x_m(pos,:);
	b = demodulate(modem.pskdemod(Ax), x_mp);
	b = b(:);
	%Bits recibidos.
	y_mp = y_m(pos,:);
	br = demodulate(modem.pskdemod(Ax), y_mp);
	br = br(:);
	
	disp(['snr = ', num2str(snr),'(portadora numero ', num2str(pos),  ') BER = ', num2str(mean(abs(br-b)))])
 	num=num+1;
end