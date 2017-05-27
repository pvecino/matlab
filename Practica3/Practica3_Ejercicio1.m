%% Ejercicio 1.

%% Ejercicio 1a.
% Generamos canal exponencial.
n = 0:4;
pot_ultimo = 0.1;
taps = 4;
exponente_canal = log(sqrt(pot_ultimo))/taps;
h= exp(n*exponente_canal);
%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
 
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

%Multipath
h = h.';
h = sqrt(1/sum(abs(h).^2))*h; % OJO!!! Normalizamos para que la suma de las potencias de los multitrayectos valga 1.

HX_v_cp = filter(h,1,X_v_cp); % Pod�amos haberlo hecho con conv, pero as� sale directamente de las dimensiones correctas.
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
%Eliminamos el efecto del canal
h_m = fft(h,N);
h_m_rep = repmat(h_m,1,length(y_m(1,:))); 
y_m = y_m./h_m_rep;                       
%Bloques (conversor P/S)
y_v = reshape(y_m,1,length(x_v));
 
%Receptor (detecci�n) habitual
brx_v = demodulate(modem.pskdemod(Ax), y_v);
brx_v = brx_v(:);

%Estimamos la BER = 0.04996
disp(['snr = ', num2str(snr),'(canal Multipath)  BER = ', num2str(mean(abs(brx_v-b_v)))])

%% Ejercicio 1b.
% Generamos canal exponencial.
n = 0:4;
pot_ultimo = 0.1;
taps = 4;
exponente_canal = log(sqrt(pot_ultimo))/taps;
h= exp(n*exponente_canal);

%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
 
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

%Multipath
h = h.';
h = sqrt(1/sum(abs(h).^2))*h; 
HX_v_cp = filter(h,1,X_v_cp); 
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
%Eliminamos el efecto del canal
h = [h(1); zeros(4,1)];
h_m = fft(h,N);
h_m_rep = repmat(h_m,1,length(y_m(1,:))); % OJO!!! Es repmat, no reshape. Ver qu� estamos replicando h_m_rep(:,1:4)
y_m = y_m./h_m_rep;                       % Cada frecuencia (portadora) se ve afectada por un coeficiente; pero ese coeficiente es constante en el tiempo.
%Bloques (conversor P/S)
y_v = reshape(y_m,1,length(x_v));
 
%Receptor (detecci�n) habitual
brx_v = demodulate(modem.pskdemod(Ax), y_v);
brx_v = brx_v(:);

%Estimamos la BER = 0.17
disp(['snr = ', num2str(snr),'(considera solo rayo principal)  BER = ', num2str(mean(abs(brx_v-b_v)))])

%% Ejercicio 1c
% Generamos canal exponencial.
n = 0:4;
pot_ultimo = 0.1;
taps = 4;
exponente_canal = log(sqrt(pot_ultimo))/taps;
h= exp(n*exponente_canal);

%Definiciones
N = 64;
L = 8;
N_bits = 2*N*1000;
snr = 10;
 
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

%Sin multipath
HX_v_cp = X_v_cp;
%Ruido
Y_v_cp = HX_v_cp + sqrt(.5/snr) * (randn(size(HX_v_cp)) + 1i * randn(size(HX_v_cp)));

%Procesamiento OFDM
%Bloques (conversor S/P)
Y_m_cp = reshape(Y_v_cp,L+N,ceil(length(X_v_cp)/(L+N)));
%Prefijo c�clico (quitamos)
Y_m = Y_m_cp(L+1:end,:);
%DFT (directa)
y_m = sqrt(1/N)*fft(Y_m,N);
%Eliminamos el efecto del canal
h = [h(1); zeros(4,1)];
h_m = fft(h,N);
h_m_rep = repmat(h_m,1,length(y_m(1,:))); % OJO!!! Es repmat, no reshape. Ver qu� estamos replicando h_m_rep(:,1:4)
y_m = y_m./h_m_rep;                       % Cada frecuencia (portadora) se ve afectada por un coeficiente; pero ese coeficiente es constante en el tiempo.
%Bloques (conversor P/S)
y_v = reshape(y_m,1,length(x_v));
 
%Receptor (detecci�n) habitual
brx_v = demodulate(modem.pskdemod(Ax), y_v);
brx_v = brx_v(:);

%Estimamos la BER = 0.000867
disp(['snr = ', num2str(snr),'(AWGN)  BER = ', num2str(mean(abs(brx_v-b_v)))])