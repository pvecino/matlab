clear all
close all
clc

% Par�metros simulaci�n
Nb = 52000; % Numero de bits a transmitir.
vector_SNR = (0:2:10); % Relaciones se�al a ruido

%Transmisor habitual
M = 2; % Orden de la constelaci�n


% Generaci�n de los bits
b = round(rand(1,Nb));

%CODIFICACION DE CANAL
% señal sin codificar
Ax = modem.pskmod('M',2,'SymbolOrder','gray', 'InputType', 'bit');

%señal codificada Hamming r=3
r_r3 = 3;
[H_r3, G_r3, n_ham_r3, k_r3] = hammgen(r_r3);
Tabla_Decod_r3 = syndtable(H_r3);
Hamm_r3_Tx = encode(b, n_ham_r3, k_r3, 'linear', G_r3)


%señal codificada Hamming r=5
r_r5 = 5;
[H_r5, G_r5, n_ham_r5, k_r5] = hammgen(r_r5);
Tabla_Decod_r5 = syndtable(H_r5);
Hamm_r5_Tx = encode(b, n_ham_r5, k_r5, 'linear', G_r5)

%señal codificada bloque repeticion 3 bits
n_rep_3 = 3;
k_rep_3 = 1;
G_rep_3 = [1 1 1];
Rep3_Tx = encode(b ,n_rep_3 , k_rep_3, 'linear', G_rep_3)

%señal codificada bloque repeticion 5 bits
n_rep_5 = 5;
k_rep_5 = 1;
G_rep_5 = [1 1 1 1 1];
Rep5_Tx = encode(b ,n_rep_5 , k_rep_5, 'linear', G_rep_5)

N = Nb / log2(M); % N�mero de s�mbolos a transmitir para transmitir Nb bits.

berBPSK = zeros(size(vector_SNR));
berBPSK_hamm_r3 = zeros(size(vector_SNR));
berBPSK_hamm_r5 = zeros(size(vector_SNR));
berBPSK_rep_3 = zeros(size(vector_SNR));
berBPSK_rep_5 = zeros(size(vector_SNR));
m = 1;


for SNR = vector_SNR,
m
    % Generaci�n de los bits
    %b = round(rand(1,Nb));

    % Modulaci�n
    x = modulate(Ax, b);
    x_r3 = modulate(Ax, Hamm_r3_Tx)
    x_r5 = modulate(Ax, Hamm_r5_Tx)
    x_rep_3 = modulate(Ax, Rep3_Tx)
    x_rep_5 = modulate(Ax, Rep5_Tx)

    % Generaci�n del ruido
    snr = 10^(SNR/10);
    n = (randn(1,N) + j.*randn(1,N)) * (1/sqrt(2)) * sqrt(1/snr);
    n_r3 = (randn(1,length(x_r3)) + j.*randn(1,length(x_r3))) * (1/sqrt(2)) * sqrt(1/snr);
    n_r5 = (randn(1,length(x_r5)) + j.*randn(1,length(x_r5))) * (1/sqrt(2)) * sqrt(1/snr);
    n_rep3 = (randn(1,length(x_rep_3)) + j.*randn(1,length(x_rep_3))) * (1/sqrt(2)) * sqrt(1/snr);
    n_rep5 = (randn(1,length(x_rep_5)) + j.*randn(1,length(x_rep_5))) * (1/sqrt(2)) * sqrt(1/snr);

    % Introducimos la se�al en el canal.
    y = x + n;
    y_r3 = x_r3 + n_r3;
    y_r5 = x_r5 + n_r5;
    y_rep3 = x_rep_3 + n_rep3;
    y_rep5 = x_rep_5 + n_rep5;
    
    % Demodulador
    % sin codificacion de canal
    br = demodulate(modem.pskdemod(Ax),y);
    
    %con condificacion de canal Hamming r = 3
    br_r3 = demodulate(modem.pskdemod(Ax),y_r3);
    
    %con condificacion de canal Hamming r = 5
    br_r5 = demodulate(modem.pskdemod(Ax),y_r5);
    
    %con condificacion de canal Repeticion 3 bits
    br_rep3 = demodulate(modem.pskdemod(Ax),y_rep3);
    
    %con condificacion de canal Repeticion 5 bits
    br_rep5 = demodulate(modem.pskdemod(Ax),y_rep5);
    
    % DECODIFICACION DE CANAL
    Hamm_r3_RX = decode(br_r3, n_ham_r3, k_r3, 'linear', G_r3, Tabla_Decod_r3)
    Hamm_r5_RX = decode(br_r5, n_ham_r5, k_r5, 'linear', G_r5, Tabla_Decod_r5)
    Rep_3_bits_RX = decode(br_rep3, n_rep_3, k_rep_3, 'linear', G_rep_3)
    Rep_5_bits_RX = decode(br_rep5, n_rep_5, k_rep_5, 'linear', G_rep_5)
    

    % C�lculo de la BER
    berBPSK(m) = sum((b~=br))./Nb;
    berBPSK_hamm_r3(m) = sum((b~=Hamm_r3_RX))./Nb;
    berBPSK_hamm_r5(m) = sum((b~=Hamm_r5_RX))./Nb;
    berBPSK_rep_3(m) = sum((b~=Rep_3_bits_RX))./Nb;
    berBPSK_rep_5(m) = sum((b~=Rep_5_bits_RX))./Nb;
   
    m=m+1;    
end

semilogy(berBPSK)
hold on
semilogy(berBPSK_hamm_r3, 'r')
hold on
semilogy(berBPSK_hamm_r5, 'g')
hold on
semilogy(berBPSK_rep_3, 'k')
hold on
semilogy(berBPSK_rep_5, 'c')
legend('Sin codificacion', 'Hamming r = 3', 'Hamming r = 5', 'Cod. Repeticion 3 bits', 'Cod. Repeticion 5 bits')
grid on
