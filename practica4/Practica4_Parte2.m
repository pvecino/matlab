%% HITO 2: CIS
clear all
close all
clc

Vector_SNR = 0:13;
m = 1;
for SNR = Vector_SNR
	snr = 10^(SNR/10);
	BER(m) = qfunc(sqrt(2*snr));
	m = m + 1;
end

semilogy(Vector_SNR, BER)
N_muestras = 100./min(BER);   % N = 100/p  ???? - coincide pe = 10^-v --> N~= 10^v+1

Nb = 1e6;
Vector_SNR =0:13;

bits = round(rand(Nb, 1));
Ax = modem.pskmod('M',2, 'SymbolOrder', 'binary', 'InputType', 'bit');
M = 2;   
N = Nb / log2(M); 

berBPSK = zeros(size(Vector_SNR));
berBPSK_CIS = zeros(size(Vector_SNR));

% MC
m = 1; 
for SNR = Vector_SNR
    %QPSK
    Simb_BPSK = modulate(Ax, bits);
    %RUIDO
    snr = 10^(SNR./10);
    ruido = sqrt(1/(2*snr)).*randn(1,Nb);
    y = Simb_BPSK + ruido';
    %----------------
    %   RECEPCION    
    %----------------
    bits_r = demodulate(modem.pskdemod(Ax),y);
    berBPSK(m) = sum((bits~=bits_r))./Nb;
    m=m+1; 
end
hold on
semilogy(Vector_SNR,berBPSK, 'r*')


%CIS
Vector_SNR = 0:13;
m = 1;
var = 10;                             %se introduce el error relativo
for SNR = Vector_SNR
	Simb_BPSK = modulate(Ax, bits);
	snr = 10^(SNR/10);
	ruido = sqrt(var/snr).*randn(1,Nb);
    W = peso(ruido, snr, var);

    y_cis = Simb_BPSK + ruido';
    %----------------
    %   RECEPCION    
    %----------------
	bits_r_cis = demodulate(modem.pskdemod(Ax),y_cis);
	berBPSK_CIS(m) = mean((bits ~= bits_r_cis).*W);
	m = m + 1;
end
hold on
semilogy(Vector_SNR,berBPSK_CIS, 'go')
legend('BER Teorica', 'BER MC', 'BER CIS')
title('BER vs SNR');
xlabel('SNR')
ylabel('BER')

