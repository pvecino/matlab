%% HITO 3.1
close all
clear all
clc
N = 1e5;
Px = 1;
vector_SNR = (0:2:20); % Relaciones se�al a ruido
%Generamos los vectores de canal, ruido, y se�al despu�s del canal
    canal = [1, -0.3, 0.1];
    canal_iid = sqrt(1/2)*(randn(1,N/2) + 1i*randn(1,N/2));
    n = zeros(1,N);
    y = zeros(1,N);
    ber = zeros(size(vector_SNR));
    bits = round(rand(N,1));
    % Bucle sobre SNR.
    pruebas = [1 2 3];
for m = pruebas,
    if m == 1 
      k = 1;
        for SNR = vector_SNR,
            
            % Generaci�n de se�al (QPSK, Potencia Unidad)
            h = modem.pskmod('M', 4, 'SymbolOrder' ,'gray', 'InputType' ,'bit');
            x_mod = modulate(h,bits);
            %Generamos el canal 
            % Generamos el ruido
            snr = 10^(SNR/10);
            n = (sqrt(1/(2*snr)).*((randn(1,N/2) + 1i*randn(1,N/2))));
            pn = mean(abs(n).^2);
            %Pasamos la se�al por el canal
            y = x_mod + n';
            % C�lculo de la BER
            h_demod = modem.pskdemod('M',4,'SymbolOrder', 'gray', 'OutputType', 'Bit');
            y_demod = demodulate(h_demod, y);
            ber(k) = mean(y_demod~=bits);
        k=k+1;    
        end
    
        semilogy(vector_SNR,ber,'--mo');
        hold on;
        grid on;
        xlabel('SNR(dB)');
        ylabel('BER');
    
    elseif m == 2
        k = 1;
        for SNR = vector_SNR,
            
            % Generaci�n de se�al (QPSK, Potencia Unidad)
            h = modem.pskmod('M', 4, 'SymbolOrder' ,'gray', 'InputType' ,'bit');
            x_mod = modulate(h,bits);
            %Generamos el canal 
            % Generamos el ruido
            snr = 10^(SNR/10);
            n = (sqrt(1/(2*snr)).*((randn(1,N/2) + 1i*randn(1,N/2))));
            %Pasamos la se�al por el canal
            y = filter(canal,1,x_mod);
            y_n = y + n';

            % Demodulo la se�al para volver a obtener bits
            h_demod = modem.pskdemod('M',4,'SymbolOrder', 'gray', 'OutputType', 'Bit');
            y_demod = demodulate(h_demod, y_n);
            %Miro los errores que he cometido
            ber(k) = mean(y_demod~=bits);
            k=k+1;    
        end
        
        semilogy(vector_SNR,ber,'--ro');
        grid on;
        xlabel('SNR(dB)');
        ylabel('BER')
        hold on;
        title('QPSK en diferentes escenarios');
        
    elseif m == 3
        k = 1;
        for SNR = vector_SNR,
            
            % Generaci�n de se�al (QPSK, Potencia Unidad)
            h = modem.pskmod('M', 4, 'SymbolOrder' ,'gray', 'InputType' ,'bit');
            x_mod2 = modulate(h,bits);
            %Generamos el canal 
            % Generamos el ruido
            snr = 10^(SNR/10);
            n = (sqrt(1/(2*snr)).*((randn(1,N/2) + 1i*randn(1,N/2))));
            %Pasamos la se�al por el canal
            y = filter(canal,1,x_mod2);
            y_n = y + n';
            %Creo un ecualizador ZF de 3 etapas
            h1 = [-0.3 1 0;
                 0.1 -0.3 1;
                  0 0.1 -0.3];
            y_deseada = [0;1;0];
            coeff_zf = h1\y_deseada;
            %Uso el ecualizador
            y_ecualizada = filter(coeff_zf,1,y_n);
            % Demodulo la se�al para volver a obtener bits
            h_demod = modem.pskdemod('M',4,'SymbolOrder', 'gray', 'OutputType', 'Bit');
            y_demod = demodulate(h_demod, y_ecualizada);
            %Miro los errores que he cometido
            ber(k) = mean(y_demod~=bits);
            k=k+1;    
        end
    
        semilogy(vector_SNR,ber,'--bo');
        legend('Canal AWGN','Canal Multipath + AWGN sin ecualizar','Canal Multipath + AWGN con ZF 3 etapas');
        grid on;
        hold on;
    end
   
    
end
clear all
Nb = 5200;
Vector_SNR =0:2:20;
bits = round(rand(Nb, 1));
%GENERO LOS SIMBOLOS DE LA QPSK
Ax = modem.pskmod('M',4, 'SymbolOrder', 'binary', 'InputType', 'bit');
M = 4;    %M=4
N = Nb / log2(M); 
berQPSK = zeros(size(Vector_SNR));

%Canal Doopler con desvanecimiento Rayleigh sin ecualizador
m = 1; 
for SNR = Vector_SNR
    
    %MODULAMOS LA QPSK
    Simb_QPSK = modulate(Ax, bits);
    %GENERO EL RUIDO
    snr = 10^(SNR/10);
    ruido =(randn(1,N) + j*randn(1,N))/sqrt(2*snr);
    %CANAL RAYLEIGH
    h = (randn(1,N)+ j.*randn(1,N)) * sqrt(1/2);
    %CREO LA SEÑAL REAL
    y = Simb_QPSK.*h' + ruido';
    %----------------
    %   RECEPCION    
    %----------------
    % DEMODULO LA QPSK
   
    bits_r = demodulate(modem.pskdemod(Ax),y);
    % CALCULO LA BER
    berQPSK(m) = sum((bits~=bits_r))./Nb;   
    m=m+1; 
end
semilogy(Vector_SNR,berQPSK,'--go')



% Canal Doopler con desvanecimiento Rayleigh con ecualizador
m = 1;
for SNR = Vector_SNR
    
    %MODULAMOS LA QPSK
    Simb_QPSK = modulate(Ax, bits);
    %Simb_QPSK = pskmod(bits,M);
    %GENERO EL RUIDO
    snr = 10^(SNR/10);
    ruido =(randn(1,N) + j*randn(1,N))/sqrt(2*snr);
    %CANAL RAYLEIGH
    h =  (randn(1,N)+ j.*randn(1,N)) * sqrt(1/2); 
    %CREO LA SEÑAL REAL
    y = Simb_QPSK.*h' + ruido';
   
    %----------------
    %   RECEPCION    
    %----------------
    %ECUALIZADOR
    y = y./h';
    %DEMODULO LA QPSK
    bits_r = demodulate(modem.pskdemod(Ax),y);
    % CALCULO LA BER
    berQPSK(m) = sum((bits~=bits_r))./Nb;
    m=m+1; 
end

hold on
semilogy(Vector_SNR,berQPSK,'--ko')
legend('Canal AWGN','Canal Multipath + AWGN sin ecualizar','Canal Multipath + AWGN con ZF 3 etapas','Rayleigh sin Eq','Rayleigh con Eq')

xlabel('SNR(dB)');
ylabel('BER');



%% HITO 3.2
clear all
close all
clc

% Par�metros simulaci�n
Nb = 52000; % Numero de bits a transmitir.
vector_SNR = (0:2:20); % Relaciones se�al a ruido

%Transmisor habitual
M = 2; % Orden de la constelaci�n

% Generaci�n de los bits
b = round(rand(1,Nb));

%CODIFICACION DE CANAL
% señal sin codificar
Ax = modem.pskmod('M',2,'SymbolOrder','gray', 'InputType', 'bit');

Bx = modem.dpskmod('M',2,'SymbolOrder','gray', 'InputType', 'bit');

N = Nb / log2(M); % N�mero de s�mbolos a transmitir para transmitir Nb bits.

berBPSK = zeros(size(vector_SNR));
berBPSK_dbpsk = zeros(size(vector_SNR));

m = 1;


for SNR = vector_SNR,

    
    % Modulaci�n
    x_bpsk = modulate(Ax, b);
    x_dbpsk = modulate(Ax, b);
    

    % Generaci�n del ruido
    snr = 10^(SNR/10);
    n = (randn(1,N) + j.*randn(1,N)) * (1/sqrt(2)) * sqrt(1/snr);
    n_dbpsk = (randn(1,length(x_dbpsk)) + j.*randn(1,length(x_dbpsk))) * (1/sqrt(2)) * sqrt(1/snr);
    

    % Introducimos la se�al en el canal.
    y = x_bpsk + n;
    y_dbpsk = x_dbpsk + n_dbpsk;
    
    
    % Demodulador
    % sin codificacion de canal
    br = demodulate(modem.pskdemod(Ax),y);
    
    br_dbpsk = demodulate(modem.dpskdemod(Bx),y_dbpsk); 
    
    % C�lculo de la BER
    berBPSK(m) = sum((b~=br))./Nb;
    berBPSK_dbpsk(m) = sum((b~=br_dbpsk))./Nb;
   
   
    m=m+1;    
end

semilogy(vector_SNR, berBPSK, '--o')
hold on
semilogy(vector_SNR, berBPSK_dbpsk, '--ro')
legend('BPSK', 'DBPSK')
xlabel('SNR(dB)');
ylabel('BER');
grid on
