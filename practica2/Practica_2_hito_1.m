%% Hito 1: Modulaciones digitales: Diagramas de Ojo y COnstelaciones
%% HITO 1.1
%-------------------------
%       B/Q/8 PSK
%-------------------------
clc 
clear all
close all

N = 1000;
M = [2 4 8];
roff=0.5;
Fd=1;
Fs=10;
delay = 3;
propd = delay.*(Fs/Fd)+1;

for k = 1:length(M)
    %GENERO LA SEÑAL B/Q/8 PSK
    Simb_PSK = randi(M(k),1,N)-1;
    Mod_PSK = pskmod(Simb_PSK,M(k));
    %NORMALIZO LA SEÑAL
    E_PSK = mean(abs(Mod_PSK).^2);
    Mod_PSK = Mod_PSK./(sqrt(E_PSK));
    %FILTRO LA SEÑAL CON COSENO ALZADO
    Filt_PSK = rcosflt(Mod_PSK,Fd,Fs,'fir/normal',roff,delay);
    %CON ESTO SE QUITAN LAS COLAS DEL FILTRO DEL COSENO ALZADO
    Filt_PSK = Filt_PSK(propd:end-2*(propd-1),:); 
    eyediagram(Filt_PSK,Fs)
end

%-------------------------
%         64QAM
%-------------------------
K = 64;
%GENERO LA SEÑAL 64QAM
Simb_64QAM = randi(K,1,N)-1; 
Mod_64QAM = qammod(Simb_64QAM,K);
%NORMALIZAMOS LA SEÑAL
E_64QAM = mean(abs(Mod_64QAM).^2);
Mod_64QAM = Mod_64QAM./(sqrt(E_64QAM));
%FILTRO CON COSENO ALZADO
Filt_64QAM = rcosflt(Mod_64QAM,Fd,Fs,'fir/normal',roff,delay);
%CON ESTO SE QUITAN LAS COLAS DEL FILTRO DEL COSENO ALZADO
Filt_64QAM =Filt_64QAM(propd:end-2*(propd-1),:);
%REPRESENTO EL DIAGRAMA DEL OJO
eyediagram(Filt_64QAM,Fs)

%% HITO 1.2
clc
close all
clear all
%-------------------------
%         QPSK
%-------------------------
N = 1000;
M = 4; %QPSK
roff=0.5;
Fd=1;
Fs=10;
delay = 3;
propd = delay.*(Fs/Fd)+1;
%GENERO DE LA SEÑAL QPSK
Simb_QPSK = randi(M,1,N)-1;
Mod_QPSK = pskmod(Simb_QPSK,M);
%normalizo la señal
E_QPSK = mean(abs(Mod_QPSK).^2);
Mod_QPSK = Mod_QPSK./(sqrt(E_QPSK));
%FILTRO LA SEÑAL CON COSENO ALZADO
Filt_QPSK = rcosflt(Mod_QPSK,Fd,Fs,'fir/normal',roff,delay);
%CON ESTO SE QUITAN LAS COLAS DEL FILTRO DEL COSENO ALZADO
Filt_QPSK = Filt_QPSK(propd:end-2*(propd-1),:); 
%UNIDADES NATURALES 20, 10 Y 3
SNR =[100,10,2]; 
for n=1:length(SNR)
    %GENERO EL RUIDO Y LE SUMO A LA SEÑAL FILTRADA
    ruido = (randn(size(Filt_QPSK)) + j*randn(size(Filt_QPSK)))/sqrt(2*SNR(n));
    S_QPSK = Filt_QPSK + ruido;
    eyediagram(S_QPSK,Fs)
end
%-------------------------
%         64QAM
%-------------------------
M = 64;
%GENERO LA SEÑAL 64QAM
Simb_64QAM = randi(M,1,N)-1;
Mod_64QAM = qammod(Simb_64QAM,M);
%NORMALIZO LA SEÑAL
E_64QAM = mean(abs(Mod_64QAM).^2);
Mod_64QAM = Mod_64QAM./(sqrt(E_64QAM));
%FILTRO CON COSENO ALZADO
Filt_64QAM = rcosflt(Mod_64QAM,Fd,Fs,'fir/normal',roff,delay);
%CON ESTO SE QUITAN LAS COLAS DEL FILTRO DEL COSENO ALZADO
Filt_64QAM =Filt_64QAM(propd:end-2*(propd-1),:);
% UNIDADES NATURALES!!!
SNR =[100,10,2]; 
for n=1:length(SNR) 
    %GENERO EL RUIDO Y LE SUMO A LA SEÑAL FILTRADA
    ruido = (randn(size(Filt_64QAM)) + j*randn(size(Filt_64QAM)))/sqrt(2*SNR(n));
    S_64QAM = Filt_64QAM + ruido;
    eyediagram(S_64QAM,Fs)
end

%% HITO 1.3
close all
clear all
clc 
%-------------------------
%       B/Q/8 PSK
%-------------------------
N = 1000;
roff=0.5;
Fd=1;
Fs=10;
delay = 3;
propd = delay.*(Fs/Fd)+1;
M = [2,4,8];
%UNIDADES NATURALES DE LA SNR DE 20, 10 Y 3 db
SNR = [100,10];

fig=1;%para poder hacer hold on en scaterplott
for m = 1:length(M)
    for n = 1:length(SNR)
        %GENERO LA SEÑAL B/Q/8 PSK
        Simb_PSK = randi(M(m),1,N)-1;
        Mod_PSK = pskmod(Simb_PSK,M(m));
         %NORMALIZO LA SEÑAL
        E_PSK = mean(abs(Mod_PSK).^2);
        Mod_PSK = Mod_PSK./(sqrt(E_PSK));

        %GENERO Y AÑADO RUIDO
        ruido = (randn(size(Mod_PSK)) + j*randn(size(Mod_PSK)))/sqrt(2*SNR(n));
        S_PKS = Mod_PSK + ruido;
        
        scatterplot(S_PKS,1,0,'*b');
        hold on; 
        scatterplot(Mod_PSK,1,0,'*r',fig)
        fig=fig+1;
    end   
end
%-------------------------
%         64QAM
%-------------------------
M = 64;
%GENERO LA SEÑAL 64QAM
Simb_64QAM = randi(M,1,N)-1;
Mod_64QAM = qammod(Simb_64QAM,M);
%NORMALIZO LA SEÑAL
E_64QAM = mean(abs(Mod_64QAM).^2);
Mod_64QAM = Mod_64QAM./(sqrt(E_64QAM));
% UNIDADES NATURALES!!!
SNR =[1000,10]; 
for n=1:length(SNR) 
    %GENERO EL RUIDO Y LE SUMO A LA SEÑAL FILTRADA
    ruido = (randn(size(Mod_64QAM)) + j*randn(size(Mod_64QAM)))/sqrt(2*SNR(n));
    S_64QAM = Mod_64QAM + ruido;
    %rEPRESENTAR LA CONSTELACION
    scatterplot(S_64QAM,1,0,'*b');
    hold on; 
    scatterplot(Mod_64QAM,1,0,'*r',fig)
    fig=fig+1;
end
%% HITO 1.4
close all
clear all
clc
%-------------------------
%         QPSK
%-------------------------
N = 1000;
M = 4; %QPSK
roff=0.5;
Fd=1;
Fs=10;
delay = 3;
propd = delay.*(Fs/Fd)+1;
SNR = 100;
%GENERO DE LA SEÑAL QPSK
Simb_QPSK = randi(M,1,N)-1;
Mod_QPSK = pskmod(Simb_QPSK,M);
%VECTOR DE MULTITRAYECTO
h2 =[1 -0.4 0.25];
h1 = [1 0.3];
%CONVOLUCIONO EL MULTITRAYECTO Y LA QPSK
Filt_QPSK_1 = filter(h1,1,Mod_QPSK);
scatterplot(Filt_QPSK_1,1,0,'*b');
Filt_QPSK_2 = filter(h2,1,Mod_QPSK);
scatterplot(Filt_QPSK_2,1,0,'*b');
%-------------------------
%         64QAM
%-------------------------
M = 64;
%GENERO LA SEÑAL 64QAM
Simb_64QAM = randi(M,1,N)-1;
Mod_64QAM = qammod(Simb_64QAM,M);
%NORMALIZO LA SEÑAL
E_64QAM = mean(abs(Mod_64QAM).^2);
Mod_64QAM = Mod_64QAM./(sqrt(E_64QAM));
%CONVOLUCIONO EL MULTITRAYECTO Y LA QPSK
Filt_64QAM_1 = filter(h1,1,Mod_64QAM);
scatterplot(Filt_64QAM_1,1,0,'*b');
Filt_64QAM_2 = filter(h2,1,Mod_64QAM);
scatterplot(Filt_64QAM_2,1,0,'*b');