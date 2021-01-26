%% Ejercicio 7
clear all
clc
close all
%% Ejercicio 7a
clc
clear all
close all
% Par�metros simulaci�n
N = 10000;
snr = 10.^1.5;
M = 6; %M es el n�m. de usuarios
Umbral = 10e-4;
% Inicializaci�n vectores
reg_usuarios = zeros(1,N);
Ax1 = modem.qammod('M',64,'SymbolOrder','gray', 'InputType', 'bit');
Ax2 = modem.qammod('M',16,'SymbolOrder','gray', 'InputType', 'bit');
Ax3 = modem.qammod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
nb_v = zeros(M,1); %Num. de bits tx por cada usuario

%Generamos info
btx = round(rand(M,2*N));

%Generamos canal
h = sqrt(.5)*(randn(M,N)+1i*randn(M,N));
bool_64qam_ac = zeros(N,1);
bool_16qam_ac = zeros(N,1); 
bool_4qam_ac = zeros(N,1);
bool_none_ac = zeros(N,1);

for n = 1:N
    %Identificamos el usuario que accede al canal
    m_opt = mod(n-1,M) + 1; % Acceso fijo (round robin)
    reg_usuarios(n) = m_opt;
    nb = nb_v(m_opt);
    
    BER1(n) = .2*exp(-snr*abs(h(m_opt,n))^2/(64-1));
    BER2(n) = .2*exp(-snr*abs(h(m_opt,n))^2/(16-1));
    BER3(n) = .2*exp(-snr*abs(h(m_opt,n))^2/(4-1));
    % Comprobamos con que modulaci�n va a transmitir.
    if BER1(n) < Umbral % 64-QAM
        nb_v(m_opt) = nb_v(m_opt) + 6;
        bool_64qam_ac(n) = 1;
    elseif BER2(n) < Umbral % 16-QAM
        nb_v(m_opt) = nb_v(m_opt) + 4;
        bool_16qam_ac(n) = 1;
    elseif BER3(n) < Umbral % 4-QAM
        nb_v(m_opt) = nb_v(m_opt) + 2;
        bool_4qam_ac(n) = 1;
    else % No Transmite
        nb_v(m_opt) = nb_v(m_opt);
        bool_none_ac(n) = 1;
    end    
end
for m_opt = 1:M
    disp( ['   Tasa usuario ', num2str(m_opt), ' = ',num2str(nb_v(m_opt)/N) ] )
end
no_transmitido = sum(bool_none_ac)/100
qam4 = sum(bool_4qam_ac)/100
qam16 = sum(bool_16qam_ac)/100
qam64 = sum(bool_64qam_ac)/100

%% Ejercicio 7b
clc
close all
clear all

% Par�metros simulaci�n
N = 10000;
snr = 10.^1.5;
M = 6; %M es el n�m. de usuarios
Umbral = 10e-4;
% Inicializaci�n vectores
y = zeros(1,N); 
reg_usuarios = zeros(1,N);
Ax1 = modem.qammod('M',64,'SymbolOrder','gray', 'InputType', 'bit');
Ax2 = modem.qammod('M',16,'SymbolOrder','gray', 'InputType', 'bit');
Ax3 = modem.qammod('M',4,'SymbolOrder','gray', 'InputType', 'bit'); 
nb_v = zeros(M,1); %Num. de bits tx por cada usuario

%Generamos info
btx = round(rand(M,2*N));
h = sqrt(.5)*(randn(M,N)+1i*randn(M,N));
time_us=zeros(M,1);
bool_64qam_ac = zeros(N,1);
bool_16qam_ac = zeros(N,1); 
bool_4qam_ac = zeros(N,1);
bool_none_ac = zeros(N,1);
Mod = [64 16 4];
for n = 1:N
    %Identificamos el usuario que accede al canal
    %calculo de la BER inst.
    for user = 1:6
        for mod = 1:3
            BER(user,mod) = .2*exp(-snr*abs(h(user,n))^2/(Mod(mod)-1));
        end
    end
    aux = BER < Umbral;
    if sum(aux(:,1)) == 1
        % Transmite con 64
            [m_opt] = find(aux(:,1) == 1);
            modul = 1;      
    elseif sum(aux(:,1)) >= 2
        % Se elige a uno para transmitir con 64
            [user] = find(aux(:,1) == 1);
            modul = 1;
            [max,pos] = max(time_us(user));
            m_opt = user(pos);
    elseif sum(aux(:,2)) == 1
        % Transmite con 16
            [m_opt] = find(aux(:,2) == 1);
            modul = 2;
    elseif sum(aux(:,2)) >= 2
        % Se elige a uno para transmitir con 16
            [user] = find(aux(:,2) == 1);
            modul = 2;
            [max,pos] = max(time_us(user));
            m_opt = user(pos);
    elseif sum(aux(:,3)) == 1
        % Se transmite con 4
            [m_opt] = find(aux(:,3) == 1);
            modul = 3;
    elseif sum(aux(:,3)) >= 2
        % Se elige a uno para transmitir con 4
            [user] = find(aux(:,3) == 1);
            modul = 3;
            [max,pos] = max(time_us(user));
            m_opt = user(pos);
    else
        modul = 0;
        m_opt = 0;
    end
    % Actualizo timer
    if m_opt == 0
        m_opt = 1;
        reg_usuarios(n) = m_opt;
    else
        time_us = time_us +1;
        time_us(m_opt) = 1;
        reg_usuarios(n) = m_opt;
        nb = nb_v(m_opt);       
    end
    clearvars user max pos
    
    if modul == 1
        nb_v(m_opt) = nb_v(m_opt) + 6;
        bool_64qam_ac(n) = 1;
    elseif modul == 2
        nb_v(m_opt) = nb_v(m_opt) + 4;
        bool_16qam_ac(n) = 1;
    elseif modul == 3 
        nb_v(m_opt) = nb_v(m_opt) + 2;
        bool_4qam_ac(n) = 1;
    else
        nb_v(m_opt) = nb_v(m_opt);
        bool_none_ac(n) = 1;
    end
    clearvars m_opt
end
for m_opt = 1:M
    disp( ['   Tasa usuario ', num2str(m_opt), ' = ',num2str(nb_v(m_opt)/N) ] )
end
no_transmitido = sum(bool_none_ac)/100
qam4 = sum(bool_4qam_ac)/100
qam16 = sum(bool_16qam_ac)/100
qam64 = sum(bool_64qam_ac)/100