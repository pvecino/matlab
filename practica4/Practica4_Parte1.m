%% Parte 2
clear all
close all
clc
load Muestras_Hito1.mat
Wo = 25; 
%% ESTIMADOR W0
%% CASO 1 k = 10
    k_10=10;
    x_10 = Observaciones1(:,1,:);
    y_10 = Observaciones1(:,2,:);
    % Estimador Wo 
        Wo_est_10 = sum(y_10-5 .* x_10)/k_10;
        sesgo_Wo_10 = mean(Wo_est_10) - Wo;
        var_Wo_10 = var(Wo_est_10);
%% CASO 2 k = 100
    k_100=100;
    x_100 = Observaciones2(:,1,:);
    y_100 = Observaciones2(:,2,:);
    % Estimador Wo 
        Wo_est_100 = sum(y_100-5 .* x_100)/k_100;
        sesgo_Wo_100 = mean(Wo_est_100) - Wo;
        var_Wo_100 = var(Wo_est_100);
%% CASO 3 k = 1000
    k_1000=1000;
    x_1000 = Observaciones3(:,1,:);
    y_1000 = Observaciones3(:,2,:);
    % Estimador Wo 
        Wo_est_1000 = sum(y_1000-5 .* x_1000)/k_1000;
        mean(Wo_est_1000)
        sesgo_Wo_1000 = mean(Wo_est_1000) - Wo;
        var_Wo_1000 = var(Wo_est_1000);
%% VALORES DE CONFIANZA
vc= [90, 95, 99]; % (%)
    % Ordenamos Wo
    Wo_est_10 = sort(Wo_est_10);
    Wo_est_100 = sort(Wo_est_100);
    Wo_est_1000 = sort(Wo_est_1000);
intervalos = [50,25,5];
conf_w0 = zeros(9,2); 
n = 1;
disp('Intervalos de confianza de Wo')
for k = 1:3
    if k ==1
        w0 = Wo_est_10;
        disp('Para k = 10')
    elseif k == 2
        w0 = Wo_est_100;
        disp('Para k = 100')
    else
        w0 = Wo_est_1000;
        disp('Para k = 1000')
    end
    r = 1;
    for m = intervalos
        w0_ML= w0(m+1:end-m);
        conf_w0(n,1) = w0_ML(1);
        conf_w0(n,2) = w0_ML(end);       
        disp([ '    Al ' num2str(vc(r)) '% : ' num2str(conf_w0(n,1)) '-' num2str(conf_w0(n,2))])
        n =n+1;
        r = r+1;
    end
end
