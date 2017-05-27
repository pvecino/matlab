%% HITO 3: BOOTSTRAP

clear all; clc;close all;
X = [-2.4100 4.8600 6.0600 9.1100 10.2000 12.8100 13.1700 14.1000 15.7700 15.7900];
B = 1000;

for k = 1:1000
    muestreo_reemplazo = datasample(X,10);       %muestras bootstrap
    y(k) = mean(muestreo_reemplazo);  
    varianza_boot(k) = var(muestreo_reemplazo);
end

figure
subplot(2,1,1)
hist(y, 100)
title('Histograma media Bootstrap')
subplot(2,1,2)
hist(varianza_boot, 100)
title('Histograma varianza Bootstrap')

alpha = [0.1 0.05 0.01];
X_boot_ordenado = sort(y);
X_boot_ordenado_varianza = sort(varianza_boot);
for k = alpha
       primer_valor = (k*B*0.5);
       ultimo_valor = B - primer_valor;
       disp([ 'Intervalo de confianza (media)' num2str((100 - k*100)) '% : ' num2str(X_boot_ordenado(primer_valor)) '-' num2str(X_boot_ordenado(ultimo_valor))])
       disp([ 'Intervalo de confianza (varianza)' num2str((100 - k*100)) '% : ' num2str(X_boot_ordenado_varianza(primer_valor)) '-' num2str(X_boot_ordenado_varianza(ultimo_valor))])
end
