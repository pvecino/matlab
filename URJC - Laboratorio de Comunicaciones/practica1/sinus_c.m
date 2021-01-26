function [s,t]= sinus_c(A, f0, phi, fs, ti, tf)
%A= amplitud
% f0 = frecuencia en Herzios
% phi = fase inicial en radianes
% fs = frecuencia de muestreo en Herzios
% ti = tiempo inicial en segundos
% tf = tiempo final en segundos
Ts = 1/fs;
t = ti:Ts:tf;
s = A*cos(2*pi*f0*t + phi);

end