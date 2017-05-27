function [x, n] = sinus(A, omega, fase, n0, n1)
% Funcion seno que genrena una sinusoide de longitud finita
% Recibe 5 parametros: amplitud, omega0, fase, y dos numeros que
% especifiquen el primer y el ultimo valor de n para los que se genera la
% sinusoide.
n= n0:0.000001:n1;
x=A*cos(omega*n+fase);
end