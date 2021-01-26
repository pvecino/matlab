function d = delta(n)
%%%%%%%%%%%%%%%%%%%%%
%Subfuncion que implementa la Delta de Kronecker
%%%%%%%%%%%%%%%%%%%%%
 
%Inicializamos la secuencia de salida "d" como ceros
d = zeros(length(n),1);
 
%Buscamos el punto donde el indice es cero (delta no nula) y le asignamos el varor 1
d(n == 0) = 1;
