function y = mifuncion(x, b, opcion)

% x es un vector
% b es una constante
%opcion es un numero entero del 1 al 5
if nargin == 2
        opcion = 1;
end

switch opcion
    case 1
% en caso de que opcion sea 1, y sera 1 si algun elemento de x
% es igual a b.
        for n = 1:length(x)
            if x(n) == b
                y = 1;
            else 
                y = -1;
            end
        end
    case 2
% en caso de que opcion sea 2, y sera un vector con los elementos
%de x que sean mayores que b
        m = 1;
        for n = 1:length(x)
            if x(n) > b
                y(m) = x(n);
                m = m + 1;
            else 
                y = -1;
            end
        end
    case 3
 % en caso de que opcion sea 3, y sera un vector igual a x, pero 
 % con los elementos menores a b, igual a 0.
        for n = 1:length(x)
            if x(n) < b
                y(n)= 0;
            else
                y(n) = x(n);
            end
        end
    case 4
 % en caso de que opcion sea 4, y es el vector x, pero ordenado
 % de menos a mayor.
        y = sort(x);
    case 5
  % en caso de que opcion sea 5, y debe ser una matriz formada por 
  % b columnas y que esas columnas sea el vector x
      aux=size(x);
            if aux(1) > 1
                x=x';
            end
            
            y=x';
            for n=2:b
                y = [y,x'];
            end
    otherwise
            y = -1;
        
end