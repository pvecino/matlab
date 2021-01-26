


function [y n] = expande(x, L)

    k = 0:length(x)-1;
    y = zeros(1, L*length(x));
    y(k*L+1) = x(k+1);
    n = 0:length(y)-1;
end