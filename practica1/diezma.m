function [y, n] = diezma(x, N)

    l=0:length(x)-1;
    y=x(mod(l,N)==0);
    n = 0:length(y)-1;
end