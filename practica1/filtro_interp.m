function [h n] = filtro_interp(L)
    k = -2*L : 2*L;
    h = sinc(k/L)
    n = 0:length(h)-1;
end