function [x, n]=dsf_sintesis(a, N)

n = 0:N-1;
x = zeros(1,length(N));

for  k = n
    
    x(k + 1) = sum(a.*exp(1j*k*2*pi*n/N));
end

end