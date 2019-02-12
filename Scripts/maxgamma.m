function [k,thmax] = maxgamma(Kmin,Kmax,N)

theta = (-1:2/N:1)*pi;
k = (0:1/N:1)*(Kmax-Kmin) + Kmin;

for i = 1:(N+1)
    for j = 1:(N+1)
        g(j) = gamma2(k(i),theta(j));
        [maxg, maxj] = max(abs(g));
    end
    thmax(i) = theta(maxj);
end