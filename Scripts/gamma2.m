function g = gamma2(k,th)

om = sqrt(1+k^2);
g = k^2/2/(2*om^2)*((1+om^2)*sin(2*th) - sin(th) + i*om*(2*cos(2*th) - cos(th)))/sqrt(4*k^2*sin(th/2)^2 + 1);