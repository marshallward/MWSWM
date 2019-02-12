function gammaplot(Kmax,N)
[k,th] = meshgrid((0:1/N:1)*Kmax, (-1:1/(N/2):1)*pi);
[x,y] = pol2cart(th,k);
for i = 1:(N+1)
    for j = 1:(N+1)
        g(i,j) = gamma2(k(i,i),th(j,j));
    end
end

figure(1);
hold off;
polar([0 2*pi], [0 Kmax]);
title(strcat('\Gamma^{\pm\pm0} up to K = ',num2str(Kmax)));
hold on;
[C,h] = contourf(x,y,abs(g'));
colorbar;