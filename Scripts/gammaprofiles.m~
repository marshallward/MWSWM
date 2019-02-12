clear k a G000 G00p G0pp Gp00 Gp0p Gppp;

K1 = 1/3;   L1 = 1;
K2 = 1;     L2 = 0;
K3 = -2/3;  L3 = 1;

a_min = -10; a_max = 10;
a = [a_min:1:a_max];
base = 2;

for j = 1:length(a)
k(j) = base^a(j);
end

kL = k(1:round(length(a)/2)-1);
kH = k(round(length(a)/2)+1:length(a));

for j = 1:length(a)
    G000(j) = gamma(k(j)*K1, k(j)*L1, 0, k(j)*K2, k(j)*L2, 0, k(j)*K3, k(j)*L3, 0);
    G00p(j) = gamma(k(j)*K1, k(j)*L1, 0, k(j)*K2, k(j)*L2, 0, k(j)*K3, k(j)*L3, 1);
%    G0pp(j) = gamma(k(j)*K1, k(j)*L1, 0, k(j)*K2, k(j)*L2, 1, k(j)*K3, k(j)*L3, 1);
    Gp00(j) = gamma(k(j)*K1, k(j)*L1, 1, k(j)*K2, k(j)*L2, 0, k(j)*K3, k(j)*L3, 0);
    Gp0p(j) = gamma(k(j)*K1, k(j)*L1, 1, k(j)*K2, k(j)*L2, 0, k(j)*K3, k(j)*L3, 1);
    Gppp(j) = gamma(k(j)*K1, k(j)*L1, 1, k(j)*K2, k(j)*L2, 1, k(j)*K3, k(j)*L3, 1);
end

figure(1); hold off;
loglog(k,abs(G000),'k');
hold on;
plot(kL,kL.^4,'k-.');
plot(kH,0.5*kH.^1,'k-.');
title('|\Gamma^{000}| vs K');
xlabel('Mode, K');
ylabel('Interaction Coefficient, |\Gamma^{000}|');


figure(2); hold off;
loglog(k,abs(G00p),'k');
hold on;
plot(kL,0.6*kL.^1,'k-.');
plot(kH,0.4*kH.^1,'k-.');
title('|\Gamma^{00+}| vs K');
xlabel('Mode, K');
ylabel('Interaction Coefficient, |\Gamma^{00+}|');


figure(3); hold off;
loglog(k,abs(Gp00),'k');
hold on;
plot(kL,2*kL.^3,'k-.');
plot(kH,1.2*kH.^1,'k-.');
title('|\Gamma^{+00}| vs K');
xlabel('Mode, K');
ylabel('Interaction Coefficient, |\Gamma^{+00}|');


figure(4); hold off;
loglog(k,abs(Gp0p),'k');
hold on;
plot(kL,1.1*kL.^2,'k-.');
plot(kH,0.6*kH.^1,'k-.');
title('|\Gamma^{++0}| vs K');
xlabel('Mode, K');
ylabel('Interaction Coefficient, |\Gamma^{++0}|');


figure(5); hold off;
loglog(k,abs(Gppp),'k');
hold on;
plot(kL,0.6*kL,'k-.');
plot(kH,0.12*kH,'k-.');
title('|\Gamma^{+++}| vs K');
xlabel('Mode, K');
ylabel('Interaction Coefficient, |\Gamma^{+++}|');


%figure(2); loglog(k,abs(G00p),'.');
%figure(3); loglog(k,abs(G0pp),'.');
%figure(4); loglog(k,abs(Gp00),'.');
%figure(5); loglog(k,abs(Gp0p),'.');
%figure(6); loglog(k,abs(Gppp),'.');