function g = gamma(k1, l1, s1, k2, l2, s2, k3, l3, s3)

o1 = s1*sqrt(1+k1^2+l1^2); o2 = s2*sqrt(1+k2^2+l2^2); o3 = s3*sqrt(1+k3^2+l3^2);

N1 = sqrt((k1^2 + l1^2) * (o1^2 + 1) + (o1^2 - 1)^2);
N2 = sqrt((k2^2 + l2^2) * (o2^2 + 1) + (o2^2 - 1)^2);
N3 = sqrt((k3^2 + l3^2) * (o3^2 + 1) + (o3^2 - 1)^2);

k1ck2 = k1*l2 - k2*l1;
k1ck3 = k1*l3 - k3*l1;
k2ck3 = k2*l3 - k3*l2;
k1dk2 = k1*k2 + l1*l2;
k1dk3 = k1*k3 + l1*l3;
k2dk3 = k2*k3 + l2*l3;

G1 = (o1^2 - 1)*(o3^2 - 1)*(k1ck2 - i*k1dk2*o2) + (o1^2 - 1)*(o2^2 - 1)*(k1ck3 - i*k1dk3*o3);
G2 = (k2ck3 - i*k2dk3*o3)*((1+o1*o2)*k1dk2 + i*(o1+o2)*k1ck2) + (-k2ck3 - i*k2dk3*o2)*((1+o1*o3)*k1dk3 + i*(o1+o3)*k1ck3);

g = (G1 + G2)/(2*N1*N2*N3);