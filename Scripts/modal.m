function [q, gp, gm] = modal(u, v, h, Lsize)

[NumModes, dummy, NumT] = size(u);

Lx = Lsize;
Ly = Lsize;
nM = NumModes;
nN = NumModes;
nT = NumT;

q = zeros(nM, nN, nT);
p = zeros(nM, nN, nT);
m = zeros(nM, nN, nT);

for k = 1:nM
  kmode(k) = 2*pi * ((k-1) - nM*floor((k-1)/(nM/2 + 1))) / Lx;
end
  

for t = 1:nT
    uF = fft2(u(:,:,t));
    vF = fft2(v(:,:,t));
    hF = fft2(h(:,:,t));
    for l = 1:nN
        for k = 1:nM
            kmode = 2*pi * ((k-1) - nM*floor((k-1)/(nM/2 + 1))) / Lx;
            lmode = 2*pi * ((l-1) - nN*floor((l-1)/(nN/2 + 1))) / Ly;
      
            Km = sqrt(kmode^2 + lmode^2);
            omega = sqrt(1 + Km^2);

            q(k,l,t) = (-i*lmode*uF(k,l) + i*kmode*vF(k,l) - hF(k,l))...
                / omega / (nM*nN);

            if(k == 1 && l == 1)
                gp(k,l,t) = (uF(1,1) + i*vF(1,1)) / sqrt(2) / (nM*nN);
            else
                gp(k,l,t) = ((omega*kmode - i*lmode)/Km*uF(k,l)...
                    + (i*kmode + omega*lmode)/Km*vF(k,l)...
                    + Km*hF(k,l))/(sqrt(2)*omega) / (nM*nN);
            end

            if(k == 1 && l == 1)
                gm(k,l,t) = (-uF(1,1) + i*vF(1,1)) / sqrt(2) / (nM*nN);
            else
                gm(k,l,t) = ((-omega*kmode - i*lmode)/Km*uF(k,l)...
                    + (i*kmode - omega*lmode)/Km*vF(k,l)...
                    + Km*hF(k,l))/(sqrt(2)*omega) / (nM*nN);
            end
      end
    end
end