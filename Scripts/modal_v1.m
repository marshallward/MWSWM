function [q, gp, gm] = modal(u, v, h, Lsize, NumModes, NumT)

Lx = Lsize;
Ly = Lsize;
nM = NumModes;
nN = NumModes;
nT = NumT;

for t = 1:nT
    uF = fft2(squeeze(u(t,:,:)));
    vF = fft2(squeeze(v(t,:,:)));
    hF = fft2(squeeze(h(t,:,:)));
    for k = 1:nM
        for l = 1:nN
            lmode = 2*pi * ((k-1) - nM*floor((k-1)/(nM/2 + 1))) / Lx;
            kmode = 2*pi * ((l-1) - nN*floor((l-1)/(nN/2 + 1))) / Ly;
      
            Km = sqrt(kmode^2 + lmode^2);
            omega = sqrt(1 + Km^2);

            q(t,k,l) = (-i*lmode*uF(k,l) + i*kmode*vF(k,l) - hF(k,l))...
                / omega / (nM*nN);

            if(k == 1 && l == 1)
                gp(t,k,l) = (uF(1,1) + i*vF(1,1)) / sqrt(2) / (nM*nN);
            else
                gp(t,k,l) = ((omega*kmode - i*lmode)/Km*uF(k,l)...
                    + (i*kmode + omega*lmode)/Km*vF(k,l)...
                    + Km*hF(k,l))/(sqrt(2)*omega) / (nM*nN);
            end

            if(k == 1 && l == 1)
                gm(t,k,l) = (-uF(1,1) + i*vF(1,1)) / sqrt(2) / (nM*nN);
            else
                gm(t,k,l) = ((-omega*kmode - i*lmode)/Km*uF(k,l)...
                    + (i*kmode - omega*lmode)/Km*vF(k,l)...
                    + Km*hF(k,l))/(sqrt(2)*omega) / (nM*nN);
            end
      end
    end
end