i = sqrt(-1);
for j = 1:(5000)
  uF = fft2(squeeze(u(j,:,:)))/(nM*nN);
  vF = fft2(squeeze(v(j,:,:)))/(nM*nN);
  hF = fft2(squeeze(h(j,:,:)))/(nM*nN);
  for k = 1:nM
    for l = 1:nN
      if(k <= nM/2)
        km = 2*pi*(k - 1);
      else
        km = 2*pi*(k - 1 - nM);
      end
      if(l <= nN/2)
        lm = 2*pi*(l - 1);
      else
        lm = 2*pi*(l - 1 - nN);
      end
      omega = sqrt(1 + km^2 + lm^2);
      aPV(j,k,l) = (i*km*vF(k,l) - i*lm*uF(k,l) - hF(k,l))/omega;
      if (km == 0 && lm == 0)
        aGP(j,k,l) = ( uF(k,l) + i*vF(k,l)) / sqrt(2);
        aGM(j,k,l) = (-uF(k,l) + i*vF(k,l)) / sqrt(2);
      else
        aGP(j,k,l) = (( omega*km - i*lm)*uF(k,l) + (i*km + omega*lm)*vF(k,l) + (omega^2 - 1)*hF(k,l))/(sqrt(2*(km^2 + lm^2))*omega);
        aGM(j,k,l) = ((-omega*km - i*lm)*uF(k,l) + (i*km - omega*lm)*vF(k,l) + (omega^2 - 1)*hF(k,l))/(sqrt(2*(km^2 + lm^2))*omega);
      end
    end
  end
end