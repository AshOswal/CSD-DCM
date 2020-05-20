% demo and checks for CSP paper
%==========================================================================

% simple test
%--------------------------------------------------------------------------
clear Y
V.dt    = 0.002;
M.Hz    = 1:64;
T       = spm_vec(0:V.dt:8);
d       = 6;
D       = fix(d/1000/V.dt);
Y(:,1)  = sin(2*pi*16*T) + randn(length(T),1)/8;
Y(:,2)  = sin(2*pi*16*(T - 0.006)) + randn(length(T),1)/8;
Y(:,1)  = randn(length(T),1)*2;
Y(:,2)  = [Y((D + 1):end,1); randn(D,1)] + randn(length(T),1)/1024;

 
% VAR estimate
%--------------------------------------------------------------------------
mar     = spm_mar(Y,12);
mar     = spm_mar_spectra(mar,M.Hz,1/V.dt);
CSD     = mar.P;
csd     = CSD;
 
% Equivalent representations
%--------------------------------------------------------------------------
[ccf pst] = spm_csd2ccf(csd,M.Hz);
[CCF PST] = spm_csd2ccf(CSD,M.Hz);

[coh fsd] = spm_csd2coh(csd,M.Hz);
[COH FSD] = spm_csd2coh(CSD,M.Hz);
 
i   = 2;
j   = 1;
t   = 32;
 
% graphical results
%--------------------------------------------------------------------------
figure;
subplot(4,2,1)
plot(1000*T(1:16),Y(1:16,:))
title('predicted')
xlabel('time series')


% subplot(4,2,1)
% plot(M.Hz,real(csd(:,j,i)),M.Hz,imag(csd(:,j,i)),'-.')
% title('predicted')
% xlabel('cross-spectra')


subplot(4,2,2)
plot(M.Hz,real(CSD(:,j,i)),M.Hz,imag(CSD(:,j,i)),':')
title('estimated')
xlabel('cross-spectra')
 
subplot(4,2,3)
plot(1000*pst,ccf(:,j,i))
title('predicted')
set(gca,'XLim',[-t t])
xlabel('cross-correlation')
 
subplot(4,2,4)
plot(1000*PST,CCF(:,j,i))
title('estimated')
set(gca,'XLim',[-t t])
xlabel('cross-correlation')

subplot(4,2,5)
plot(M.Hz,coh(:,i,j))
title('predicted')
xlabel('coherence')
 
subplot(4,2,6)
plot(M.Hz,COH(:,i,j))
title('estimated')
xlabel('coherence')

subplot(4,2,7)
plot(M.Hz,1000*fsd(:,i,j))
title('predicted')
xlabel('delay')

subplot(4,2,8)
plot(M.Hz,1000*fsd(:,i,j))
title('estimated')
xlabel('delay')
 
