% demo and checks for CSD paper: here, conduction delays are modeled in
% terms of phase delays following the Fourier transform of the system's
% kernels.
%==========================================================================
clear; clf
phase = 1;
 
 
% frequencies of interest
%--------------------------------------------------------------------------
N  = 128;
dt = 1/N;
Hz = (1:N/2)';
 
W  = [0 1:N/2 -(N/2 - 1):-1]';
W  = conj(2*pi*sqrt(-1)*W);
t  = (1:N)'*dt;
d  = [0 31;
      4  0]/1000;
 
 
% kernels
%--------------------------------------------------------------------------
if phase
    K1(:,1,1)   = exp(-t*4).*sin(2*pi*8*t);
    K1(:,2,1)   = exp(-t*4).*sin(2*pi*8 *(t - d(2,1)))/4;
 
    K1(:,1,2)   = exp(-t*4).*sin(2*pi*32*(t - d(1,2)))/4;
    K1(:,2,2)   = exp(-t*4).*sin(2*pi*32*t)/2;
 
    % Transfer functions
    %----------------------------------------------------------------------
    S   = fft(K1);
 
    % [cross]-spectral density
    %----------------------------------------------------------------------
    [N,nc,nu] = size(K1);
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = 1:nc
 
            % cross-spectral density from neuronal interactions
            %--------------------------------------------------------------
            for k = 1:nu
                Gij      = S(:,i,k).*conj(S(:,j,k));
                Gij      = Gij([1:N/2] + 1);
                G(:,i,j) = G(:,i,j) + Gij;
            end
        end
    end
 
else
 
    K1(:,1,1)   = exp(-t*4).*sin(2*pi*8*t);
    K1(:,2,1)   = exp(-t*4).*sin(2*pi*8*t)/4;
 
    K1(:,1,2)   = exp(-t*4).*sin(2*pi*32*t)/4;
    K1(:,2,2)   = exp(-t*4).*sin(2*pi*32*t)/2;
 
    % Transfer functions
    %----------------------------------------------------------------------
    S   = fft(K1);
 
    % [cross]-spectral density
    %----------------------------------------------------------------------
    [N,nc,nu] = size(K1);
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = 1:nc
            
            % cross-spectral density with delays
            %--------------------------------------------------------------
            for k = 1:nu
                Si       = S(:,i,k).*exp(W*d(i,k));
                Sj       = S(:,j,k).*exp(W*d(j,k));
                Gij      = Si.*conj(Sj);
                Gij      = Gij([1:N/2] + 1);
                G(:,i,j) = G(:,i,j) + Gij;
            end
        end
    end
 
end
 
% plot cross spectrum densities
%----------------------------------------------------------------------
subplot(2,2,1)
plot(Hz,real(G(:,1,1)),Hz,real(G(:,1,2)),Hz,real(G(:,2,2)))
title('real csd')
xlabel('Hz')
axis square
 
subplot(2,2,2)
plot(Hz,imag(G(:,1,1)),Hz,imag(G(:,1,2)),Hz,imag(G(:,2,2)))
title('imag csd')
xlabel('Hz')
axis square
drawnow
 
return
 
 
% A simple test of the signs of imaginary cross-section density
%==========================================================================
clear Y
dt      = 1/256;
Hz      = 1:64;
T       = spm_vec(0:dt:2);
d       = 8;
Y(:,1)  = sin(2*pi*16*T) + randn(length(T),1)/256;
Y(:,2)  = sin(2*pi*16*(T - d/1000)) + randn(length(T),1)/256;
 
 
% VAR estimate
%--------------------------------------------------------------------------
mar     = spm_mar(Y,12);
mar     = spm_mar_spectra(mar,Hz,1/dt);
csd     = mar.P;
 
% Equivalent representations
%--------------------------------------------------------------------------
[ccf pst] = spm_csd2ccf(csd,Hz);
[coh fsd] = spm_csd2coh(csd,Hz);
 
% graphical results
%--------------------------------------------------------------------------
subplot(4,1,1)
plot(T(1:64),Y(1:64,:))
title('predicted')
xlabel('time series')
legend('first','second')
 
subplot(4,2,3)
plot(Hz,real(csd(:,2,1)))
title('real csd')
xlabel('Hz')
 
subplot(4,2,4)
plot(Hz,imag(csd(:,2,1)),'-.')
title('imag csd')
xlabel('Hz')
 
subplot(4,2,6)
plot(Hz,1000*fsd(:,2,1))
title('predicted')
xlabel('Hz')
