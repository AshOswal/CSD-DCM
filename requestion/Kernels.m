% High-order approximation to delay operator
%==========================================================================
clear; clf

% frequencies of interest
%--------------------------------------------------------------------------
N  = 1024;
dt = 1/N;
Hz = (1:N/2)';

W  = [0 1:N/2 -(N/2 - 1):-1]';
W  = conj(2*pi*sqrt(-1)*W);
t  = (1:N)'*dt;

% number of channels and exogenous (neuronal) inputs
%--------------------------------------------------------------------------
nc   = 3;  % number of channels
nu   = 3;  % number of inputs

hz   = 2*pi*64;
P.A  = [-8 hz  1;
      -hz -8  2;
       8  0 -8];


P.C  = [1 0 0;
       0 1 0;
       0 0 1];

% P.A  = [-10 50;
%         -50 -10];
% 
% P.C  = [1 0;
%         0 0];

M.f  = inline('P.A*x + P.C*u','x','u','P','M'); % the state equation 
M.g  = inline('x','x','u','P','M');             % the observer function which in this case is just 1
M.x  = sparse(nu,1);
M.pE = P;
M.n  = nc;
M.m  = nu;
M.l  = nc;
M.nodelay = 0;

%P.x  = M.x;

% get (extrinsic) delays and (intrinsic) delay operator
%--------------------------------------------------------------------------
% tj   = find(t < 64/1000);
% hj   = find(Hz < 128);
% T    = t(tj)*1000;

tj = find(t);
hj = find(Hz);
T  = t;


d    =  1:200:800;
D    = [1 4 4;
        4 1 4;
        4 4 1]/1000;

%D    = [1 1;...
%        1 1]/1e3;

for id = 1:length(d)
    
    % delay operator
    %----------------------------------------------------------------------
    P.D       = d(id)*D;
    P.m       = M.m;
    Q         = spm_dcm_delay(P,M);
    M.D       = Q;
    % augment and bi-linearise (with intrinsic delays)
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,P);
    
    % kernels
    %----------------------------------------------------------------------
    [K0,K1,K2]   = spm_kernels(M0,M1,L,N,dt);
    
    % plot kernels
    %----------------------------------------------------------------------
    subplot(2,1,1)
    plot(T,K1(tj,2,1)),hold on
    title('kernel')
    xlabel('time')
    axis square
    
    
    % Transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S   = fft(K1);
    
    
    % [cross]-spectral density from neuronal innovations
    %----------------------------------------------------------------------
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = 1:nc
            for k = 1:nu
                Gij      = S(:,i,k).*conj(S(:,j,k));
                Gij      = Gij((1:N/2) + 1);
                G(:,i,j) = G(:,i,j) + Gij;
            end
        end
    end
    
    
    % plot cross spectrum densities
    %----------------------------------------------------------------------
    subplot(2,2,3)
    plot(Hz(hj),abs(G(hj,2,1))),hold on
    title('real csd')
    xlabel('Hz')
    axis square
    
    subplot(2,2,4)
    plot(Hz(hj),imag(G(hj,2,1))),hold on
    title('imag csd')
    xlabel('Hz')
    axis square
    drawnow
    
    %pause
    
end
