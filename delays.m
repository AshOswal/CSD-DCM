% simulation of delays taken from Friston 2011 - DCM for CSD paper

clear all;close all;

N = 1e5;

% lambda is spectral density of innovation
lambda(1,:) = [exp(0) exp(0.125) exp(0.25) exp(0.5) exp(1) exp(1.5) exp(2) exp(2.5) exp(3) exp(3.5) exp(4) exp(4.5) exp(5) exp(5.5) exp(6) exp(6.5) exp(7) exp(7.5) exp(8)];
lambda(2,:) = [exp(0) exp(0.125) exp(0.25) exp(0.5) exp(1) exp(1.5) exp(2) exp(2.5) exp(3) exp(3.5) exp(4) exp(4.5) exp(5) exp(5.5) exp(6) exp(6.5) exp(7) exp(7.5) exp(8)];
lambda      = sqrt(lambda); 

% transfer function
K(1,1)  = 2;
K(1,2)  = exp(1i*(pi/4));
K(2,1)  = 2*exp(1i*(pi/8));
K(2,2)  = 2;

S       = zeros(2,N,2,size(lambda,2));

for LL = 1:2
    for  L = 1:size(lambda,2)
        if LL == 1
            I_lambda = [lambda(1,1), lambda(1,1);lambda(2,L), lambda(2,L)];
        else
            I_lambda = [lambda(1,L), lambda(1,L);lambda(2,1), lambda(2,1)];
        end
        
        for n = 1:N
            Uk = sum(randn(2,2).*I_lambda.*[1,1i;1,1i],2);
            for s = 1:2
                for in = 1:size(Uk,1)
                    S(s,n,LL,L) = S(s,n,LL,L) + (abs(K(in,s)) * abs(Uk(in,:)) * exp(1i*(angle(Uk(in,:)) + angle(K(in,s)))));
                end
            end
        end
        
    end
end
%%
% we can get the phase difference from the 
S(2,:) = ctranspose(S(2,:));
delta_phi = squeeze(angle(S(1,:,:,:) .* (S(2,:,:,:))));
%%
figure;
iv = 1:4:19;
for L = 1:size(delta_phi,2)
    for LL = 1:numel(iv)%size(delta_phi,3)
        d_ij = delta_phi(:,L,iv(LL));
        subplot(2,numel(iv),LL+numel(iv)*(L-1));
        histogram(d_ij,100,'Normalization','probability','FaceColor','k');xlim([-1 1])
        xlabel('phase difference dij');
        ylabel('density');
        if L == 1
            title(['log lambda1/lambda2 =' num2str(log(1/(lambda(2,iv(LL)).^2)))]);
        else
            title(['log lambda1/lambda2 =' num2str(log((lambda(2,iv(LL)).^2))/1)]);
        end
        
    end
end
%%

% Generate the second plot
E_mean_csd = angle(mean(squeeze(S(1,:,:,:) .* (S(2,:,:,:))),1)); % phase delay of mean
E_dij      = squeeze(mean(delta_phi,1)); % mean phase - delay
figure;plot(E_dij(:),E_mean_csd(:),'ks');hold on;
plot([-0.8,0.4],[-0.8 0.4],'r-'); 
xlabel('<dij> mean of phase difference');
ylabel('arg<gij> phase difference of mean');

%S = unwrap(angle(S),[],1);
%delta_phi2                =   S(1,:) - S(2,:);
% delta_phi(delta_phi>0)   = mod(delta_phi(delta_phi>0) ,pi);
% delta_phi(delta_phi<0)   = mod(delta_phi(delta_phi<0) ,-pi);

figure; histogram(delta_phi,500,'Normalization','probability')
% xlabel('phase difference dij');
% ylabel('density')
%figure; histogram(delta_phi2,400)