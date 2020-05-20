clear all
close all

% specify model
%==========================================================================

% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
n     = 2; % number of sources
Nc    = 2; % number of channels

M.dipfit.type = 'LFP';
M.dipfit.model(1).source='MMC';
M.dipfit.model(2).source='BGT';
M.dipfit.Ns = n;
M.dipfit.Nc = Nc;

% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = [0 1;0 0];                         % a forward connection
A{2}  = [0 0;1 0];                         % a backward connection
% A{1}  = [0 0;0 0];                         % a forward connection
% A{2}  = [0 0;0 0];
A{3}  = sparse(n,n);                       % lateral connections
B     = {};                                % trial-specific modulation
C     = speye(n,n);                        % sources receiving innovations

% get priors
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(A,B,C,M.dipfit.model); % neuronal priors
[pE,pC] = spm_ssr_priors(pE,pC);           % spectral priors
[pE,pC] = spm_L_priors(M.dipfit,pE,pC);    % spatial  priors

% take priors from posterio inversion
% load('bestmodels_nobadtrials/DCM_OFF0_ON1_full_base14_priorsONOFF3_c_noDSnoJAr_nobadtrials_rescaled_he16_hc4_deepmiddlemmc_indnoise_mmc512_bgc512_it5_nobadtrials_rescaled_OFFON_DF_left');
% load('bestmodels_nobadtrials/DCM_OFFonly_optimizepriors_14_priorsONOFF3_he16_hc4_deepmiddlemmc_indnoise_mmc512_bgc512_c_noDSnoJAr_nobadtrials_rescaled_OFFON_DF_left');
% pE=DCM.Ep;
% pC=DCM.Cp;

% pE.A{1}=-32*ones(2,2);
% pE.A{2}=-32*ones(2,2);
% pE.A{3}=-32*ones(2,2);
% pE.A{4}=-32*ones(2,2);

% Suppress channel noise
%--------------------------------------------------------------------------
pE.b  = pE.b - 16;
pE.c  = pE.c - 16;

% create LFP model
%--------------------------------------------------------------------------

[x,f]    = spm_dcm_x_neural(pE,M.dipfit.model);

M.IS = 'spm_csd_mtf';
%M.FS = 'spm_fs_csd';
M.g  = 'spm_gx_erp';
M.f  = f;
M.x  = x;
M.n  = length(spm_vec(x));
M.pE = pE;
M.pC = pC;
M.m  = n;
M.l  = Nc;
M.Hz = (1:90);
M.N  = []; % order of delay approximation
ord  = {0 256 'no delay'};
%M.u  = sparse(M.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('priors_gransaverage_OFF.mat');

M.MMC_T     = MMC_T;
M.MMC_G     = MMC_G;
%M.BGC_T     = BGC_T;
%M.BGC_G     = BGC_G;
%M.BGC_E     = BGC_E;
M.MMC_E     = MMC_E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = {'D','T','S','G'};

% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


D     = 4;

%USE obs=1 FOR CORTEX AND obs=2 FOR STN
cortex = 1;
stn    = 2;

cols   = {'r','b','k','y','g'};

% determine how extrinsic delays influence the linearly approximated
% spectra as these are the major source of nonlinearity


for k = 1%:length(P)
    
    % check parameter exists
    %----------------------------------------------------------------------
    % sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
    %         spm_figure('GetWin',sfig);
    
    Q = getfield(pE,P{k});
    
    
    
    % line search (with solution for steady state)
    %----------------------------------------------------------
    if strcmp(P{k},'D')
        for ii = 1:numel(Q)
            dQ(ii,:) = linspace(Q(ii) - D, Q(ii) + D,4);
        end
    end
    dQ = [dQ(:,1:2) zeros(size(dQ,1),1) dQ(:,3:4)];
    leg = {};
    for l = 1:size(dQ,2)
        leg = [leg, ['exp ' num2str(dQ(1,l))]];
    end
    
    for o = 1:numel(ord)
        M.N = ord{o};
        figure(2*o-1);
        figure(2*o)
        for q = 3% 1:length(dQ)
            qE     = pE;
            qE     = setfield(qE,P{k},reshape(dQ(:,q),size(Q)));
            
            
            %[y,w,s,g,K,tt] = spm_csd_mtf_ash(qE,M,[]);
            
            
            % integrate with no inputs using local linearisation
            %--------------------------------------------------------------------------
            N    = 1001;
            U.dt = 1/N;
            U.u  = [zeros(2,20)';1/U.dt 0 ;zeros(2,N)'];
            LFP  = spm_int_L(qE,M,U);
            U.u  = [zeros(2,20)';0 1/U.dt ;zeros(2,N)'];
            LFP1 = spm_int_L(qE,M,U);
            
            figure(2*o-1);

            subplot(2,length(dQ),q);
            plot(U.dt*(1:N+21),U.u,'k:');
            xlabel('time(s)')
            ylabel('Amplitude');
            xlim([0 0.2]);
            set(gca,'FontSize',12);

            subplot(2,length(dQ),q+length(dQ));
            plot((1:N+21).*U.dt,LFP(:,1),'r');hold on;
            plot((1:N+21).*U.dt,LFP(:,2),'b:');hold on;
            plot((1:N+21).*U.dt,LFP1(:,1),'r:');hold on;
            plot((1:N+21).*U.dt,LFP1(:,2),'b');hold on;
            xlabel('time(s)')
            ylabel('Amplitude');
            title(['IR delay exp' num2str(dQ(1,q))]);
            xlim([0 0.2]);
            set(gca,'FontSize',12);
            if isequal(q,length(dQ))
                l = legend({'M1 1','STN 2','M1 2', 'STN 1'});
                set(l,'box','off');               
            end

            
            %---------extract cross (auto) spectra
            G     = zeros(numel(w),Nc,Nc);
            for i = 1:numel(w)
                Si       = shiftdim(s{:}(i,:,:),1);
                G(i,:,:) = Si*Si';
            end
            
            % extract auto(cross) spectra
            GW_stn(:,q) = squeeze(G(:,stn,stn));
            GW_m1(:,q)  = squeeze(G(:,cortex,cortex));
            GW_csd(:,q) = squeeze(G(:,cortex,stn));
            
            % extract time domain kernels
            K_stn       = squeeze(K(:,2,:));
            K_m1        = squeeze(K(:,1,:));
            
            ind = 4*length(dQ);
            is  = 4*length(dQ)/4;
            
            figure(2*o);
            subplot(2,ind,1:is-1);
            plot(w,abs(GW_stn(:,q)),cols{q});hold on;
            xlabel('frequency Hz');title('STN power');set(gca,'FontSize',12);
            
            subplot(2,ind,is+1:is*2-1);
            plot(w,abs(GW_m1(:,q)),cols{q});hold on;
            xlabel('frequency Hz');title('M1 power');set(gca,'FontSize',12);
            
            subplot(2,ind,2*is+1:is*3-1);
            plot(w,abs(GW_csd(:,q)),cols{q});hold on;
            xlabel('frequency Hz');title('abs CSD');set(gca,'FontSize',12);
            
            subplot(2,ind,is*3+1:is*4-1);
            plot(w,imag(GW_csd(:,q)),cols{q});hold on;
            xlabel('frequency Hz');title('imag CSD');set(gca,'FontSize',12);
            h = legend(leg);
            set(h,'box','off');
            
            slim = is*4+1+ 4*(q-1);
            elim = is*4+3+ 4*(q-1);
            subplot(2,20,slim:elim);
            plot(tt,K_stn(:,1),'b:',tt,K_stn(:,2),'b');hold on;
            plot(tt,K_m1(:,1),'r',tt,K_m1(:,2),'r:');xlim([0 0.2]);
            if isequal(q,length(dQ))
                h1 = legend({'STN 2','STN 1','M1 1', 'M1 2'});
                set(h1,'box','off');
            end
            set(gca,'FontSize',12);
            xlabel('time');
            title(['kernel exp' num2str(dQ(1,q))]);
            
            
        end
    end
    
    
end



