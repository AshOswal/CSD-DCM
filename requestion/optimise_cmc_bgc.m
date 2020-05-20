clear all, close all

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
M.FS = 'spm_fs_csd';
M.g  = 'spm_gx_erp';
M.f  = f;
M.x  = x;
M.n  = length(spm_vec(x));
M.pE = pE;
M.pC = pC;
M.m  = n;
M.l  = Nc;
M.Hz = (1:90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('priors_gransaverage_OFF.mat');

M.MMC_T     = MMC_T;
M.MMC_G     = MMC_G;
% M.BGC_T     = BGC_T;
% M.BGC_G     = BGC_G;
% M.BGC_E     = BGC_E;
M.MMC_E     = MMC_E;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = {'G','T','S','D'};

% solve for steady state
%--------------------------------------------------------------------------
M.x     = spm_dcm_neural_x(pE,M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iplot = 1;
ifig  = 1;
D     = 2;

%USE obs=1 FOR CORTEX AND obs=2 FOR STN
obs=2;


for p = 1:length(M.dipfit.model);
    
    for k = 1:length(P)
        
        % check parameter exists
        %----------------------------------------------------------------------
        sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
        %         spm_figure('GetWin',sfig);
        figure, set(gcf,'position',[396    12   460   664])
        
        if isempty(strmatch(P{k},'D'));
            
            Q = getfield(pE.int{p},P{k});
            
            if isnumeric(Q)
                for i = 1:size(Q,1)
                    for j = 1:size(Q,2);
                        
                        if iplot > 4
                            iplot = 1;
                            ifig  = ifig + 1;
                            sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
                            %                             spm_figure('GetWin',sfig);
                            figure, set(gcf,'position',[396    12   460   664])
                        end
                        
                        % line search (with solution for steady state)
                        %----------------------------------------------------------
%                         dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);
                                                dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);%6
                        for q = 1:length(dQ)
                            if dQ(q)<0; cl='b'; else cl='r'; end
                            tmp     = pE.int{p};
                            tmp     = setfield(tmp,P{k},{i,j},dQ(q));
                            qE      = pE; qE.int{p}=tmp;
                            [G w]   = spm_csd_mtf(qE,M,[]);
                            GW(:,q) = squeeze(G{1}(:,obs,obs));
                            
                            subplot(4,2,2*iplot - 1), hold on
                            plot(w,abs(GW(:,q)),cl)
                            
                        end
                        
                        % plot
                        %----------------------------------------------------------
                        xlabel('Frequency {Hz}')
                        title(sprintf('%s  Param: %s(%i,%i)',M.dipfit.model(p).source,P{k},i,j),'FontSize',12)
                        set(gca,'XLim',[0 w(end)]);
                        set(gca,'xtick',[0:10:w(end)])
                        set(gcf,'PaperPositionMode','auto')
                        
                        subplot(4,2,2*iplot - 0)
                        imagesc(dQ,w,log(abs(GW)))
                        title('Transfer functions','FontSize',12)
                        ylabel('Frequency')
                        xlabel('(log) parameter scaling')
                        axis xy; drawnow; colormap jet
                        
                        % update graphics
                        %----------------------------------------------------------
                        iplot     = iplot + 1;
                        
                    end
                end
            end
            
        elseif strmatch(P{k},'D') & strmatch(M.dipfit.model(p).source,'BGC');
            
            Q = getfield(pE,P{k});
            
            for i = 1:size(Q,1)
                for j = 1:size(Q,2);
                    
                    if iplot > 4
                        iplot = 1;
                        ifig  = ifig + 1;
                        sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
                        %                         spm_figure('GetWin',sfig);
                        figure, set(gcf,'position',[396    12   460   664])
                    end
                    
                    % line search (with solution for steady state)
                    %----------------------------------------------------------
%                     dQ    = linspace(Q(i,j) - D,Q(i,j) + D,32);
                                        dQ    = linspace(Q(i,j) - D,Q(i,j) + D,6);
                    for q = 1:length(dQ)
                        if dQ(q)<0; cl='b'; else cl='r'; end
                        qE     = pE;
                        qE     = setfield(qE,P{k},{i,j},dQ(q));
                        [G w]   = spm_csd_mtf(qE,M,[]);
                        GW(:,q) = squeeze(G{1}(:,obs,obs));
                        
                        subplot(4,2,2*iplot - 1), hold on
                        plot(w,abs(GW(:,q)),cl)
                    end
                    
                    % plot
                    %----------------------------------------------------------
                    xlabel('frequency {Hz}')
                    title(sprintf('%s Param: %s(%i,%i)',M.dipfit.model(p).source,P{k},i,j),'FontSize',12)
                    set(gca,'XLim',[0 w(end)]);
                    set(gca,'xtick',[0:10:w(end)])
                    set(gcf,'PaperPositionMode','auto')
                    
                    subplot(4,2,2*iplot - 0)
                    imagesc(dQ,w,log(abs(GW)))
                    title('Transfer functions','FontSize',12)
                    ylabel('Frequency')
                    xlabel('(log) parameter scaling')
                    axis xy; drawnow; colormap jet
                    
                    iplot     = iplot + 1;
                end
            end
        end
    end
    ifig=ifig+1; iplot=1;
    
end

%%%% extrinsic connections

P={'A{1}(1,2)';'A{2}(1,2)';'A{3}(2,1)';'A{4}(2,1)'};
Pname={'tha to ctx/ss';'tha to ctx/dp';'ctx to str';'ctx to stn'};

sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
spm_figure('GetWin',sfig);

for k = 1:length(P)
    
    eval(['Q = pE.',P{k},';']);
    
    if iplot > 4
        iplot = 1;
        ifig  = ifig + 1;
        sfig = sprintf('%s: Parameter dependency - %i',M.dipfit.model(obs).source,ifig);
        spm_figure('GetWin',sfig);
    end
    
    % line search (with solution for steady state)
    %----------------------------------------------------------
    dQ    = linspace(Q - D,Q + D,32);
    for q = 1:length(dQ)
        if dQ(q)<0; cl='b'; else cl='r'; end
        qE     = pE;
        eval(['qE.',P{k},'= dQ(q);']);
        
        [G w]   = spm_csd_mtf(qE,M,[]);
        GW(:,q) = squeeze(G{1}(:,obs,obs));
        
        subplot(4,2,2*iplot - 1), hold on
        plot(w,abs(GW(:,q)),cl)
    end
    
    % plot
    %----------------------------------------------------------
    xlabel('frequency {Hz}')
    title(sprintf('Param: %s',Pname{k}),'FontSize',12)
    set(gca,'XLim',[0 w(end)]);
    
    
    subplot(4,2,2*iplot - 0)
    imagesc(dQ,w,log(abs(GW)))
    title('Transfer functions','FontSize',12)
    ylabel('Frequency')
    xlabel('(log) parameter scaling')
    axis xy; drawnow; colormap jet
    
    % update graphics
    %----------------------------------------------------------
    iplot     = iplot + 1;
    
end

%look at activity single populations
qE=pE;

qE.J{1}=[1 0 0 0 0 0 0 0]; %ss
qE.J{2}=[1 0 0 0 0 0 0 0 0 0]; %str
[G w]   = spm_csd_mtf(qE,M,[]);
G_ss=abs(G{1}(:,1,1));
G_str=abs(G{1}(:,2,2));

qE.J{1}=[0 0 1 0 0 0 0 0]; %sp
qE.J{2}=[0 0 1 0 0 0 0 0 0 0]; %gpe
[G w]   = spm_csd_mtf(qE,M,[]);
G_sp=abs(G{1}(:,1,1));
G_gpe=abs(G{1}(:,2,2));

qE.J{1}=[0 0 0 0 1 0 0 0]; %ii
qE.J{2}=[0 0 0 0 1 0 0 0 0 0]; %stn
[G w]   = spm_csd_mtf(qE,M,[]);
G_ii=abs(G{1}(:,1,1));
G_stn=abs(G{1}(:,2,2));

qE.J{1}=[0 0 0 0 0 0 1 0]; %dp
qE.J{2}=[0 0 0 0 0 0 1 0 0 0]; %gpi
[G w]   = spm_csd_mtf(qE,M,[]);
G_dp=abs(G{1}(:,1,1));
G_gpi=abs(G{1}(:,2,2));

qE.J{1}=[1 0 0 0 0 0 0 0]; %ss
qE.J{2}=[0 0 0 0 0 0 0 0 1 0]; %tha
[G w]   = spm_csd_mtf(qE,M,[]);
G_tha=abs(G{1}(:,2,2));

% figure,subplot(121),plot(w,[G_ss,G_sp,G_ii,G_dp]),legend('ss','sp','ii','dp'),title('CMC')
% subplot(122),plot(w,[G_str,G_gpe,G_stn,G_gpi,G_tha]),legend('str','gpe','stn','gpi','tha'),title('BGC')

figure,subplot(242),plot(w,[G_ss,G_sp,G_ii,G_dp]),legend('ss','sp','ii','dp'),title('MMC')
subplot(245),plot(w,[G_ss]),legend('ss'),title('MMC')
subplot(246),plot(w,[G_sp]),legend('sp'),title('MMC')
subplot(247),plot(w,[G_ii]),legend('ii'),title('MMC')
subplot(248),plot(w,[G_dp]),legend('dp'),title('MMC')

figure,set(gcf,'position',[396    12   460   664])
subplot(242),plot(w,[G_str,G_gpe,G_stn]),legend('str','gpe','stn'),title('BGC')
xlabel('Frequency [Hz]'), set(gca,'ylim',[0 .3],'ytick',[0:.1:1],'xlim',[0 50],'xtick',[0:10:w(end)])
set(gcf,'PaperPositionMode','auto')
subplot(245),plot(w,[G_str]),legend('str'),title('BGC')
subplot(246),plot(w,[G_gpe]),legend('gpe'),title('BGC')
subplot(247),plot(w,[G_stn]),legend('stn'),title('BGC')
%%%%%%

return

qE.J{1}=[0 0 1 0 0 0 1 0]; %ctx
[G w]   = spm_csd_mtf(qE,M,[]);
G_ctx=abs(G{1}(:,1,1));
% figure,subplot(121),plot(w,[G_sp,G_dp,G_ctx]),legend('sp','dp','ctx')


qE=pE;

M.ons=100;
M.dur=2;
U.dt=.001;
t=0:U.dt:.3;
[U.u] = spm_erp_u(t,qE,M);
% U.u(:,2)=U.u(:,1);
U.u(:,2)=zeros(size(U.u(:,1)));


qE.J{1}=[1 0 0 0 0 0 0 0]; %ss
qE.J{2}=[1 0 0 0 0 0 0 0 0 0]; %str
[y,pst] = spm_gen_erp(qE,M,U);
ERP_SS=y{1}(:,1);
ERP_STR=y{1}(:,2);

qE.J{1}=[0 0 1 0 0 0 0 0]; %sp
qE.J{2}=[0 0 1 0 0 0 0 0 0 0]; %gpe
[y,pst] = spm_gen_erp(qE,M,U);
ERP_SP=y{1}(:,1);
ERP_GPE=y{1}(:,2);

qE.J{1}=[0 0 0 0 1 0 0 0]; %ii
qE.J{2}=[0 0 0 0 1 0 0 0 0 0]; %stn
[y,pst] = spm_gen_erp(qE,M,U);
ERP_II=y{1}(:,1);
ERP_STN=y{1}(:,2);

qE.J{1}=[0 0 0 0 0 0 1 0]; %dp
qE.J{2}=[0 0 0 0 0 0 1 0 0 0]; %gpi
[y,pst] = spm_gen_erp(qE,M,U);
ERP_DP=y{1}(:,1);
ERP_GPI=y{1}(:,2);

qE.J{1}=[1 0 0 0 0 0 0 0]; %ss
qE.J{2}=[0 0 0 0 0 0 0 0 1 0]; %tha
[y,pst] = spm_gen_erp(qE,M,U);
ERP_THA=y{1}(:,2);

figure(33),
subplot(211),plot(t,U.u(:,1)./2,'k:',t,ERP_SS,'b',t,ERP_SP,'g',t,ERP_II,'r',t,ERP_DP,'m'); legend('stim','ss','sp','ii','dp'),title('ERP response - cortex')
subplot(212),plot(t,U.u(:,2)./2,'k:',t,ERP_STR,'b',t,ERP_GPE,'g',t,ERP_STN,'r',t,ERP_GPI,'m',t,ERP_THA,'c'); legend('stim','str','gpe','stn','gpi','tha'),title('ERP response - bgc')

figure,plot(t,U.u(:,1)./8,'k:',t,ERP_SS,'b',t,ERP_SP,'g',t,ERP_II,'r',t,ERP_DP,'m'); legend('input u /8','ss','sp','ii','dp'),title('ERP response - CMC model'),xlabel('Time [s]'),ylabel('Membrane potential')




%%%%% look at delay and Jacobian
% figure,subplot(221),imagesc(J),title('J'),colorbar, subplot(222),imagesc(double(Delay)),title('Delay'); colorbar; subplot(223),imagesc(DelayJ),title('DelayJ');colorbar; subplot(224),imagesc(DelayOp),title('DelayOp'),colorbar

% %%%%% look at noise
% [Gu,Gs,Gn,Hz] = spm_csd_mtf_gu(pE,M);
%
%   % plot spectral density of innovations
%     % ---------------------------------------------------------------------
%     figure,
%     subplot(2,1,1)
%     plot(M.Hz,Gu)
%     xlabel('frequency (Hz)')
%     title('Spectrum of innovations','FontSize',16)
%     axis square, grid on, spm_axis scale
%
%
%     % plot spectral density of noise
%     % ---------------------------------------------------------------------
%     subplot(2,2,3)
%     plot(M.Hz,Gs)
%     xlabel('frequency (Hz)')
%     title('Channel-specific noise')
%     axis square, grid on, spm_axis scale
%
%     % plot spectral density of noise
%     % ---------------------------------------------------------------------
%     subplot(2,2,4)
%     plot(M.Hz,Gn)
%     xlabel('frequency (Hz)')
%     title('Non-specific noise')
%     axis square, grid on, spm_axis scale
%



%%%%%%%%%%%%
%%

figure,
plot(M.Hz,exp(0)*M.Hz.^(-exp(0)),'k');
xlim([0 50]),set(gcf,'position',[363   596   131   102])
xlabel('Frequency [Hz]')
set(gca,'box','off')


%%%%%%%%%%%%
%%

delay  = 20;
scale  = 2;
U      = exp(-(M.Hz - delay).^2/(2*scale^2));
figure,
plot(M.Hz,U,'k')
xlim([0 50]),set(gcf,'position',[363   596   131   102])
xlabel('Frequency [Hz]')
set(gca,'box','off')
