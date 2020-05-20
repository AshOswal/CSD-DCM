% Figures and analysis for paper
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
 
% load empirical DCM and set connectivity for simulations
%==========================================================================
clear all
 
load DCM_CSD
 
% invert model using real (E)mpirical data
%--------------------------------------------------------------------------
DCM = spm_dcm_csd(DCM);
SIM = DCM;

% spm_dcm_csd_results(DCM);
 
% set extrinsic connections (suppress backward connections) - SIM
%--------------------------------------------------------------------------
P           = DCM.Ep;
pE          = DCM.M.pE;
pC          = DCM.M.pC;
 
P.A{1}      = pE.A{1};
P.A{2}      = pE.A{2};
P.A{3}      = pE.A{3};
 
P.A{1}(3:4,1:2) =  log(4);
P.A{2}(1:2,3:4) = -log(2);
 
% P.a         = pE.a + 1;
% P.b         = pE.b - 1;
% P.c         = pE.c - 1;
 

% simulate data
%--------------------------------------------------------------------------
SIM.xY.y = spm_csd_mtf(P,SIM.M,SIM.xU);
 
% adjust starting estimate
%--------------------------------------------------------------------------
SIM.M.P  = spm_dcm_csd_priors(SIM.M,SIM.xU,SIM.xY);
 
% EM: inversion
%--------------------------------------------------------------------------
[Qp,Cp]  = spm_nlsi_GN(SIM.M,SIM.xU,SIM.xY);
 
 
% create DCM structure for display 
%==========================================================================
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc     = spm_csd_mtf(Qp,SIM.M,SIM.xU);                    % prediction
Rc     = spm_unvec(spm_vec(SIM.xY.y) - spm_vec(Hc),Hc);   % prediction error
 

% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
qp     = Qp;
qp.L   = qp.L - qp.L;
qp.b   = qp.b - 32;
qp.c   = qp.c - 32;
Hs     = spm_csd_mtf(qp,DCM.M,DCM.xU);

SIM.Hc = Hc;
SIM.Hs = Hs;
SIM.Rc = Rc;

SIM.Ep = Qp;
SIM.Cp = Cp;
 
% spm_dcm_csd_results(SIM);
 

% repeat but with abs(y) - SSR
%==========================================================================
SSR       = SIM;
SSR.xY.y  = spm_lfp_mtf(P,SSR.M,SSR.xU);
SSR.M.IS  = 'spm_lfp_mtf';
 
% EM: inversion
%--------------------------------------------------------------------------
[sQp,sCp] = spm_nlsi_GN(SSR.M,SSR.xU,SSR.xY);
 
SSR.Ep = sQp;
SSR.Cp = sCp;

% compare priors and posteriors
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'),clf
 
vpC   = spm_vec(pC);
vCp   = spm_vec(diag(Cp));
vQp   = spm_vec(Qp);
vPp   = spm_vec(P);
vpE   = spm_vec(pE);
vsQp  = spm_vec(sQp);
vsCp  = spm_vec(diag(sCp));
 
i     = spm_fieldindices(Qp,'A','D');
i     = i(find(vpC(i) > 1/32));
 
subplot(3,1,1); hold off
spm_plot_ci(vpE(i),vpC(i)), hold on
plot(vPp(i),'*'), hold off
title('True and prior estimates','FontSize',16)
axis square
set(gca,'YLim',[-2 2]);
 
subplot(3,1,2); hold off
spm_plot_ci(vQp(i),vCp(i)), hold on
plot(vPp(i),'*'), hold off
title('Posterior estimates (complex)','FontSize',16)
axis square
set(gca,'YLim',[-2 2]);

subplot(3,1,3); hold off
spm_plot_ci(vsQp(i),vsCp(i)), hold on
plot(vPp(i),'*'), hold off
title('Posterior estimates (modulus)','FontSize',16)
axis square
set(gca,'YLim',[-2 2]);


% Compare conditional estimates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
i     = spm_fieldindices(Qp,'A','S','G','J','C','T','H','a','b','c');
j     = spm_fieldindices(Qp,'D');
 
subplot(2,1,1); hold off
plot(vCp(i),vsCp(i),'.', 'MarkerSize',16), hold on
plot(vCp(j),vsCp(j),'.r','MarkerSize',16), hold on
plot([0 1/8],[0 1/8],':'), hold on
title('Posterior uncertainty','FontSize',16)
xlabel('under complex')
ylabel('under modulus')
axis image
 
 
% Focus on delays (with and without complex part)
%==========================================================================
i1    = 1;  % to
i2    = 3;  % from
 
D   = exp(P.D(i1,i2))*16;
Vp  = spm_unvec(vCp,Qp);
sVp = spm_unvec(vsCp,Qp);
 
x   = linspace(1/128,32,128);
qp  = spm_Npdf(log(x/16), Qp.D(i1,i2), Vp.D(i1,i2));
sp  = spm_Npdf(log(x/16),sQp.D(i1,i2),sVp.D(i1,i2));
pp  = spm_Npdf(log(x/16), pE.D(i1,i2), pC.D(i1,i2));
 
subplot(2,1,2),hold off
plot(x,qp,x,sp,'-.',x,pp,':',[D D],[0 max(qp)])
title('Conditional Delay','FontSize',16)
xlabel('delay (ms)')
ylabel('conditional density')
axis square
legend({'posterior (complex)','posterior (modulus)','prior','true'})


% Conditional density over coherence etc.
%==========================================================================
spm_figure('GetWin','Figure 3a'); clf

i1    = 4;  % to
i2    = 1;  % from

spm_dcm_csd_plot(SIM,i1,i2,1);

% true functions
%----------------------------------------------------------------------
hc        = spm_csd_mtf(P,DCM.M,DCM.xU);
[ccf pst] = spm_csd2ccf(hc,DCM.M.Hz);
[coh fsd] = spm_csd2coh(hc,DCM.M.Hz);
    
    
% Overlay true values
%--------------------------------------------------------------------------
HZ  = DCM.Hz(:);
PST = DCM.pst(:)*1000;
D   = exp(P.D(i1,i2))*16;
 
subplot(2,2,1); hold on
i   = abs(PST) < 128;
plot(PST(i),ccf{1}(i,i1,i2),'LineWidth',4)
 
subplot(2,2,2); hold on
plot(HZ,coh{1}(:,i1,i2),'LineWidth',4)
 
subplot(2,2,3); hold on
plot(HZ,1000*fsd{1}(:,i1,i2),'LineWidth',4)

subplot(2,2,4),hold on
plot([D D],[0 1],'LineWidth',4)


% Now repeat but for sources
%==========================================================================
spm_figure('GetWin','Figure 3b'); clf

i1    = 4;  % to
i2    = 1;  % from

spm_dcm_csd_plot(SIM,i1,i2);

% true functions
%----------------------------------------------------------------------
qp        = P;
qp.L      = qp.L - qp.L;
qp.b      = qp.b - 32;
qp.c      = qp.c - 32;
hc        = spm_csd_mtf(qp,DCM.M,DCM.xU);

[ccf pst] = spm_csd2ccf(hc,DCM.M.Hz);
[coh fsd] = spm_csd2coh(hc,DCM.M.Hz);
    
    
% Overlay true values
%--------------------------------------------------------------------------
HZ  = DCM.Hz(:);
PST = DCM.pst(:)*1000;
D   = exp(P.D(i1,i2))*16;
 
subplot(2,2,1); hold on
i   = abs(PST) < 128;
plot(PST(i),ccf{1}(i,i1,i2),'LineWidth',4)
 
subplot(2,2,2); hold on
plot(HZ,coh{1}(:,i1,i2),'LineWidth',4)
 
subplot(2,2,3); hold on
plot(HZ,1000*fsd{1}(:,i1,i2),'LineWidth',4)

subplot(2,2,4),hold on
plot([D D],[0 1],'LineWidth',4)




% Conditional density over coherence for real data
%==========================================================================
spm_figure('GetWin','Figure 4'); clf

i1    = 4;  % to
i2    = 1;  % from

spm_dcm_csd_plot(DCM,i1,i2);


 
save Figures
