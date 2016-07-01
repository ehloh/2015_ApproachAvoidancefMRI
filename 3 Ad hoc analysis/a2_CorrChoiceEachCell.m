% Load data + perform adhoc analysis
clear all; close all hidden; clc

% Request
log.specificsubjects={};

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI'; else where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; end
    addpath(where.where); where.data=[where.where filesep '1 Behavioural data']; 

    % Load subjects    
    log.w=load([where.data filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp('=======================================================')
    
end

%% (1) Load data (trialstats)
% 	subjdata: col 2=cF trialstats, col 3= ct trialstats, 
%       col 4= all trialstats, col 5=behaviour stats (pre-calculated)

for o1=1:1 % Data specifications (columns)
    
    col.Task=5; % Event classifiers
    col.Resp1=8; 
    col.Resp2=10;
    col.TrialValid=13; 
    col.OutcomePresented=14;
    col.OutcomeMagnitude=15;
    col.RT1=9;
    col.Block=17;
    col.Trialnum=7;
    
    col.EnvThreat=3; % Design columns
    col.NTokens=2;
    col.pLoss=36;
    col.Entropy=37;
    col.VExplore=38;
    col.EV=41;
    col.OutcomeMean=39;
    col.OutcomeVariance=40;
    col.TrialType=1;
end 

% Load data (trialstats + stats) 
subjdata=cell(log.n_subjs,4);
for s=1:log.n_subjs
    ws.dd=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_taskfMRI.mat']);
    subjdata{s,1}=log.subjects{s};
    subjdata{s,2}=ws.dd.conflict;
    subjdata{s,3}=ws.dd.control;
    subjdata{s,4}=sortrows([ws.dd.conflict; ws.dd.control], col.Trialnum);
    subjdata{s,5}.cF=ws.dd.cF_stats;
    subjdata{s,5}.ct=ws.dd.ct_stats;
    ws=[];
end
    
%% Correlate choices across all types of trials

% Set up data
d_phoicemat_cF=repmat({nan(log.n_subjs,3)}, 6,6);  % pChoice for each subject (6x6)
d_phoicemat_ct=repmat({nan(log.n_subjs,3)}, 6,6);
d_pchoice_incluster4=nan(log.n_subjs,6); % pChoice in in-cluster (4 cell): pAcceptcF, pRejectcF,pExplorecF, pNoBomb, pBombct,pExplorect
d_pchoice_incluster6=nan(log.n_subjs,6); % pChoice in in-cluster (6 cell)
log.incluster_min=[]; % 0.05; % Empty to NOT exclude
for s=1:log.n_subjs
    for e=1:6
        for n=1:6
            d_phoicemat_cF{e,n}(s,:)=[subjdata{s,5}.cF{e,n}.p_Accept subjdata{s,5}.cF{e,n}.p_Reject subjdata{s,5}.cF{e,n}.p_Explore];
            d_phoicemat_ct{e,n}(s,:)=[subjdata{s,5}.ct{e,n}.p_NoBomb subjdata{s,5}.ct{e,n}.p_Bomb subjdata{s,5}.ct{e,n}.p_Explore];
        end
    end
    
    % InCluster 4
    ws.in4_Accept_cF= subjdata{s,5}.cF{3,5}.n_Accept+ subjdata{s,5}.cF{3,6}.n_Accept+ subjdata{s,5}.cF{4,5}.n_Accept+subjdata{s,5}.cF{4,6}.n_Accept;
    ws.in4_Reject_cF= subjdata{s,5}.cF{3,5}.n_Reject+ subjdata{s,5}.cF{3,6}.n_Reject+ subjdata{s,5}.cF{4,5}.n_Reject+subjdata{s,5}.cF{4,6}.n_Reject;
    ws.in4_Explore_cF= subjdata{s,5}.cF{3,5}.n_Explore+ subjdata{s,5}.cF{3,6}.n_Explore+ subjdata{s,5}.cF{4,5}.n_Explore+subjdata{s,5}.cF{4,6}.n_Explore;
    ws.in4_ncF=ws.in4_Accept_cF+ws.in4_Reject_cF+ws.in4_Explore_cF;
    ws.in4_NoBomb_ct= subjdata{s,5}.ct{3,5}.n_NoBomb+ subjdata{s,5}.ct{3,6}.n_NoBomb+ subjdata{s,5}.ct{4,5}.n_NoBomb+subjdata{s,5}.ct{4,6}.n_NoBomb;
    ws.in4_Bomb_ct= subjdata{s,5}.ct{3,5}.n_Bomb+ subjdata{s,5}.ct{3,6}.n_Bomb+ subjdata{s,5}.ct{4,5}.n_Bomb+subjdata{s,5}.ct{4,6}.n_Bomb;
    ws.in4_Explore_ct= subjdata{s,5}.ct{3,5}.n_Explore+ subjdata{s,5}.ct{3,6}.n_Explore+ subjdata{s,5}.ct{4,5}.n_Explore+subjdata{s,5}.ct{4,6}.n_Explore;
    ws.in4_nct=ws.in4_NoBomb_ct+ws.in4_Bomb_ct+ws.in4_Explore_ct;
    
    % InCluster 6
    ws.in6_Accept_cF= subjdata{s,5}.cF{3,4}.n_Accept+ subjdata{s,5}.cF{3,5}.n_Accept+ subjdata{s,5}.cF{3,6}.n_Accept+ subjdata{s,5}.cF{4,4}.n_Accept+subjdata{s,5}.cF{4,5}.n_Accept+subjdata{s,5}.cF{4,6}.n_Accept;
    ws.in6_Reject_cF= subjdata{s,5}.cF{3,4}.n_Reject+ subjdata{s,5}.cF{3,5}.n_Reject+ subjdata{s,5}.cF{3,6}.n_Reject+ subjdata{s,5}.cF{4,4}.n_Reject+subjdata{s,5}.cF{4,5}.n_Reject+subjdata{s,5}.cF{4,6}.n_Reject;
    ws.in6_Explore_cF= subjdata{s,5}.cF{3,4}.n_Explore+subjdata{s,5}.cF{3,5}.n_Explore+ subjdata{s,5}.cF{3,6}.n_Explore+ subjdata{s,5}.cF{4,4}.n_Explore+subjdata{s,5}.cF{4,5}.n_Explore+subjdata{s,5}.cF{4,6}.n_Explore;
    ws.in6_ncF=ws.in6_Accept_cF+ws.in6_Reject_cF+ws.in6_Explore_cF;
    ws.in6_NoBomb_ct= subjdata{s,5}.ct{3,4}.n_NoBomb+ subjdata{s,5}.ct{3,5}.n_NoBomb+ subjdata{s,5}.ct{3,6}.n_NoBomb+ subjdata{s,5}.ct{4,4}.n_NoBomb+subjdata{s,5}.ct{4,5}.n_NoBomb+subjdata{s,5}.ct{4,6}.n_NoBomb;
    ws.in6_Bomb_ct= subjdata{s,5}.ct{3,4}.n_Bomb+ subjdata{s,5}.ct{3,5}.n_Bomb+ subjdata{s,5}.ct{3,6}.n_Bomb+ subjdata{s,5}.ct{4,4}.n_Bomb+subjdata{s,5}.ct{4,5}.n_Bomb+subjdata{s,5}.ct{4,6}.n_Bomb;
    ws.in6_Explore_ct= subjdata{s,5}.ct{3,4}.n_Explore+ subjdata{s,5}.ct{3,5}.n_Explore+ subjdata{s,5}.ct{3,6}.n_Explore+ subjdata{s,5}.ct{4,4}.n_Explore+subjdata{s,5}.ct{4,5}.n_Explore+subjdata{s,5}.ct{4,6}.n_Explore;
    ws.in6_nct=ws.in6_NoBomb_ct+ws.in6_Bomb_ct+ws.in6_Explore_ct;
    
    % Calculate stats for in-cluster (if minimum n applies)
    d_pchoice_incluster6(s,:) =[ws.in6_Accept_cF/ws.in6_ncF ws.in6_Reject_cF/ws.in6_ncF ws.in6_Explore_cF/ws.in6_ncF      ws.in6_NoBomb_ct/ws.in6_nct ws.in6_Bomb_ct/ws.in6_nct ws.in6_Explore_ct/ws.in6_nct];
    d_pchoice_incluster4(s,:) =[ws.in4_Accept_cF/ws.in4_ncF ws.in4_Reject_cF/ws.in4_ncF ws.in4_Explore_cF/ws.in4_ncF      ws.in4_NoBomb_ct/ws.in4_nct ws.in4_Bomb_ct/ws.in4_nct ws.in4_Explore_ct/ws.in4_nct];
    if isempty(log.incluster_min)==0
        if sum(d_pchoice_incluster6(s,:)<log.incluster_min)>1; d_pchoice_incluster6(s,:)=nan; end
        if sum(d_pchoice_incluster4(s,:)<log.incluster_min)>1; d_pchoice_incluster4(s,:)=nan; end
    end
    
end 
if isempty(log.incluster_min)==0
    d_pchoice_incluster6=d_pchoice_incluster6(isnan(d_pchoice_incluster6(:,1))==0,:); % remove subjects without
    d_pchoice_incluster4=d_pchoice_incluster4(isnan(d_pchoice_incluster4(:,1))==0,:);
end

% Get correlation stats (print r if sig)
log.corrthresholdp=0.05;
rc_ar_cF=nan(6,6); rc_ae_cF=nan(6,6); rc_re_cF=nan(6,6);
rp_ar_cF=nan(6,6); rp_ae_cF=nan(6,6); rp_re_cF=nan(6,6);
rc_ar_ct=nan(6,6); rc_ae_ct=nan(6,6); rc_re_ct=nan(6,6);
rp_ar_ct=nan(6,6); rp_ae_ct=nan(6,6); rp_re_ct=nan(6,6);
for e=1:6
    for n=1:6
        
        % cF
        [wc.r wc.p]=corrcoef(d_phoicemat_cF{e,n});
        wc.r=wc.r.*(wc.p<log.corrthresholdp); % apply p threshold filtering
        rc_ar_cF(e,n) =wc.r(2,1);
        rc_ae_cF(e,n) =wc.r(3,1);
        rc_re_cF(e,n) =wc.r(2,3);
        rp_ar_cF(e,n) =wc.p(2,1); % pvals
        rp_ae_cF(e,n) =wc.p(3,1);
        rp_re_cF(e,n) =wc.p(2,3);
        
        % ct
        [wc.r wc.p]=corrcoef(d_phoicemat_ct{e,n});
        wc.r=wc.r.*(wc.p<log.corrthresholdp); % apply p threshold filtering
        rc_ar_ct(e,n) =wc.r(2,1);
        rc_ae_ct(e,n) =wc.r(3,1);
        rc_re_ct(e,n) =wc.r(2,3);
        rp_ar_ct(e,n) =wc.p(2,1); % pvals
        rp_ae_ct(e,n) =wc.p(3,1);
        rp_re_ct(e,n) =wc.p(2,3);
        
        
        
        wc=[];
    end
end

% Plot InClusters
figure('Name', ['In cluster correlation r (sig only), min % choice=' num2str(log.incluster_min)], 'Color','w');
subplot(2,2,1);  % In cluster 4
[r p]=corrcoef(d_pchoice_incluster4(:,1:3));  


r=r.*(p<=log.corrthresholdp);
[rr pp]=corrcoef(d_pchoice_incluster4(:,4:6)); rr=rr.*(pp<=log.corrthresholdp);
bar([r(1,2) r(1,3) r(2,3); rr(1,2) rr(1,3) rr(2,3)]); axis([0 3   -1 1]); title(['InCluster4 (n=' num2str(size(d_pchoice_incluster4,1)) ')'])
legend(gca,['A v R';'A v E'; 'R v E'],'Location', 'Best'); set(gca,'XTick', [1 2], 'XTickLabel', {'cF';'ct'})
subplot(2,2,3);
bar(1-[p(1,2) p(1,3) p(2,3); pp(1,2) pp(1,3) pp(2,3)]); axis([0 3   0.9 1]);  title('InCluster4 (1-p)')




subplot(2,2,2);   % In cluster 6
[r p]=corrcoef(d_pchoice_incluster6(:,1:3));  r=r.*(p<=log.corrthresholdp);
[rr pp]=corrcoef(d_pchoice_incluster6(:,4:6)); rr=rr.*(pp<=log.corrthresholdp);
bar([r(1,2) r(1,3) r(2,3); rr(1,2) rr(1,3) rr(2,3)]); axis([0 3   -1 1]);   title(['InCluster6 (n=' num2str(size(d_pchoice_incluster6,1)) ')'])
legend(gca,['A v R';'A v E'; 'R v E'],'Location', 'Best'); set(gca,'XTick', [1 2], 'XTickLabel', {'cF';'ct'})
subplot(2,2,4);
bar(1-[p(1,2) p(1,3) p(2,3); pp(1,2) pp(1,3) pp(2,3)]); axis([0 3   0.9 1]);  title('InCluster6 (1-p)')



% Plot: Correlations across entire design
figure('Name', 'Correlation r (sig only)', 'Color','w');
m=2; n=3; p=1; range=[-1 1];
log.st_gapHor  =0.001; log.st_gapVer =0.05; log.st_margH =0.01; log.st_margV =0.01;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(rc_ar_cF, range); title('[cF] A v R'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV);imagesc(rc_ae_cF, range); title('[cF] A v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(rc_re_cF, range); title('[cF] R v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(rc_ar_ct, range); title('[ct] A v R'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(rc_ae_ct, range); title('[ct] A v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(rc_re_ct, range); title('[ct] R v E'); axis square; axis off; colorbar;p=p+1;
figure('Name', 'Correlation 1-p vals (hotter = > sig)', 'Color','w');
m=2; n=3; p=1; range=[0.95 1];
log.st_gapHor  =0.001; log.st_gapVer =0.05; log.st_margH =0.01; log.st_margV =0.01;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(1-rp_ar_cF, range); title('[cF] A v R'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV);imagesc(1-rp_ae_cF, range); title('[cF] A v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(1-rp_re_cF, range); title('[cF] R v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(1-rp_ar_ct, range); title('[ct] A v R'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(1-rp_ae_ct, range); title('[ct] A v E'); axis square; axis off; colorbar;p=p+1;
subtightplot(m,n, p, [log.st_gapHor  log.st_gapVer], log.st_margH, log.st_margV); imagesc(1-rp_re_ct, range); title('[ct] R v E'); axis square; axis off; colorbar;p=p+1;




%%


