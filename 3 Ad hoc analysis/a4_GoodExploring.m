% Load data + perform adhoc analysis
clear all; close all hidden; clc

% Request
log.specificsubjects={};

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; else where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; end
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

%% (0) Load data (trialstats)
% 	subjdata: col 2=Collated all data   
%  col 3=cF fmri, col 4=ct fmri
%  col 5=cF training, col 6=ct training
%  col 7=cF training (2nd half), col 8=ct training (2nd half)
col.ExptSect=1;
col.Choice1=2;
col.Choice2=3;
col.ExpBombShown=4;
col.OutcomeMag=5;
col.OutcomeShown=6;
col.FollowedExploreInfo=7; % 1/0/nan
    
for o1=1:1 % [Disused] Data specifications (columns)
    
    % From fMRI raw trialstats
%     col.Task=5; % Event classifiers
%     col.Resp1=8; 
%     col.Resp2=10;
%     col.TrialValid=13; 
%     col.OutcomePresented=14;
%     col.OutcomeMagnitude=15;
%     col.RT1=9;
%     col.Block=17;
%     col.Trialnum=7;
%     
%     col.EnvThreat=3; % Design columns
%     col.NTokens=2;
%     col.pLoss=36;
%     col.Entropy=37;
%     col.VExplore=38;
%     col.EV=41;
%     col.OutcomeMean=39;
%     col.OutcomeVariance=40;
%     col.TrialType=1;
end 

% Load data 
subjdata=cell(log.n_subjs,4);
for s=1:log.n_subjs
    subjdata{s,1}=log.subjects{s};    
%     ws.dl=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_1learnenv.mat']);
    ws.tcf=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_2taskconflict.mat']); ws.tcf=ws.tcf.taskconflict.data;
    ws.tct=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_3taskcontrol.mat']); ws.tct=ws.tct.taskcontrol.data;
    ws.fmrid=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_taskfMRI.mat']);
    ws.cf=ws.fmrid.conflict;  ws.ct=ws.fmrid.control;
    
    % fMRI stage 
    subjdata{s,3}=ws.cf(:, [8 10 18 15 14]);
    subjdata{s,3}=[1+0.*subjdata{s,3}(:,1)   subjdata{s,3}];
    subjdata{s,4}=ws.ct(:, [8 10 18 15 14]);
    subjdata{s,4}=[2+0.*subjdata{s,4}(:,1)   subjdata{s,4}];
    
    % Training stages
    subjdata{s,5}=ws.tcf(:, [8 10 14 15]);
    subjdata{s,5}=[3+0.*subjdata{s,5}(:,1)   subjdata{s,5}   1+0.*subjdata{s,5}(:,1) ];
    subjdata{s,6}=ws.tct(:, [8 10 14 15]);
    subjdata{s,6}=[4+0.*subjdata{s,6}(:,1)   subjdata{s,6}   1+0.*subjdata{s,6}(:,1) ];
    
    % Mark: Did 2nd choice follow exploration?
    for c=3:6
        wc=subjdata{s,c};
        wc(:,col.FollowedExploreInfo)=nan;
        wc(wc(:, col.Choice1)==3,col.FollowedExploreInfo)=0;
%         wc(  wc(:, col.Choice1)==3 & wc(:, col.ExpBombShown)==1 &  wc(:, col.Choice2)==2  ,  col.FollowedExploreInfo)=1;
        wc(  wc(:, col.Choice1)==3 & wc(:, col.ExpBombShown)==0 &  wc(:, col.Choice2)==1  ,  col.FollowedExploreInfo)=1;
        subjdata{s,c}=wc;
        wc=[];
    end
    
    % Training stages 2nd half
    subjdata{s,7}=subjdata{s,5}(floor(size(subjdata{s,5},1)/2):end, :);
    subjdata{s,7}(:,1)=5;
    subjdata{s,8}=subjdata{s,6}(floor(size(subjdata{s,6},1)/2):end, :);
    subjdata{s,8}(:,1)=6;
    
    
    % Combine (i.e. 2nd half of training counted twice!!)
    subjdata{s,2}=[subjdata{s,3}; subjdata{s,4}; subjdata{s,5}; subjdata{s,6}; subjdata{s,7}; subjdata{s,8};];
    
    % Checks?
    ws.expbomb=subjdata{s,2}(subjdata{s,2}(:, 4)==1,:);
    if sum(ws.expbomb(:, col.Choice1)~=3)~=0; error('Explored bombs when subject didnt explore!'); end
    if sum(ws.expbomb(:, col.OutcomeShown)~=1)~=0; error('Explored bombs when exploration not shown!'); end
    
    ws=[];
end

%% (1) Does 2nd step choice follow exploration information? 

doexplore=1;
if doexplore
    d_scores=[[{'Subjects'}; log.subjects] cell(log.n_subjs+1,10)];
    for s=1:log.n_subjs
        ws.cf= subjdata{s,3};
        ws.ct= subjdata{s,4};
        ws.tcf= subjdata{s,5};
        ws.tct= subjdata{s,6};
        ws.t2cf= subjdata{s,7};
        ws.t2ct= subjdata{s,8};
        
        % Overall money earned
        k=2; d_scores{1,k}='o.all';   d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.tcf(:, col.OutcomeMag)) + sum(ws.ct(:, col.OutcomeMag)) + sum(ws.tct(:, col.OutcomeMag)) ;
        k=k+1; d_scores{1,k}='o.allfMRI';    d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.ct(:, col.OutcomeMag)) ;
        k=k+1; d_scores{1,k}='o.cF_all'; d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.tcf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.ct_all';  d_scores{s+1,k}= sum(ws.ct(:, col.OutcomeMag)) + sum(ws.tct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.cF';   d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.ct';   d_scores{s+1,k}= sum(ws.ct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.tcF';   d_scores{s+1,k}= sum(ws.tcf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.tct';   d_scores{s+1,k}= sum(ws.tct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.t2cF';   d_scores{s+1,k}= sum(ws.t2cf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.t2ct';   d_scores{s+1,k}= sum(ws.t2ct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='o.cF_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'o.cF'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'o.tcF'))};
        k=k+1; d_scores{1,k}='o.ct_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'o.ct'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'o.tct'))};
        
        % Only explore trials
        ws.cf= ws.cf(ws.cf(:, col.Choice1)==3, :);
        ws.ct= ws.ct(ws.ct(:, col.Choice1)==3, :);
        ws.tcf= ws.tcf(ws.tcf(:, col.Choice1)==3, :);
        ws.tct= ws.tct(ws.tct(:, col.Choice1)==3, :);
        ws.t2cf= ws.t2cf(ws.t2cf(:, col.Choice1)==3, :);
        ws.t2ct= ws.t2ct(ws.t2ct(:, col.Choice1)==3, :);
%         error
%         
%         %         wc(  wc(:, col.Choice1)==3 & wc(:, col.ExpBombShown)==1 &  wc(:, col.Choice2)==2  ,  col.FollowedExploreInfo)=1;
%         wc(  wc(:, col.Choice1)==3 & wc(:, col.ExpBombShown)==0 &  wc(:, col.Choice2)==1  ,  col.FollowedExploreInfo)=1;
% 
%         col.ExpBombShown
        
        ws.cf= ws.cf(ws.cf(:, col.ExpBombShown)==0, :);   % Only explore trials in which null information shown?
        ws.ct= ws.ct(ws.ct(:, col.ExpBombShown)==0, :);
        ws.tcf= ws.tcf(ws.tcf(:, col.ExpBombShown)==0, :);
        ws.tct= ws.tct(ws.tct(:, col.ExpBombShown)==0, :);
        ws.t2cf= ws.t2cf(ws.t2cf(:, col.ExpBombShown)==0, :);
        ws.t2ct= ws.t2ct(ws.t2ct(:, col.ExpBombShown)==0, :);
        
        
        
        
        % Money earned on explore trials only
        k=k+1; d_scores{1,k}='oe.all';  d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.ct(:, col.OutcomeMag)) + sum(ws.tcf(:, col.OutcomeMag)) + sum(ws.tct(:, col.OutcomeMag)) ;
        k=k+1; d_scores{1,k}='oe.allfMRI';  d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.ct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.cF_all'; d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag)) + sum(ws.tcf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.ct_all';  d_scores{s+1,k}= sum(ws.ct(:, col.OutcomeMag)) + sum(ws.tct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.cF';   d_scores{s+1,k}= sum(ws.cf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.ct';   d_scores{s+1,k}= sum(ws.ct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.tcF';   d_scores{s+1,k}= sum(ws.tcf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.tct';   d_scores{s+1,k}= sum(ws.tct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.t2cF';   d_scores{s+1,k}= sum(ws.t2cf(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.t2ct';   d_scores{s+1,k}= sum(ws.t2ct(:, col.OutcomeMag));
        k=k+1; d_scores{1,k}='oe.cF_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'oe.cF'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'oe.tcF'))};
        k=k+1; d_scores{1,k}='oe.ct_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'oe.ct'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'oe.tct'))};
        
        % Following explore info
        k=k+1; d_scores{1,k}='info.all';   d_scores{s+1,k}= nanmean([ws.cf(:, col.FollowedExploreInfo);  ws.ct(:, col.FollowedExploreInfo); ws.tcf(:, col.FollowedExploreInfo);  ws.tct(:, col.FollowedExploreInfo)]);
        k=k+1; d_scores{1,k}='info.allfMRI';   d_scores{s+1,k}= nanmean([ws.cf(:, col.FollowedExploreInfo);  ws.ct(:, col.FollowedExploreInfo)]);
        k=k+1; d_scores{1,k}='info.cF_all';   d_scores{s+1,k}= nanmean([ws.cf(:, col.FollowedExploreInfo);  ws.tcf(:, col.FollowedExploreInfo)]);
        k=k+1; d_scores{1,k}='info.ct_all';   d_scores{s+1,k}= nanmean([ws.ct(:, col.FollowedExploreInfo);  ws.tct(:, col.FollowedExploreInfo)]);
        k=k+1; d_scores{1,k}='info.cF';  d_scores{s+1,k}=  nanmean(ws.cf(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.ct';  d_scores{s+1,k}=  nanmean(ws.ct(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.tcF';  d_scores{s+1,k}=  nanmean(ws.tcf(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.tct';  d_scores{s+1,k}=  nanmean(ws.tct(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.t2cF';  d_scores{s+1,k}=  nanmean(ws.t2cf(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.t2ct';  d_scores{s+1,k}=  nanmean(ws.t2ct(:, col.FollowedExploreInfo));
        k=k+1; d_scores{1,k}='info.cF_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'info.cF'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'info.tcF'))};
        k=k+1; d_scores{1,k}='info.ct_train2fMRI';  d_scores{s+1,k} = d_scores{s+1, find(strcmp(d_scores(1,:), 'info.ct'))}- d_scores{s+1, find(strcmp(d_scores(1,:), 'info.tct'))};
        
        %
        ws=[];
    end
    
    % Plot post-explore beh 
    a= cell2mat(d_scores(2:end, [find(strcmp(d_scores(1,:), 'info.cF')) find(strcmp(d_scores(1,:), 'info.ct'))])); 
    barwitherr(std(a)/sqrt(log.n_subjs), mean(a))
    
end
openvar d_scores
error('Done! :)')



%% (2) Is value only negative when participants Reject?
% In a RejectOr x vBUvalence (2x2), how many items in each cell?

% Load value 
d_mod=load(['C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\2 Behavioural model details\bpm16bpm11\Model values (07-Sep-2014).mat']);
mcol= d_mod.scol;

% 2x2x2, task x choice (rejector) x valence (+/-)
d_scores=[[{'Subjects'}; log.subjects]  [ {'cF_RejPos' 'cF_RejNeg' 'cF_NonrejPos' 'cF_NonrejNeg' 'ct_RejPos' 'ct_RejNeg' 'ct_NonrejPos' 'ct_NonrejNeg' }; cell(log.n_subjs,8)]]; 
for s=1:log.n_subjs
    for t=1:2        
        wt.d=d_mod.d_trialstats{s,t+1};
        d_scores{s+1, (t-1)*4+2}= sum(wt.d(:, mcol.Choice)==2 & wt.d(:, mcol.vBestUnchosen)>=0); % Reject Pos
        d_scores{s+1, (t-1)*4+3}= sum(wt.d(:, mcol.Choice)==2 & wt.d(:, mcol.vBestUnchosen)<0); % Reject Neg
        d_scores{s+1, (t-1)*4+4}= sum(wt.d(:, mcol.Choice)~=2 & wt.d(:, mcol.vBestUnchosen)>=0); % NonReject Pos
        d_scores{s+1, (t-1)*4+5}= sum(wt.d(:, mcol.Choice)~=2 & wt.d(:, mcol.vBestUnchosen)<0); % NonReject Neg
    end
end
openvar d_scores



mean(cell2mat(d_scores(2:end, 2:end)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    