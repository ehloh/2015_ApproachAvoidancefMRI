% Does choice differ by task? 2x6x6 ANOVAs

clear all; close all hidden; clc

% 
where.where='C:\Users\e.loh\Dropbox\SCRIPPS\5 Explore fMRI';
 
log.specificsubjects={};

for o1=1:1 % General settings and specifications
    
    % Add paths
    addpath(where.where);
    where.data=[where.where filesep '1 Behavioural data']; 

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

%% (1) Normal analysis, all data 

for o1=1:1
% % Load data 
% subjdata=cell(log.n_subjs,2);
% for s=1:log.n_subjs
%     subjdata{s,1}=log.subjects{s};
%     subjdata{s,2}=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_taskfMRI.mat']);
% end
% 
% % Write to a txt file - easier to do in SPSS
% d=cell(log.n_subjs+1, 2*6*6); d(2:end,1)=log.subjects;
% d_a=d; d_r=d; d_e=d; k=1;
% for t=1:2
%     switch t
%         case 1
%             wc.t='cF';
%             wc.a='Accept';
%             wc.b='Reject';
%         case 2
%             wc.t='ct';
%             wc.a='NoBomb';
%             wc.b='Bomb';
%     end
%     for e=1:6
%         for n=1:6
%             % Headers
%             d_a{1,1+k}=[wc.t '_e' num2str(e) '_n' num2str(n)];
%             d_r{1,1+k}=[wc.t '_e' num2str(e) '_n' num2str(n)];
%             d_e{1,1+k}=[wc.t '_e' num2str(e) '_n' num2str(n)];
%             
%             % Load subject data
%             for s=1:log.n_subjs
%                 eval(['ws=subjdata{s,2}.' wc.t '_stats;'])
%                 eval(['d_a{s+1,k+1}=ws{7-e,n}.p_' wc.a ';'])
%                 eval(['d_r{s+1,k+1}=ws{7-e,n}.p_' wc.b ';'])
%                 d_e{s+1,k+1}=ws{7-e,n}.p_Explore; 
%                 ws=[];
%             end
%             
%             k=k+1;
%         end
%     end
% end
% print2txt([where.data filesep 'Group behaviour files'], 'Accept_TaskxEnvxNTok', d_a)
% print2txt([where.data filesep 'Group behaviour files'], 'Reject_TaskxEnvxNTok', d_r)
% print2txt([where.data filesep 'Group behaviour files'], 'Explore_TaskxEnvxNTok', d_e)
end

%% (2) Analyse data split within block (looking for cross-task contamination)
% Adapt this script to perform this analysis purely from the trialstats

splitblock=1;
task='conflict';
prefix='cF1_';

for o1=1:1 % Data specifications (columns)
    
    col.Task=5; % Event classifiers
    col.Resp1=8; 
    col.Resp2=10;
    col.TrialValid=13; 
    col.OutcomePresented=14;
    col.OutcomeMagnitude=15;
    col.RT1=9;
    col.Block=17;
    
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

% Load data (trialstats) - alter here to analyse alternative data
subjdata=cell(log.n_subjs,1);
p.accept=cell(6,6); p.reject=cell(6,6); p.explore=cell(6,6);
for s=1:log.n_subjs
    ws.dd=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_splitblock.mat']);
    eval(['ws.d= ws.dd.splitblockfMRI.' task num2str(splitblock) 'data;']);
    subjdata{s}=ws.d;
    
    % Format trialstats into 6x6
    for e=1:6
        for n=1:6
            wg=ws.d(ws.d(:,col.EnvThreat)==e & ws.d(:,col.NTokens)==n, col.Resp1);
            %
            p.accept{7-e, n}(s,1)=sum(wg==1)/length(wg);
            p.reject{7-e, n}(s,1)=sum(wg==2)/length(wg);
            p.explore{7-e, n}(s,1)=sum(wg==3)/length(wg);
            wg=[];
        end
    end
end

% Print to SPSS-format columns
d=cell(log.n_subjs, 6*6+1); d{1,1}='Subject';
for s=1:log.n_subjs 
    d{s+1,1}=log.subjects{s};
end
d_a=d; d_r=d; d_e=d; k=2;
for e=1:6
    for n=1:6
        d_a{1,k}=[prefix 't' num2str(e) '_' num2str(n)'];
        d_r{1,k}=[prefix 't' num2str(e) '_' num2str(n)'];
        d_e{1,k}=[prefix 't' num2str(e) '_' num2str(n)'];
        %
        d_a(2:end, k)= num2cell(p.accept{7-e,n});
        d_r(2:end, k)= num2cell(p.reject{7-e,n});
        d_e(2:end, k)= num2cell(p.explore{7-e,n});
        %
        k=k+1; 
    end
end

doprint=0;
if doprint
    print2txt([where.data filesep 'Group behaviour files' filesep 'SplitWithinBlock'], [prefix 'Accept_TaskxEnvxNTok'], d_a)
    print2txt([where.data filesep 'Group behaviour files' filesep 'SplitWithinBlock'], [prefix 'Reject_TaskxEnvxNTok'], d_r)
    print2txt([where.data filesep 'Group behaviour files' filesep 'SplitWithinBlock'], [prefix 'Explore_TaskxEnvxNTok'], d_e)
end

%%
 