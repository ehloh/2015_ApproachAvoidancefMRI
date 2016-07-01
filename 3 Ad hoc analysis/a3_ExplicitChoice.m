% Load data + perform adhoc analysis
clear all; close all hidden; clc

% Request
log.specificsubjects={};

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI'; else 
        where.where='C:\Users\e.loh\Dropbox\SCRIPPS\5 Explore fMRI'; end
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
    subjdata{s,1}=log.subjects{s};
%     ws.dl=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_1learnenv.mat']);
%     ws.dcf=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_2taskconflict.mat']);
    ws.dct=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_3taskcontrol.mat']);
    disp(' '); disp(log.subjects{s})
    disp([ws.dct.taskcontrol.settings.disp.contingencycolours ws.dct.taskcontrol.settings.disp.contingencycolours_original(:,4)])
    subjdata{s,2}=input('Recheck?   (0=Learning ok, 1=Learning wrong)    ');
    ws=[];
end
    


% CHOICE TASK ================
%     Col 1: Trial number
%     Col 2: Left Colour #
%     Col 3: Right Colour #
%     Col 4: Better colour (correct choice - 1=Left, 2=Right)
%     Col 5: Choice (1=Left, 2=Right)
%     Col 6: Accuracy