% First-level Contrasts - Flexible Models only


this script is disused!! has been integrated with general flexmod contrast procedure


clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={};

% Model details -----------------------
log.firstlevel_contraststype='f1_TaskChoiceTrialtype';

for o1=1:1 % General settings and specifications
    
    % Add paths
    where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];     
    addpath(where.where);  addpath(genpath([where.where filesep '4 Set up models']));
    
    % Which onsets model goes with the requested first-level (contrast) model?
    switch log.firstlevel_contraststype(1)
        case 'f'; log.onsetsmodel='m_f1_ChoicexTrialtype';
        case 't'; log.onsetsmodel='m_t1_Trialtype';
    end
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Settings that don't really change
    log.prefix='swubf'; 
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' ');disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp(['Model (onsets): ' log.onsetsmodel]); 
    disp(['Model (contrasts): ' log.firstlevel_contraststype]); disp(' ')
    disp(' '); input('Hit Enter to start      ')
    disp('======================================================='); disp(' ');disp(' ')
    
end

%% (1) Run contrasts for ALL available regressors (subject-specific)

log.FLest=['2 First level' filesep log.onsetsmodel ' Estimated' filesep];
log.FLcon=['2 First level' filesep log.onsetsmodel(1:4) ' Contrasted   ' log.firstlevel_contraststype filesep];

for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ' (' log.subjects{s}  ')  ---------'])
    ws.FLcon=[where.data_brain filesep log.subjects{s} filesep log.FLcon];
    if isdir(ws.FLcon)==0 ;  copyfile([where.data_brain filesep log.subjects{s} filesep log.FLest],  ws.FLcon); end
    ws.s=load([ws.FLcon 'SPM.mat']);
    ws.regs=ws.s.SPM.xX.name'; c=1; cons={};
    
    % Set up contrasts: Weight every (non-nuisance) regressor
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.FLcon 'SPM.mat']);
    for r=1:size(ws.regs,1)
        if isempty(strfind(ws.regs{r}, '*bf(1)'))==0    &     strcmp(ws.regs{r}(7), 'n')==0  % Criteria for weighting a regressor
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=ws.regs{r}(7:strfind(ws.regs{r}, '*')-1);
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=zeros(1,size(ws.regs,1));
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec(r)=1;
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
            %
            cons{c,1}=matlabbatch{1}.spm.stats.con.consess{c}.tcon.name;
            cons{c,2}=r;
            c=c+1;
        end
    end
    if s==1;  disp('Contrasts for subjects #1: '); disp(cons); input('Check for no nuisance regressors weighted. OK to continue?   '); end
    
    disp(['Running ' num2str(c-1) ' contrasts -----------------'])
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];ws=[];
end