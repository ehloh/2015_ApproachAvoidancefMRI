% First-level Contrasts - Flexible & Trialtype models only
%       Each FL model has only its own specified contrasts. Branches from
%       the same Estimated model, and existing FL contrasts are deleted.
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={}    ;
request.execute_firstlevelcontrasts=1;
request.execute_secondlevelmodels=0;

% Model details -----------------------
% log.firstlevel_contraststype='f1_1_ChoicexTrialtype';                   % ##### Contrasts type/model ######## 
% log.firstlevel_contraststype='f1_2_Trialtype';
% log.firstlevel_contraststype='f1_3_ExploreOr_AllQualCells';
% log.firstlevel_contraststype='f1_4_ExploreOr_FixWind4';
% log.firstlevel_contraststype='f1_5_ExploreOr_MovWind';
% log.firstlevel_contraststype='f2_1_ChunkChoicexTrialtype';     
% log.firstlevel_contraststype='f3_1_Event';     
% log.firstlevel_contraststype='t1_2_RLvars';
% log.firstlevel_contraststype='t2_2_RLvars';
% log.firstlevel_contraststype='t3_2_RLvars';
% log.firstlevel_contraststype='t4_2_RLvars';

% MODELS TO RUN
% log.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% log.onsetsmodel='m_c13g_ChoiceFull_ULPEN';
log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
log.onsetsmodel='m_c20g_ChoicePredChoice_bpji08bpji11';

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==0;
        where.where='D:\Dropbox\SANDISK\5 Explore fMRI';  where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
        where.param_scripts='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
    else;  where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI'; where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
        where.param_scripts='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
    end
    addpath(where.where); addpath(where.param_scripts)
    addpath(genpath([where.where filesep '4 Set up models']));
    
    % Which onsets model?
    switch log.firstlevel_contraststype(1:2)
        case 'f1';  log.onsetsmodel='m_f1_ChoicexTrialtype';
        case 'f2';  log.onsetsmodel='m_f2_ChunkChoicexTrialtype';
        case 'f3';  log.onsetsmodel='m_f3_Event';
        case 't1';  log.onsetsmodel='m_t1_Trialtype';
        case 't2';  log.onsetsmodel='m_t2_TrialtypeNc';
        case 't3';  log.onsetsmodel='m_t3_ChunkTrialtype';
        case 't4';  log.onsetsmodel='m_t4_ChunkTrialtypeNc';
        otherwise; error('Onsets models for this FL contrast model has not yet been specified!');
    end
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Settings that don't really change
    log.prefix='swubf'; 
    log.execute_contrasts=request.execute_firstlevelcontrasts;
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.onsetsmodel ' - ' log.firstlevel_contraststype ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' ');disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp(['Model (onsets): ' log.onsetsmodel]); 
    disp(['Model (contrasts): ' log.firstlevel_contraststype]); disp(' ')
    disp(' '); disp('[ CHECK HERE ] WHICH PROCESSES TO EXECUTE?:'); disp(request);
    disp(' '); input('Hit Enter to start      ')
    disp('======================================================='); disp(' ');disp(' ')
    
end

%% (1) Run First-level contrasts function

% 2nd-level Results thread for this 1st-level model
where.firstlevel_resultsfolder=[where.experiment_folder filesep '2 Second level results' filesep log.firstlevel_contraststype];
if isdir(where.firstlevel_resultsfolder)==0; mkdir(where.firstlevel_resultsfolder); end

% Signpost: Execute contrasts or not?
if request.execute_firstlevelcontrasts==1; input('STEP 1: Set up AND execute contrasts   ' )
else input('STEP 1: Fetch first-level contrasts, to prep 2nd level. Contrasts will NOT be executed  ' ); 
end

% error('Pause')

% Run the specified contrast type (execute contrasts or set up 2nd level)
eval(['[Con log] =' log.firstlevel_contraststype '(where,log);']); disp(Con) % Output 'con' is not anything systematic

%% (2) Run all 2nd level models that go with this first-level model

% input('Continue to Run Second-level models?     ');

% error('Pause - ready to run 2nd level')

if request.execute_secondlevelmodels==1
    disp('Running all Second level models requested ####################')
    
    % Run all possible models
    Batch=cell(size(Con.SecondLevelAnalysis,1),1);
    for m=1: size(Con.SecondLevelAnalysis,1)
%         m=2
        disp(['Running 2nd level model # ' num2str(m) '   (' Con.SecondLevelAnalysis{m,1} ') -----------------'])
        eval(['Batch{m}= m_' Con.SecondLevelAnalysis{m,2} '(where,log,Con.SecondLevelAnalysis{m,3});']);
        %
        matlabbatch=Batch{m};
        
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
    end
end

% save('tryflexfac1','matlabbatch')

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed:'); disp(' ')
disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
disp(['Onsets model: ' log.onsetsmodel]); disp(' ')
disp(['Contrasts type: ' log.firstlevel_contraststype]); disp(' ')
disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (' mfilename '  -  '  log.firstlevel_contraststype  ')']), ' ',1);
end
