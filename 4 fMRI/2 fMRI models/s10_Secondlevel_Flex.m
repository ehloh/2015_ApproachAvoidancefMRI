% Set up Multiple regression at 2nd level 
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={}; 

% Second level model details (Name of model, FL contrasts) ########################
log.SecondLevelSetup={
%     'ExplorePar'   {'cF_EntropyNTok';'ct_Entropy'}  
%     'vExplore'   {'cF_vExplore';'ct_vExplore'}  
%     'BestUnchosenRev_negMpos'   { 'vBestUnchosen_posrev'; 'vBestUnchosen_negrev'};
%     'cF_BestUnchosenRev_negMpos'   { 'cF_vBestUnchosen_posrev'; 'cF_vBestUnchosen_negrev'};
    'Task_x_RejOrVal'   { 'cF_vBUrev_Rej';'cF_vBUrev_Or'; 'ct_vBUrev_Rej';'ct_vBUrev_Or'};
    
    };

% Which FL model?  ################################################
log.FLthread=' s4Ants';   log.AntsType='_Basic';
% log.FLthread=' s6'; log.AntsType=[];  ------------------------------
% log.FLmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11';
% log.FLmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11';
log.FLmodel='m_v5c_ChoiceXvChosenAnd_bpm16bpmi11';
% log.FLmodel='m_c4_Choice_OUPEN';
% log.FLmodel='m_v1c_vChoice_bpm16bpmi11';
for o1=1:1 % Un-used onsets models 
    % log.FLmodel='m_c2_CompeteBasicVExplore';
    % log.FLmodel='m_c3_CompeteFull_XUVPEN';
    % log.FLmodel='m_c4_CompeteFull_XUPEN';
    % log.FLmodel='m_c5_CompeteFullRT_XUVPEN';
    % log.FLmodel='m_c6_Cluster4CompeteFull_XUVPEN';
    % log.FLmodel='m_c7_Cluster6CompeteFull_XUVPEN';
    % log.FLmodel='m_c8_Cluster4MovCompeteFull_XUVPEN';
    % log.FLmodel='m_c9_Cluster6MovCompeteFull_XUVPEN';
end

for o1=1:1 %  General settings and specifications  
        
    % Things that don't change much
    % where.where='/Users/EleanorL/Dropbox/sandisk/5 Explore fMRI'; where.data_brain=[where.exp_folder filesep '1 MRI data'];
    where.where='D:\Dropbox\SANDISK\5 Explore fMRI';  where.exp_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.exp_folder filesep '1 Brain data'];
    addpath(where.where);  addpath([where.where filesep '4 Set up models' filesep '4 Regression models']); 

    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Change first-level model addresses
    if sum(strcmp(log.FLmodel(1),{'f';'t'}))==1;  log.FLfolname=['m_' log.FLmodel(1:2) ' Contrasted   ' log.FLmodel]; % Flexible or trial-type models
    else log.FLfolname=[log.FLmodel ' Contrasted']; % Parametized models       
    end
    where.subFLfols=cellfun(@(x)[where.data_brain filesep x filesep '2 First level'  log.FLthread filesep log.FLmodel log.AntsType ' Contrasted' filesep], log.subjects, 'UniformOutput', 0);    
    log.allFLcons=load([where.subFLfols{1} 'SPM.mat']); log.allFLcons=cellstr(char(log.allFLcons.SPM.xCon(:).name));
    
    
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_Cluster4CompeteFull_XUVPEN';'m_c7_Cluster6CompeteFull_XUVPEN';'m_c8_Cluster4MovCompeteFull_XUVPEN';'m_c9_Cluster6MovCompeteFull_XUVPEN';'m_c10_Cluster6CompeteRT_XUVPEN'};
    if sum(strcmp(log.FLmodel, w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.FLmodel);
    end
    
    % Second level model folder
    where.secondlevelfol=[where.exp_folder filesep '2 Second level results' log.FLthread filesep log.FLmodel log.AntsType filesep log.SecondLevelSetup{1}  filesep];
    mkdir(where.secondlevelfol)
    
    % Log
%     diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.FLmodel ' - ' log.RegressionModel ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             First level model:              ' log.FLmodel]);
    disp(['             Second level model:              ' log.SecondLevelSetup{1}]);
    disp('               FL contrasts:');  disp(log.SecondLevelSetup{2})
    disp('                 (See log.allFLcons for all listed contrasts)')    
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%%

% Identify FL contrasts (log.SecondLevelSetup{3}{c})
for c=1:length(log.SecondLevelSetup{2})
    if sum(strcmp(log.allFLcons, log.SecondLevelSetup{2}{c}))~=1; disp('Available FL contrast names:');  disp(log.allFLcons); error(['Could not find requested FL contrast:  ' log.SecondLevelSetup{2}{c}]); end
    log.SecondLevelSetup{3}{c}=find(strcmp(log.allFLcons, log.SecondLevelSetup{2}{c}));    
    if log.SecondLevelSetup{3}{c}>9; log.SecondLevelSetup{3}{c}=['con_00' num2str(log.SecondLevelSetup{3}{c})];
    else log.SecondLevelSetup{3}{c}=['con_000' num2str(log.SecondLevelSetup{3}{c})];
    end
end
save([where.secondlevelfol 'details_2ndlevelsetup.mat'], 'log')

%% Manually-specified batches

for o1=1:1 % [1] Paired t-test
% input('Continue with paired t test? ');
% matlabbatch{1}.spm.stats.factorial_design.dir = {where.secondlevelfol};
% matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
% matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
% matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
% matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
% matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
% matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
% matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
% matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
% for s=1:log.n_subjs
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans ={[where.subFLfols{s} log.SecondLevelSetup{3}{1} '.img,1'];   [where.subFLfols{s} log.SecondLevelSetup{3}{2} '.img,1'];};
% end
% spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];
end


% 2x3 ANOVA: Task x Choice (Accept/Reject/Explore, NoBomb/Bomb/Explore)
matlabbatch{1}.spm.stats.factorial_design.dir = {where.secondlevelfol};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Task';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'RejOrVal';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


%%

% Specify Contrast files + Factorial cell 
contrastfiles=cell(size(choices,1),1);
for i=1:size(choices,1)
    % Assign level
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = choices{i,2};
    
    % Assign contrast files
    contrastfiles{i}=cell(length(subjectlist),1);
    for s=1:length(subjectlist)
        contrastfiles{i}{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep choices{i,3} ',1'];
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans =contrastfiles{i};
end

% Include subjects as covariates
for s=1:length(subjectlist)
    sub=zeros(1,length(subjectlist)); sub(s)=1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(choices,1)])';
   matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
   matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
end


%%

% Estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[where.secondlevelfol filesep 'SPM.mat']}; 
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch); matlabbatch=[];


%% END
cd(where.secondlevelfol)
disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]);disp(' ')
disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end; disp(' ');
disp(['Data location (brain): ' where.data_brain]); disp(' ');
disp('Analysis complete for 2nd-level regression model:'); disp(' ')
disp(['             First level model:              ' log.FLmodel]);
disp(['             Second level model:              ' log.SecondLevelSetup{1}]);
disp('               FL contrasts:');  disp(log.SecondLevelSetup{2})
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (' mfilename ' -  '  log.FLmodel ' - ' log.RegressionModel ')']), ' ',1);
catch
end


