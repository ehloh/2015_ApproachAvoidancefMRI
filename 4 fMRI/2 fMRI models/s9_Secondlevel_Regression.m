% Set up Multiple regression at 2nd level 
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={}; % {'p04_MP'; 'p02_YY'};

% Which regression model? (Regression, Firstlevel, PPI) ########################
log.BehVariables={   % Name of model, behavioural parameters included, Mean-centreing or not (1/0)
%     'TaskExplore'   {'per.cF_Explore'; 'per.ct_Explore'}   [1 1];                    % Behavioural variable ---------------
%     'TraitAnxiety'        {'st.Trait'}   1;          
    'StateAnxiety'        {'st.State'}   1;      
%     'perReject'   {'per.cF_Reject'}   [1];      
    };

% Which FL model?  ################################################
log.FLthread=' s4Ants';   log.AntsType='_Basic';
% log.FLthread=' s6'; log.AntsType=[];  ------------------------------
% log.FLmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11';  log.FLcontrast='cF_vBestUnchosen';
log.FLmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11';
% log.FLcontrast='vBestUnchosen_negrev';
% log.FLcontrast='vBestUnchosen_posrev';
log.FLcontrast='cF_vBestUnchosen_posrev';


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
    log.behprofile='Behavioural profile (23-Aug-2013)';
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Change first-level model addresses
    if sum(strcmp(log.FLmodel(1),{'f';'t'}))==1;  log.FLfolname=['m_' log.FLmodel(1:2) ' Contrasted   ' log.FLmodel]; % Flexible or trial-type models
    else log.FLfolname=[log.FLmodel ' Contrasted']; % Parametized models       
    end
    where.subFLfols=cellfun(@(x)[where.data_brain filesep x filesep '2 First level'  log.FLthread filesep log.FLmodel log.AntsType ' Contrasted' filesep], log.subjects, 'UniformOutput', 0);
    log.FLcons=load([where.subFLfols{1} 'SPM.mat']); log.FLcons=cellstr(char(log.FLcons.SPM.xCon(:).name));
    if sum(strcmp(log.FLcons, log.FLcontrast))~=1;  disp('Available FL contrast names:');  disp(log.FLcons);  error(['Could not find requested FL contrast:  ' log.FLcontrast ]); end
    log.FLconnum= find(strcmp(log.FLcons, log.FLcontrast));    
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_Cluster4CompeteFull_XUVPEN';'m_c7_Cluster6CompeteFull_XUVPEN';'m_c8_Cluster4MovCompeteFull_XUVPEN';'m_c9_Cluster6MovCompeteFull_XUVPEN';'m_c10_Cluster6CompeteRT_XUVPEN'};
    if sum(strcmp(log.FLmodel, w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.FLmodel);
    end
    
    % Regression model details
    log.RegModName=[log.FLcontrast '_x_' log.BehVariables{1}];    
    where.secondlevelfol=[where.exp_folder filesep '2 Second level results' log.FLthread filesep log.FLmodel log.AntsType filesep log.RegModName filesep];
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
    disp(['             Regression model name :  ' log.RegModName]);
    disp('             Behavioural variables:         '); disp(log.BehVariables{2}); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% Load behaviour & pick requested variable

% Load behavioural file + subject selection
[n t r]=xlsread([where.where filesep '1 Behavioural data' filesep 'Group behaviour files' filesep 'Behavioural profile for correlation.xlsx'], 'BehInhibAndExplore');
d_beh=[t(:,1) [t(1, 2:end); num2cell(n)]]; 
[log.subjects log.n_subjs d_beh] = f_selectsubjects(d_beh, log.specificsubjects,  [d_beh(:,1) [{'All'}; num2cell(ones(size(d_beh,1)-1,1))]], 'All');  % Subject selection

% Select behavioural variables
log.behvar=nan*zeros(log.n_subjs,length(log.BehVariables{2}));
for b=1:length(log.BehVariables{2})
    if sum(strcmp(d_beh(1,:), log.BehVariables{2}{b}))==0; error(['Invalid behavioural variable requested (' log.BehVariables{2}{b} ')']); end
    log.behvar(:,b)=cell2mat(d_beh(2:end,  find(strcmp(d_beh(1,:), log.BehVariables{2}{b}), 1, 'first')));
end

%% Set up 2nd level regression model

matlabbatch{1}.spm.stats.factorial_design.dir = {where.secondlevelfol};
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans =cellfun(@(x)[x spm_select('List',  where.subFLfols{1}, ['^con_.*.0' num2str(log.FLconnum) '.img']) ',1'], where.subFLfols, 'UniformOutput', 0);  % Scans

% Behavioural covariates
for b=1: length(log.BehVariables{2})
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(b).c=log.behvar(:, b);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(b).cname=log.BehVariables{2}{b};
    if log.BehVariables{3}(b) ==1; matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(b).iCC = 1;
    else matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(b).iCC = 5;  % % 1 = mean centre, 5 = No mean centring
    end
end

matlabbatch{1}.spm.stats.factorial_design.cov.c = 1:log.n_subjs; % Subject covariates
matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
save([where.secondlevelfol 'details_secondlevel.mat'], 'log', 'd_beh')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% Estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[where.secondlevelfol filesep 'SPM.mat']}; 
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

%% (3) Add 2nd level contrasts positive and negative contrasts

disp('Adding 2nd-level contrasts to regression model ######################')
matlabbatch{1}.spm.stats.con.spmmat =  {[where.secondlevelfol filesep 'SPM.mat']}; 
matlabbatch{1}.spm.stats.con.delete = 1; k=1;

% Contrasts for all behavioural variables
for i=1: length(log.BehVariables{2})
    
    % Positive
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.name = [log.BehVariables{2}{i} '_Pos'];
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.convec =[zeros(1,i+1) 1];
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
    k=k+1;
    
    % Negative
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.name = [log.BehVariables{2}{i} '_Neg'];
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.convec =[zeros(1,i+1) -1];
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
    k=k+1;
    
end
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[]; 

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]);disp(' ')
disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end; disp(' ');
disp(['Data location (brain): ' where.data_brain]); disp(' ');
disp('Analysis complete for 2nd-level regression model:'); disp(' ')
disp(['             First level model:              ' log.FLmodel]);
disp(['             Regression model name :  ' log.RegModName]);
disp('             Behavioural variables:         '); disp(log.BehVariables{2}); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (' mfilename ' -  '  log.FLmodel ' - ' log.RegressionModel ')']), ' ',1);
catch
end


