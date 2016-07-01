% Check mask.images
clear all;close all hidden; clc

% where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.data_brain='/Volumes/SANDISK/2 EXPLORE Brain data';
where.where='I:\5 Explore fMRI'; where.data_folder= 'C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.data_folder filesep '1 Brain data']; % where.data_beh=[where.where filesep 'Behavioural data'];

% Requested analysis
log.specificsubjects={}; % BLANK to process all subjects

% Which model? --------------
log.model='m_c1_CompeteBasic';
% log.model='m_o1_OrthogBasic';
% log.model='m_t1_ChoicexTrialtype';
% log.model='m_t2_Trialtype';

for o1=1:1 % General settings and specifications
           
    % Subjects
    allsubjects=1;
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    if allsubjects==0
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep 'i8_subjectdataok_secondlevel.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.model);
    else    
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
        disp('ALL SUBJECTS FOR NOW')
    end

    % Settings that don't really change
    log.prefix='swu';
    
    % Interface
    disp('=======================================================')
    w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp('   Subset of subjects only:')
        disp(log.specificsubjects)
    end
    disp(' ')
    disp(['Data location (brain): ' where.data_brain])
    disp(' ')
    disp(['Model: ' log.model])
    disp(' ')
%     input('Hit Enter to start      ')
    disp('=======================================================')
    
end

spm('defaults', 'FMRI');spm_jobman('initcfg'); 

%% Collect scans

masks=cell(log.n_subjs,1);mask_command=cell(log.n_subjs,1); 
for s=1:log.n_subjs
    wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.model ' Contrasted' filesep];
    masks{s}=[wb.where 'mask.img,1'];
    mask_command{s}=['matlabbatch{1}.spm.util.disp.data = {''' masks{s} '''};'];
end

%% Construct global mask scan

% Assemble command to calculate global mask
com='('; cd(where.data_brain); cd ..; w.here=pwd;
for s=1:log.n_subjs
    com=[com 'i' num2str(s)];
    if s==log.n_subjs
        com=[com ')/' num2str(log.n_subjs)];
    else
        com=[com '+'];
    end
end

% Save
subjects=log.subjects; save([w.here filesep '3 Checks' filesep log.model 'Globalmask_details.mat'], 'subjects', 'masks');

% Create global mask
matlabbatch{1}.spm.util.imcalc.input = masks;
matlabbatch{1}.spm.util.imcalc.output =[log.model  ' Globalmask_n' num2str(log.n_subjs)  '.img'];
matlabbatch{1}.spm.util.imcalc.outdir = {[w.here filesep '3 Checks']};
matlabbatch{1}.spm.util.imcalc.expression =com;
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);
matlabbatch=[];

% matlabbatch{1}.spm.util.imcalc.input = {
%                                         'C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data\p01_CW\2 First level\m_ci3_ContextItemNomotor_Hit Contrasted\mask.img,1'
%                                         'C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\1 MRI data\p02_MK\2 First level\m_ci3_ContextItemNomotor_Hit Contrasted\mask.img,1'
%                                         };
% matlabbatch{1}.spm.util.imcalc.output = 'output.img';
% matlabbatch{1}.spm.util.imcalc.outdir = {'C:\Users\eloh\Desktop\1 [Context-Memory] fMRI Data\'};
% matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2)/2';
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

%% Flexible display (Single subject)

input('Continue to single subject display?')

% On-going display
disp('Starting flexible display #############################')
a=1; 
while a~=0 
    disp('Displaying mask image ------- ')
    disp(['Model = ' log.model])
    disp(['Subject=' log.subjects{a} ' (' num2str(a) ')']);
    eval(mask_command{a}) 
    spm_jobman('run',matlabbatch);
    matlabbatch=[];
    %
    disp('######################################')
    disp('Specify which subject to display next')
    disp('Input Enter to go to next; Input 0 to end displaying')
    a1=input(['Next subject no. (n=' num2str(log.n_subjs) ')    ']);
    if isempty(a1)==1
        a=a+1;
    elseif a==log.n_subjs
        a=1;
    else
        a=a1;
    end
end



