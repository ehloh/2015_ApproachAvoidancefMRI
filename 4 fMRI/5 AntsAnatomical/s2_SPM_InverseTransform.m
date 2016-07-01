% Transform group-level ROIs to native space of individual subjects 
%       (by applying inverse of deformation field used during spatial normalization)
clear all; close all hidden; clc

% Requested analysis
log.specificsubjects={};

% Which model? ########################
log.FirstLevelThread=' s4'; 
for o1=1:1 % Archive models  
% log.onsetsmodel='m_c1_Choice_ENU';
% log.onsetsmodel='m_c2_Choice_ENV'; 
% log.onsetsmodel='m_c3_ChoiceFull_OULPEN'; 
% log.onsetsmodel='m_c4_Choice_OUPEN'; 
% log.onsetsmodel='m_c5_ChoiceRTFull_OULPEN'; 
% log.onsetsmodel='m_c6_ChCluster4Full_OULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; 
% log.onsetsmodel='m_c8_ChCluster4MovFull_OULPEN'; 
% log.onsetsmodel='m_c9_ChCluster6MovFull_OULPEN'; 
% log.onsetsmodel='m_c10_ChCluster6FullRT_OULPEN'; 
% log.onsetsmodel='m_c11_Choice_VOUPEN'; 
% log.onsetsmodel='m_c12_Choice_VOULPEN'; 
% log.onsetsmodel='m_t1_1_Trialtype';
% log.onsetsmodel='m_v1_vChoice_bpm16bpm11';
% log.onsetsmodel='m_v2_vBestAnd_bpm16bpm11';
% log.onsetsmodel='m_v3_vChosenAnd_bpm16bpm11';
% log.onsetsmodel='m_v4_ChoicevChosenAnd_bpm16bpm11';
% log.onsetsmodel='m_v5_ChoiceXvChosenAnd_bpm16bpm11';   
% log.onsetsmodel='m_v6_ChoicevChosenAndposneg2_bpm16bpm11';
% log.onsetsmodel='m_v7_ChClustervChosenAnd_bpm16bpm11';
end
log.onsetsmodel='m_c3_ChoiceFull_OULPEN';     % WHICH FL and SL model?  ############ 
% log.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; 
log.secondlevelmodel='choice_2x3';                 % Second level -----
% log.secondlevelmodel='choice_cluster2x2';
% log.secondlevelmodel='Identity_1samplettest';
% log.secondlevelmodel='Identity_pairedttest';
log.ROIfol='ROI';                             % ROI details   ############ 
log.ROIs={'MEC_001all' ; 'TxC_001all'
%                 'MEC_0001all';  'TxC_0001all'; 
                };

for o1=1:1 % General settings and specifications
        
    % Append Ants Type?
    log.AntsType=[];
    log.onsetsonlymodel=log.onsetsmodel; log.onsetsmodel=[log.onsetsmodel log.AntsType];
    
    
    % Add paths
    where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
%     where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI';  where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];   
    addpath(where.where);  addpath(genpath([where.where filesep '4 Set up models']));
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_c10_';'m_v7_'};
    if sum(strcmp(log.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.onsetsonlymodel);
    end
    
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             First level model:  ' log.onsetsmodel])
    disp(['             Second level model:  ' log.secondlevelmodel]); disp(' ')
    disp('ROI images to transform to native space:'); disp(log.ROIs);
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

log.roiformat='.img';
log.whererois=[where.experiment_folder  filesep  '2 Second level results' log.FirstLevelThread filesep log.onsetsmodel filesep log.secondlevelmodel filesep log.ROIfol filesep];


%%  Batchity apply transformations
for s=1:log.n_subjs
    disp([log.subjects{s} '  ---------------'])
    
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.matname =   cellfun(@(x)[where.data_brain filesep x filesep '1 Preprocessed'  filesep spm_select('List', [where.data_brain filesep x filesep '1 Preprocessed'], '^sM.*.T1w_seg_sn.mat$')], log.subjects(s), 'UniformOutput',0);   % Segmentation files 
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox = [3 3 3];
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {'C:\Users\eloh\Desktop\2 [Explore]\3 Anatomical\AverageT1_n20.nii,1'};  % group structural



matlabbatch{1}.spm.util.defs.ofname = '';
matlabbatch{1}.spm.util.defs.fnames = cellfun(@(x)[ log.whererois x  log.roiformat ',1'], log.ROIs, 'UniformOutput',0);  % ROI file names
% matlabbatch{1}.spm.util.defs.savedir.savesrc = 1;
matlabbatch{1}.spm.util.defs.interp = 1;
matlabbatch{1}.spm.util.defs.interp = 1;

ws.savewhere=['C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\Reverse Transformed Results ROIs\' log.subjects{s}];
matlabbatch{1}.spm.util.defs.savedir.saveusr = {ws.savewhere};
%
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);


% Rename for model and results
ws.namefrom=cellfun(@(x)[ws.savewhere filesep 'w' x  log.roiformat], log.ROIs, 'UniformOutput',0);
ws.nameto=cellfun(@(x)[ws.savewhere filesep log.onsetsmodel(3:strfind(log.onsetsmodel(3:10), '_')+1) '_'  x  log.roiformat], log.ROIs, 'UniformOutput',0);
for i=1:size(ws.namefrom,1)
    ws.old=ws.namefrom{i}; ws.new=ws.nameto{i};
    eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
end
if strcmp(  log.roiformat, '.img') ==1
    for i=1:size(ws.namefrom,1)
        ws.old=[ws.namefrom{i}(1:end-4) '.hdr']; ws.new=[ws.nameto{i}(1:end-4) '.hdr'];
        eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
    end
end


matlabbatch=[]; ws=[];





end