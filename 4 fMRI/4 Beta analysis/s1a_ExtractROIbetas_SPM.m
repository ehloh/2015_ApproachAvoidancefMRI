% Extract ROI betas using SPM functions (subject-level)
clear all;close all hidden; clc
for o1=1:1 % Documentation (Read me)
    % Some variables in this script have been renamed. If bugs in future
    % scripts, see if changing the names here will help 
    %
    % instruc = 
    %                  log: 
    %         rois_extract: instruc.roi_files
    %         subcon_names: instruc.FLcons
    %        subcon_images: instruc.FLcon_imgs
    %         spec_connums: request.FLcon_nums
    %     req_subcon_names: request.FLcons
    % 
    % d_subroibetas --> d_betas
    %
% INSTRUCTIONS for using this script
%
%   (1) This script extracts ROIs from each subject x contrast image, and outputs 
%           it to the variable 'd_betas'. Betas are extracted for all
%           subjects included in the model (noted in output log.modeldetails.log.subjects)
%
%   (2) Outputs of script: 
%             - 'd_betas' is a cell array, with cells containing betas for all subjects (in order),
%                 corresponding to each ROI (row + 1) and each contrast (col + 1). Headers (ROI & Contrast)
%                 are included, thus displacing the beta vectors by 1. 
%             - 'd_vol.roi_masks':            Cell array, each cell holds the mask for one anatomical ROI
%             - 'd_vol.subcon':                 Cell array, each cell holding the full volume betas for each 
%                                                         subject (row) x contrast (col)
%             - 'instruc.FLcons':    Names of each (first-level, subjects')  contrast
%             - 'instruc.FLcon_imgs':   Contrast image name for each first-level contrast
%
% 	(2) Before running this script, ROIs must be defined within SPM. No need for importing via MarsBar. 
%
%   (3) Specify whether you would like ALL ROIs extracted, or only specific ones. 
%             - All ROIs: ROIs will be extracted for all images with the prefix 'roi', in the model's 
%                 2nd level results folder
%             - Specific ROIs: Only ROIs in the results folder's ROI folder will be extracted,
%                 again with prefix 'roi'
%   

% Other nifty ways of extracting betas
% spm_summarise('beta_1680.img',struct('def','sphere','xyz',[1 2 3]','spec',8),@mean)   %  Extracting from a sphere of 8mm at point 1,2,3
% spm_summarise('beta_1680.img','mask.img',@mean) % Extract mask.img from beta_1680.img (mean)
%   See also: spm_ROI
end

% Requested analysis
request.LoadSpecificROIs_Folder=1; % 1= Load only ROIs that are in the ROI folder, 0=Load all available ROIs for that results model
request.modeltype=1; % 1 = Par, 2= Flex

log.FLthread=' s4Ants';                 % ####### [MODEL DETAILS] #####################
log.AntsType='_Basic';
%
for o=1:1
% log.firstlevelmodel='m_c3_ChoiceFull_OULPEN';
% log.firstlevelmodel= 'm_c7_ChCluster6Full_OULPEN';   log.secondlevelmodel='choice_cluster2x2';
% log.firstlevelmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11'; log.secondlevelmodel='par vBestUnchosen';
% log.firstlevelmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11'; log.secondlevelmodel='par vBestUnchosen_pos';
% log.firstlevelmodel='m_v6e_ChoicevChosenAndposneg2_b01b01'; log.secondlevelmodel='par vBestUnchosen_pos';
% log.firstlevelmodel='m_v9c_vChosenAndposneg2_bpm16bpmi11'; log.secondlevelmodel='par vBestUnchosen_pos';
% log.firstlevelmodel='m_v8c_ChoicesubEVposneg_bpm16bpmi11'; log.secondlevelmodel='par subEVpos';
% log.firstlevelmodel='m_v12c_vModalchoiceOutcomeMagnitude_bpm16bpmi11'; log.secondlevelmodel='par CueValue';
% log.firstlevelmodel='m_v11c_vGambleOutcomeMagnitude_bpm16bpmi11'; log.secondlevelmodel='par CueValue';
% log.firstlevelmodel='m_v13c_OutcomevGambleOutcomeMagnitude_bpm16bpmi11'; log.secondlevelmodel='par CueValue';
% log.firstlevelmodel='m_v14c_OutcomevModalchoiceOutcomeMagnitude_bpm16bpmi11'; log.secondlevelmodel='par CueValue';
% log.firstlevelmodel='m_v15c_VExploreInfoVal_bpm16bpmi11'; log.secondlevelmodel='par Info';
% log.firstlevelmodel='m_v15c_VExploreInfoPE_bpm16bpmi11'; log.secondlevelmodel='par Info';
% log.firstlevelmodel='m_v17c_vChoiceOutcomeMag_bpm16bpmi11'; log.secondlevelmodel='par GamblevExplore';
% log.firstlevelmodel='m_v18c_vChoiceAtOutcomeMag_bpm16bpmi11'; log.secondlevelmodel='par ExploreInfoVal';
% log.firstlevelmodel='m_v21c_ExploreGamInfoOutcomePE_bpm16bpmi11'; log.secondlevelmodel='par ExploreGambletoOutcomePE';
% log.firstlevelmodel='m_v1c_vChoice_bpm16bpmi11'; log.secondlevelmodel='vExplore';
% log.firstlevelmodel='m_v22c_OutcomevChosenOutcomeMagnitude_bpm16bpmi11';log.secondlevelmodel='par CueValue';
% log.firstlevelmodel='m_v22f_OutcomevChosenOutcomeMagnitude_b02b01';log.secondlevelmodel='par CueValue';
end
% log.firstlevelmodel='m_v3g_vChosenAnd_bpji08bpji11';  log.secondlevelmodel='TaskVal';
% log.firstlevelmodel='m_v4g_ChoicevChosenAnd_bpji08bpji11'; log.secondlevelmodel='choice_2x3';
% log.firstlevelmodel='m_v9g_vChosenAndposneg_bpji08bpji11'; log.secondlevelmodel='vBestUnchosenPosNeg';
% log.firstlevelmodel='m_v6g_ChoicevChosenAndposneg2_bpji08bpji11'; log.secondlevelmodel='par vBestUnchosen_pos';
% log.firstlevelmodel='m_v23g_RejectOrvChosenAnd_bpji08bpji11';  log.secondlevelmodel='RejOr_vBU';
% log.firstlevelmodel='m_v27g_ChoiceRejectOrvChosenAnd_bpji08bpji11';  log.secondlevelmodel='RejOr_vBU';
% log.firstlevelmodel='m_c14_Choice'; log.secondlevelmodel='choice_2x3';
log.firstlevelmodel='m_c13_ChoiceFull_ULPEN'; log.secondlevelmodel='choice_2x3';
% log.firstlevelmodel='m_v28g_ChoicePredChoice_bpji08bpji11'; log.secondlevelmodel='choice_2x3';
% log.firstlevelmodel='m_c7_ChCluster6Full_ULPEN'; log.secondlevelmodel='choice_2x2';
% log.firstlevelmodel='m_v1g_vChoice_bpji08bpji11'; log.secondlevelmodel='par vExplore';
% log.firstlevelmodel='m_v25g_RejectOrvGamble_bpji08bpji11';  log.secondlevelmodel='TaskbyRejOr';
% log.firstlevelmodel='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11';  log.secondlevelmodel='vBestUnchosen'; 
% log.firstlevelmodel='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11';choice_2x3';
% log.firstlevelmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';log.secondlevelmodel='TaskVar_2x2';
% log.firstlevelmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';log.secondlevelmodel='TaskVar_2x2';
% log.firstlevelmodel='m_v38g_EVGainLoss_bpji08bpji11';log.secondlevelmodel='EVGainLoss';
% log.firstlevelmodel='m_v40g_vBestvWorst_bpji08bpji11';log.secondlevelmodel='TaskVar_2x2';
% log.firstlevelmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';log.secondlevelmodel='vMargCho x vBUpn';
% log.firstlevelmodel='m_v42g_pLossNTok_bpji08bpji11';log.secondlevelmodel='par pLossNTok';
% log.firstlevelmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11';  log.secondlevelmodel='TaskXChoicemargvcho_2x3';
% log.firstlevelmodel='m_c21_RejExpFull_ULPEN';  log.secondlevelmodel='choice_2x2';
% log.firstlevelmodel='m_c22_ChoiceFull_HLPEN';log.secondlevelmodel='choice_2x3';

for o1=1:1 % unused firstlevels 
    % log.firstlevelmodel='m_c1_CompeteBasic';
    % log.firstlevelmodel='m_c2_CompeteBasicVExplore';
% log.firstlevelmodel='m_c4_CompeteFull_XUPEN';
% log.firstlevelmodel='m_c5_CompeteFullRT_XUVPEN';
% log.firstlevelmodel='m_c6_Cluster4CompeteFull_XUVPEN';
% log.firstlevelmodel='m_c8_Cluster4MovCompeteFull_XUVPEN';
% log.firstlevelmodel='m_c9_Cluster6MovCompeteFull_XUVPEN';
% log.firstlevelmodel='m_c10_Cluster6CompeteRT_XUVPEN';
% log.firstlevelmodel='m_o1_OrthogBasic';
end
% log.secondlevelmodel='choice_2x3';                              % ----------- Second level-----------
% log.secondlevelmodel='vBestUnchosen';
for o1=1:1 % unused second levels

% log.secondlevelmodel='choiceRT_2x3';                              
% log.secondlevelmodel='choice_cluster2x2';
% log.secondlevelmodel='choiceRT_cluster2x2';
% log.secondlevelmodel='ExVEx_2x2';
% log.secondlevelmodel='Identity_1samplettest';
% log.secondlevelmodel='Identity_pairedttest';
% log.secondlevelmodel='RL_2xN_fakefactorial';
end
for o1=1:1 % Unused, flex models) & PPIs
% ######### [Flexible Models] ############################################
% log.firstlevelmodel='f1_TaskChoiceTrialtype';     % ----------- First level-----------
% log.firstlevelmodel='f2_ExploreOr_AllQualCells';
% log.firstlevelmodel='f3_ExploreOr_FixWind4';
% log.firstlevelmodel='f4_ExploreOr_MovWind';
% log.firstlevelmodel='t1_RLvars';
% log.rlvar='Entropy';
% log.secondlevelmodel=[log.rlvar ' TaskAnova'];

% ############ PPIs ###################

%
log.secondlevelthread=[];
% log.secondlevelthread='Comparison PPIs - HPC StrictAnat\HPC StrictAnat Choice';
% seed='Caudate_R';
% log.secondlevelmodel=[ 'compFam_' seed '_psy_cF_Rej-Exp'];
% sd='HPC_L_sac'; 
log.betasuffix=[]; 

end


for o1=1:1 % General settings and specifications
    
    % where.where='/Volumes/PENNYDISK/5 Explore fMRI'; where.experiment_folder='/Volumes/SANDISK/2 EXPLORE Brain data'; where.data_brain=[where.experiment_folder filesep '1 MRI data'];
    % where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
    where.where='E:\Dropbox\SCRIPPS\5 Explore fMRI';  where.experiment_folder='G:\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
    
    
    where.where= 'C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; 
     
%     \2 Second level results s4Ants\m_c14_Choice_Basic\choice_2x3\ROI\c13 battery
    
    % Subjects
    addpath(where.where);addpath([where.where filesep '4 Set up models'])
    log.specificsubjects={}; % Betas are always extracted for all available subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    instruc.log.subjects=log.subjects; instruc.log.n_subjs=log.n_subjs; 
    
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c17';'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_c10'};
%     w.modelsneedingsubselect = cellfun(@(x)[x log.AntsType],  w.modelsneedingsubselect, 'UniformOutput',0);
    if sum(strcmp(log.firstlevelmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
    
    
%     E:\1 Explore fMRI archive\1 First level\p01_GV\2 First level s4Ants
    % Where are FL cons?
     log.firstlevelmodel=[log.firstlevelmodel log.AntsType];
    if sum(strcmp(log.firstlevelmodel(3),{'f';'t'}))==1
        where.data_brain_FL= 'E:\1 Explore fMRI';
    elseif sum(strcmp(log.firstlevelmodel(1:6),{'m_c13g';'m_c14_'}))==1
        where.data_brain_FL= 'E:\1 Explore fMRI archive\1 First level';
    else where.data_brain_FL=[where.experiment_folder filesep '1 Brain data'];
    end
    
    % Change model names
    if sum(strcmp(log.firstlevelmodel(3),{'f';'t'}))==1 && strcmp(log.secondlevelmodel(1:4), 'ppi_')==0    % Flexible model + No PPI #########
        log.FLfolname=['m_' log.firstlevelmodel(1:2) ' Contrasted   ' log.firstlevelmodel];%       Which onsets model goes with the requested first-level (contrast) model? 
    elseif sum(strcmp(log.firstlevelmodel(3),{'f';'t'}))==1 && strcmp(log.secondlevelmodel(1:4), 'ppi_')==1  % Flexible model + PPI #########
        error('Model firstlevel folders not changed yet for flex + ppis')
    elseif strcmp(log.secondlevelmodel(1:4), 'ppi_')==0  && isempty(strfind(log.secondlevelmodel, 'psy_'))==1  % Parametized models + No PPI #########
        log.FLfolname=[log.firstlevelmodel ' Contrasted'];
    elseif strcmp(log.secondlevelmodel(1:4), 'ppi_')==1  || isempty(strfind(log.secondlevelmodel, 'psy_'))==0 % Parametized models + PPI #########        
%         log.FLfolname=[log.firstlevelmodel ' PPIs' filesep
%         log.secondlevelmodel(strfind(log.secondlevelmodel, '         ')+5:length(log.secondlevelmodel))]; % OLD specification
        log.FLfolname=[log.firstlevelmodel ' PPIs' filesep log.secondlevelmodel];
    end
    
    % Settings that don't change much
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    if strcmp(where.where(1),'/')==1; fs='/'; where.marsbar_anatomy='/Users/EleanorL/Documents/MATLAB/spm8/Anatomy masks/marsbar_imported'; else fs='\'; where.marsbar_anatomy='D:\My Documents\MATLAB\spm8\Anatomy masks\marsbar_imported'; end
    where.aPrioriROIs=[where.experiment_folder filesep '3 Checks' filesep '2 A priori ROIs' filesep 'ROI analysis'];
    
    % Model details (Full model, 1st & 2nd levels)
    where.model_res=[where.experiment_folder  filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep log.secondlevelthread log.secondlevelmodel filesep];
    m=load([where.model_res 'SPM.mat']); m=m.SPM;
    if strcmp(log.secondlevelmodel(1:4), 'ppi_')==0 || isempty(strfind(log.secondlevelmodel, 'psy_'))==0 
%         log.modeldetails=load([where.model_res 'details_2ndlevel.mat']);
%         p=load([where.data_brain filesep log.modeldetails.log.subjects{1} filesep '2 First level' log.secondlevelthread filesep log.FLfolname filesep 'SPM.mat']);p=p.SPM;
    elseif strcmp(log.secondlevelmodel(1:4), 'ppi_')==1
        % Details from m
        log.modeldetails.log.subjects=log.w.datalog(2:end,1); disp('PPI model. Assumed all subjects included.');
        log.modeldetails.log.n_subjs=length(log.modeldetails.log.subjects);  
        if size(m.xY.P,1)>log.modeldetails.log.n_subjs; error('PPI 2nd level, assumed only 1 ingoing contrast. BUT NOT TRUE. Which to load?'); end
    end
    
    log.modeldetails.log.subjects = log.subjects; 
    log.modeldetails.log.n_subjs = log.n_subjs; 
    p=load([where.data_brain_FL filesep log.modeldetails.log.subjects{1} filesep '2 First level' log.FLthread filesep log.FLfolname filesep 'SPM.mat']);p=p.SPM;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.experiment_folder])
    disp(' '); disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             Model type:               ' log.secondlevelthread])
    disp(['             First level model:      ' log.firstlevelmodel])
    disp(['                       First level contrast folder:   ' log.FLfolname])
    disp(['             Second level model:  ' log.secondlevelmodel]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% AD HOC SETTINGS * * * #################
% Change where you want to extract the betas from etc, if anything is non-standard

for o=1:1 % Template (ROI locations and names)
    for o=1:1 % Before beh modelling sorted April 15
    
% % c3 MEChoice (all)
% subfol.c3allchoice='1 ME Choice\1 All ROI images';
% rois.c3allchoice={
%     'HPC_L';
%     'HPC_R';
%     'SFG_BA46';
%     'MFG_BA9';
%     'OFC_Mid';
%     'SFG_BA11';
%     'Med_FG';
%     'Cingulate_Pos1';
%     'Cingulate_Pos2';
%     'Precuneus1';
%     'Precuneus2';
%     'PostcentralGyr';
%     'Occipit';
%     'Cereb';
%     };
% 
% % c3 Frontal only
% subfol.c3choicefrontal='1 ME Choice\2a Frontal only';
% rois.c3choicefrontal={
%     'BA46';
%     'BA9';
%     'OFC_mid';
%     'BA11';
%     'MedFG';};
% 
% % c3 MEChoice Amygdala-HPCextend
% subfol.c3choiceamyg='1 ME Choice\2e Amygdala HPC-extension';
% rois.c3choiceamyg={'Amyg_L_hpcextend';'Amyg_R_hpcextend';};
% 
% % c3 MEChoice HPC - Strict
% subfol.c3choiceHPCstrict='3 [Strict] ME Choice masked by Task x Choice';
% rois.c3choiceHPCstrict={'HPC_L_strictchoice';'HPC_R_strictchoice'};
% 
% 
% % c3 Task x Choice HPC - Strict
% subfol.c3choiceHPCstrict='4 [Strict] Task x Choice masked by Choice';
% rois.c3choiceHPCstrict={'HPC_L_strictTcf';'HPC_R_strictTcf'};
% 
% % c3 StricAnat TxC clusters
% subfol.c3hpcstrictanat='5 [StrictAnat] HPC clusters';
% rois.c3hpcstrictanat={
% 'HPC_L_strictanatChoice';
% 'HPC_L_strictanatTcf';
% 'HPC_R_strictanatChoice';
% 'HPC_R_strictanatTcf';
% };
% 
% Battery of c3 clusters (Main ROIs for beta extraction)
rois.c3battery={
    'BA10';
    'BA46';
    'HPC_L_c';
    'HPC_L_sac';
    'HPC_L_satc';
    'HPC_R_c';
    'HPC_R_sac';
    'HPC_R_satc';
    'Precuneus';
    'Striatum_L';
    'Striatum_R';
    'SupMidFG';
    };
subfol.c3battery='c3 ROI battery for beta extraction';
% 
% 
% 
% 
% % c3 MEChoice Striatum
% subfol.c3choicestriatum='1 ME Choice\2b Striatal\Anat clusters';
% rois.c3choicestriatum={
%     'Caudate_L_anat';
%     'Caudate_R_anat';
%     'Putamen_L_anat';
%     };
% % subfol.c3choicestriatum='1 ME Choice\2b Striatal\Anat clusters';
% % rois.c3choicestriatum={
% %     'Caudate_L_raw';
% %     'Caudate_R_raw';
% %     };
% 
% 
% 
% % c7 own-clusters selected
% subfol.c7ownselected='2 Selected c7 rois';
% rois.c7ownselected={'BA10_c';'BA46_EmR';'SFG_med_RmE';'HPC_L_tc'};
% 
% 
% % Battery of v4 (vBestUnchosenNeg) clusters 
% rois.v4neg={
%     'Amygdala_R';
%     'AntCingulate';    
%     'BA10';    
%     'Caudate';    
%     'HPC_aL';    
%     'HPC_aR';    
%     'vmPFC'};
% % subfol.v4neg='Neg selected rois';
% subfol.v4neg='v4c vBU selected rois';
% 
% % Selected from v6c
% rois.v6cneg={
%    'AmyL';'AmyR'; 'HPC_aL';  
%    'HPC_aR';  'HPC_aR2';  'HPC_aR3';  
% };
% subfol.v6cneg='v6c vBUneg selected rois';

    end

    
    % Battery of c13g choice clusters
    rois.c13gbattery={'DLPFC_L';'DLPFC_R';'HPC_L_sc';'HPC_L_stc';'HPC_R_stc';'HPC_R_sc';'Parietal';'OFC_L';'Striatum_R';};
    subfol.c13gbattery='c13g battery';
    
    rois.c14HPC={'HPC_aL_c';'HPC_aL_stc';'HPC_aL_tc';'HPC_aR_c';'HPC_aR_satc';'HPC_aR_stc';'HPC_aR_tc';};
    subfol.c14HPC='c14 HPC';
end
% mod='c7ownselected';   % Which template to load?
% mod='c13gbattery';
% mod='c3battery';
% mod='c14battery';
% mod='v4neg';
% mod='v6cneg'; 

% % ROI folder name (Betas will be saved here) ##############
% eval(['instruc.roi_names=rois.' mod ';'])
% eval(['roi_subfol=subfol.' mod ';']);
% where.roi_foldername=['ROI' filesep roi_subfol filesep];

% Manually define ############
% where.roi_foldername=['ROI' filesep];
% where.roi_foldername=['ROI' filesep 'Anat HPC Amyg' filesep];
% where.roi_foldername=['ROI' filesep 'Subfields binary' filesep];
% where.roi_foldername=['ROI' filesep 'Subfields ref' filesep];
% where.roi_foldername=['ROI' filesep 'Subfields 95' filesep];
% where.roi_foldername=['ROI' filesep 'c13 battery' filesep];
% where.roi_foldername=['ROI' filesep 'c13 battery and amyg' filesep];
% where.roi_foldername=['ROI' filesep 'c13 amyg' filesep];
% where.roi_foldername=['ROI' filesep 'c14 battery' filesep];
% where.roi_foldername=['ROI' filesep 'c14 HPC' filesep];
% where.roi_foldername=['ROI' filesep 'v3g HPC' filesep];
% where.roi_foldername=['ROI' filesep 'v3g hpc amyg' filesep];
% where.roi_foldername=['ROI' filesep 'v28gChoice battery' filesep];
% where.roi_foldername=['ROI' filesep 'v25g HPC' filesep];
% where.roi_foldername=['ROI' filesep 'v23g' filesep];
% where.roi_foldername=['ROI' filesep 'c20g Choice' filesep];
% where.roi_foldername=['ROI' filesep 'Subfield Choice Wrong' filesep];
% where.roi_foldername=['ROI' filesep 'Subfield fxn masked' filesep]; 
% where.roi_foldername=['ROI' filesep 'v3g Subfields ref masked' filesep]; 
% where.roi_foldername=['ROI' filesep 'v3g HPC' filesep]; 
% where.roi_foldername=['ROI' filesep 'HPC assorted' filesep]; 
% where.roi_foldername=['ROI' filesep 'v46g' filesep]; 
where.roi_foldername=['ROI' filesep 'HPC mec tc conjunc 05unc' filesep];


% instruc.roi_names=cellfun(@(x)x(1:end-4), instruc.roi_files, 'UniformOutput',0) ;  % Use file names (must run file first!!)
instruc.roi_names=[];  % Empty to use roi actual names

%% (1) Extract anatomic ROIs and volume betas 

% Load all ROIs (n-ary or binary masks)
disp('Loading ROI images ########################')
% where.model_res=[where.experiment_folder filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep log.secondlevelthread filesep log.secondlevelmodel filesep];
where.model_res=[where.experiment_folder filesep '2 Second level results' log.FLthread filesep log.firstlevelmodel filesep log.secondlevelmodel filesep];
if request.LoadSpecificROIs_Folder==1; where.model_resrois=[where.model_res where.roi_foldername filesep];  % Read only the ROIs in the ROI folder
else where.model_resrois=[where.model_res filesep];
end
instruc.roi_files=cellstr(spm_select('List', where.model_resrois, '.*.nii$'));   % Change format of files here!! 
if isempty(char(instruc.roi_files))==1;  instruc.roi_files=[]; end
d_vol.roi_masks=cell(size(instruc.roi_files,1),1);
for i=1:size(instruc.roi_files,1)
    d_vol.roi_masks{i}=spm_read_vols(spm_vol([where.model_resrois instruc.roi_files{i}]));
end

% Interface
disp('ROIs to extract:    ( !! CANNOT BE EMPTY !! )'); disp(char(instruc.roi_files)); input('OK?   '); disp('Names:'); disp(instruc.roi_names)

%% (2) Extract subjects' volume (contrast) betas & calculate mean

% Certain first-level contrast images only? ########### *** EDIT HERE *** ###########
disp('Available FL contrasts in this model:  ' ); disp([num2cell(1:size(p.xCon,2))' cellstr(char(p.xCon.name))])
request.FLcon_nums=[];
% request.FLcon_nums=[1:6]; 
if isempty(request.FLcon_nums)==0; disp('Requested contrasts:  '); disp(char(p.xCon(request.FLcon_nums).name));  else disp('Extracting betas from ALL FL contrasts');  end; input('Proceed?    '); 
    
for o1=1:1
    % Load all contrast images from all subjects (d_vol.subcon, row=sub, col=con image)
    disp('Loading subjects'' contrast images ########################')
    instruc.FLcons=cell(size(p.xCon,2),1); instruc.FLcon_imgs=cell(size(p.xCon,2),1);
    for i=1:length(instruc.FLcons) % Get names and available con images (sample 1st subject)
        instruc.FLcons{i}=p.xCon(i).name; instruc.FLcon_imgs{i}=p.xCon(i).Vcon.fname;
    end
  
    % Load first-level contrasts (see above for specification: extract from all or from some only?)
    if isempty(request.FLcon_nums)==0
        instruc.FLcons=instruc.FLcons(request.FLcon_nums);
        instruc.FLcon_imgs=instruc.FLcon_imgs(request.FLcon_nums);
    end
    
    % Extract
    d_vol.subcon=cell(log.modeldetails.log.n_subjs, length(instruc.FLcons));
    for s=1:log.modeldetails.log.n_subjs % Load volume betas for each contrast
        disp(['Subject ' num2str(s) '  -  ' log.modeldetails.log.subjects{s} ' --------------------------------'])
        ws.wheremod1st=[where.data_brain_FL filesep log.modeldetails.log.subjects{s} filesep '2 First level' log.FLthread filesep log.FLfolname filesep];
        disp(ws.wheremod1st)
        for c=1:length(instruc.FLcons)
            disp(['Contrast no. ' num2str(c) '  -  ' instruc.FLcons{c}]);
            d_vol.subcon{s,c}=spm_read_vols(spm_vol([ws.wheremod1st instruc.FLcon_imgs{c}]));
            %         v=spm_vol([ws.wheremod1st instruc.FLcon_imgs{c}]);
            %         d_vol.subcon{s,c}=spm_read_vols(v);
        end
    end
    
    % Subject level (d_betas, col+1=contrast image, row+1=roi)
    disp('Calculating ROI betas for each ROI x Subject x Contrast #################')
    d_betas=cell(length(instruc.roi_files)+1, length(instruc.FLcons)+1);
    for r=1:length(d_vol.roi_masks)
        disp(['ROI no. ' num2str(r) ': ' instruc.roi_files{r} ' -------------------------'])
        d_betas{r+1,1}=instruc.roi_files{r};
        for c=1:length(instruc.FLcons)
            disp(['Contrast no. ' num2str(c) ': '  instruc.FLcons{c} '----'])
            if r==1; d_betas{1,c+1}=instruc.FLcons{c}; end % Title
            
            % Read each subject's ROI betas
            wc.s=nan*zeros(log.modeldetails.log.n_subjs,1);
            for s=1:log.modeldetails.log.n_subjs
                disp(['Subject ' num2str(s) ': ' log.modeldetails.log.subjects{s}])
                wc.s(s)=mean(d_vol.subcon{s,c}(d_vol.roi_masks{r}~=0));
            end
            d_betas{r+1, c+1}=wc.s;
            wc=[];
        end
    end
end

% Clear for workspace 
d_vol=[];

%% (4) Export to txt (for SPSS?)

request.export2excel=1; request.calculateandprint_choicedifferencescores=0; if request.calculateandprint_choicedifferencescores; disp('Choice difference scores are requested!'); end

% Write to text file
if request.export2excel==1
    
    % ########### *** EDIT HERE *** ###########
    if isempty(instruc.roi_names); 
        % If crashing here, likely you have no ROIs
        instruc.roi_names=cellfun(@(x)x(1:strfind(x, '.')-1),  instruc.roi_files, 'UniformOutput',0); 
    end
    text.roi_names=instruc.roi_names; text.roi_names=cellfun(@(x)[x log.betasuffix], text.roi_names, 'UniformOutput', 0);
    request.FLcons=[]; if isempty(request.FLcons);  request.FLcons=instruc.FLcons; end
    
    % Confirm settings for printout
    disp('#### Printing data to excel file ########################################################')
    disp('(1) Short names for ROIs: '); disp([text.roi_names  repmat({'       ' }, length(instruc.roi_files), 1) instruc.roi_files]); disp(' ');
    disp('(2) FL contrast names:     ');disp(request.FLcons); disp(' '); input('Continue with these names?   '); disp(' '); disp(' -------------------------------' )
    
    % ########### *** EDIT HERE *** ###########
    
    t_betas= [[{'Subject'}; log.modeldetails.log.subjects] cell(log.modeldetails.log.n_subjs+1, length(instruc.roi_files)*length(request.FLcons)); ]; k=1;
    for r=1:size(instruc.roi_files,1)   % Print subject roi x contrast beta values 
        for n=1: length(request.FLcons)
            c=find(strcmp(instruc.FLcons, request.FLcons{n})); % speccific contrasts
            if isempty(c)==1; error(['Invalid first-level name requested!   ('  request.FLcons{n} ')']); end
            %
            t_betas{1,k+1}=[text.roi_names{r} '-' instruc.FLcons{c}]; % Title
%             disp(t_betas{1,k+1})
            for s=1:log.modeldetails.log.n_subjs
                t_betas{s+1, k+1}=d_betas{r+1, c+1}(s);
            end
            k=k+1;
        end
    end
    for o1=1:1  % Difference scores 
    
    % (if requested) Print difference scores for each ROI: cF_Reject-cF_Explore,  cF_Reject-ct_Bomb
    if request.calculateandprint_choicedifferencescores
        disp('Generating & printing difference scores ###########')
        
        % Instructions - specific difference-scores to compute, in order of request
        inclusterprefix=[]; 
        diffcontrasts={   {'cF_Accept'; 'ct_NoBomb'};       % Cross-task comparisons
                                  {'cF_Reject'; 'ct_Bomb'};
                                  {'cF_Explore'; 'ct_Explore'};
                                  {'cF_Reject'; 'cF_Explore'};      % Within task, compare choice
                                  {'cF_Reject'; 'cF_Accept'};
                                  {'cF_Accept'; 'cF_Explore'};
                                  {'ct_Bomb'; 'ct_Explore'};
                                  {'ct_Bomb'; 'ct_NoBomb'};
                                  {'ct_NoBomb'; 'ct_Explore'};
                              };
        
        inclusterprefix='in_'; 
        diffcontrasts={   {'cF_Reject'; 'ct_Bomb'};       % Cross-task comparisons
                                  {'cF_Explore'; 'ct_Explore'};
                                  {'cF_Reject'; 'cF_Explore'};      % Within task, compare choice
                                  {'ct_Bomb'; 'ct_Explore'};
                              };
        
        
        for r=1:size(instruc.roi_files,1)
            
            % Locate cell contrasts
%             c.cF_Accept=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Accept'])); if isempty(c.cF_Accept)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Reject=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Reject'])); if isempty(c.cF_Reject)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Explore=find(strcmp(instruc.FLcons, [inclusterprefix 'cF_Explore']));  if isempty(c.cF_Explore)==1; error('cannot find requested contrast, difference scores'); end
            %
%             c.ct_NoBomb=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_NoBomb'])); if isempty(c.ct_NoBomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Bomb=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_Bomb'])); if isempty(c.ct_Bomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Explore=find(strcmp(instruc.FLcons, [inclusterprefix 'ct_Explore']));  if isempty(c.ct_Explore)==1; error('cannot find requested contrast, difference scores'); end
            
            
            
            % Apply contrasts (instructed above in diffcontrasts)
            for jj=1:size(diffcontrasts,1)
                cons=diffcontrasts{jj};
                t_betas{1,k+1}=[text.roi_names{r} '-' inclusterprefix cons{1} '_minus_' inclusterprefix cons{2}]; disp(t_betas{1,k+1})
                eval(['a{1}=c.' cons{1} ';']); eval(['a{2}=c.' cons{2} ';']);
                for s=1:log.modeldetails.log.n_subjs
                    t_betas{s+1, k+1}=d_betas{r+1, a{1}+1}(s)-d_betas{r+1, a{2}+1}(s);
                end
                k=k+1;
            end
            
            
            
        end
    end
    end
%     t_betas=[t_betas [{'BORDER_BEHAVIOUR'}; cell(size(t_betas,1)-1,1)]];
%     xlswrite(['('  date ') Extracted betas' log.betasuffix '.xlsx'], t_betas);  movefile(['('  date ') Extracted betas' log.betasuffix '.xlsx'],  [where.model_resrois  '('  date ') Extracted betas' log.betasuffix '.xlsx'])   % Write file without the behavioural profile
    
    % Request behavioural profiles to append
    request.behfile=[where.where filesep '1 Behavioural data' filesep 'Group behaviour files' filesep 'Behavioural profile for correlation.xlsx'];
    if sum(strcmp(log.firstlevelmodel(1:5), {'m_c7_'; 'm_c17'}))==1 ; request.beh_excelsheets={'BehInhibition c7';'Explore c7'};
    else request.beh_excelsheets={'BehInhibition';'Explore'};
    end
    disp('Behavioural profile to append (new document written):'); disp(['     '  request.behfile]); disp( ' '); disp('Specific excel sheets (from beh file) to append:'); disp(request.beh_excelsheets); disp(' '); input('Continue?   '); disp(' ')
    
    % Appendity append
    for i=1:length(request.beh_excelsheets);
        t_betas=[t_betas [{['BEH_'  request.beh_excelsheets{i}]}; cell(size(t_betas,1)-1,1)]];
        [n t r]=xlsread(request.behfile,  request.beh_excelsheets{i});
        n=num2cell(n); n(isnan(cell2mat(n(:))))=repmat({[]}, sum(isnan(cell2mat(n(:)))),1); n=reshape(n, size(t)-1);
        [t_betas wb.headers wb.subs wb.n_subs excl] =f_aligntables(t_betas, [t(:,1) [t(1,2:end); n]]);
        
        if isempty(excl.table1) + isempty(excl.table2)~=2
            disp(['Subjects mismatch while appending beta table (#1) from profile ' request.beh_excelsheets{i} ' (#2)'])            
            disp('                 (subjects that are not present in both tables are excluded entirely)'); disp(excl); input('Continue?   ');
        end
    end
    [st{1} msg{2}]=xlswrite(['('  date ') Extracted betas' log.betasuffix '.xlsx'], t_betas);  [st{2} msg{2}]=movefile(['('  date ') Extracted betas' log.betasuffix '.xlsx'],  [where.model_resrois  '('  date ') Extracted betas' log.betasuffix '.xlsx']);   % Write file without the behavioural profile
    disp(['        Excel file written:                      ' num2str(st{1})])
    disp(['        File moved to correct location:  ' num2str(st{2})])
    
end

%%

error('Done :)')

b_state=cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), 'st.State'))));
b_trait=cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), 'st.Trait'))));
% b_corres=cell2mat(length(instruc.roi_names)


for r=1:length(instruc.roi_names)
    [wr.r wr.p]=corr(b_state,   cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Reject']))))- cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Explore'])))));
    [wr.r wr.p]=corr(b_trait,   cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Reject']))))- cell2mat(t_betas(2:end, find(strcmp(t_betas(1,:), [instruc.roi_names{r} '-cF_Explore'])))));
    
    
end



openvar t_betas



t_betas(1,:)











