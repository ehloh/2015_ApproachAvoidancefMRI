% Extract ROI betas using SPM functions (subject-level)
clear all;close all hidden; clc

for o1=1:1 % Documentation (Read me)
    
% INSTRUCTIONS for using this script
%
%   (1) This script extracts ROIs from each subject x contrast image, and outputs 
%           it to the variable 'd_subroibetas'. Betas are extracted for all
%           subjects included in the model (noted in output log.modeldetails.log.subjects)
%
%   (2) Outputs of script: 
%             - 'd_subroibetas' is a cell array, with cells containing betas for all subjects (in order),
%                 corresponding to each ROI (row + 1) and each contrast (col + 1). Headers (ROI & Contrast)
%                 are included, thus displacing the beta vectors by 1. 
%             - 'd_vol.roi_masks':            Cell array, each cell holds the mask for one anatomical ROI
%             - 'd_vol.subcon':                 Cell array, each cell holding the full volume betas for each 
%                                                         subject (row) x contrast (col)
%             - 'instruc.subcon_names':    Names of each (first-level, subjects')  contrast
%             - 'instruc.subcon_images':   Contrast image name for each first-level contrast
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

log.firstlevelmodel='m_c3_ChoiceFull_OULPEN';                 % ####### [MODEL DETAILS] #####################
% log.firstlevelmodel='m_c7_Cluster6CompeteFull_OULPEN';
% log.firstlevelmodel='m_v4_ChoicevChosenAnd_bpm16bpm11';
% log.firstlevelmodel='m_v6_ChoicevChosenAndposneg2_bpm16bpm11';
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
log.secondlevelmodel='choice_2x3';                              % ----------- Second level-----------
% log.secondlevelmodel='choice_cluster2x2';
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
    where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];

    % Subjects
    addpath(where.where);addpath([where.where filesep '4 Set up models'])
    log.specificsubjects={}; % Betas are always extracted for all available subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    instruc.log.subjects=log.subjects; instruc.log.n_subjs=log.n_subjs; 
    
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_Cluster4CompeteFull_XUVPEN';'m_c7_Cluster6CompeteFull_XUVPEN';'m_c8_Cluster4MovCompeteFull_XUVPEN';'m_c9_Cluster6MovCompeteFull_XUVPEN';'m_c10_Cluster6CompeteRT_XUVPEN'};
    if sum(strcmp(log.firstlevelmodel, w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
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
    where.model_res=[where.experiment_folder  filesep '2 Second level results' filesep log.firstlevelmodel filesep log.secondlevelthread log.secondlevelmodel filesep];
    m=load([where.model_res 'SPM.mat']); m=m.SPM;
    if strcmp(log.secondlevelmodel(1:4), 'ppi_')==0 || isempty(strfind(log.secondlevelmodel, 'psy_'))==0 
        log.modeldetails=load([where.model_res 'details_2ndlevel.mat']);
%         p=load([where.data_brain filesep log.modeldetails.log.subjects{1} filesep '2 First level' log.secondlevelthread filesep log.FLfolname filesep 'SPM.mat']);p=p.SPM;
    elseif strcmp(log.secondlevelmodel(1:4), 'ppi_')==1
        % Details from m
        log.modeldetails.log.subjects=log.w.datalog(2:end,1); disp('PPI model. Assumed all subjects included.');
        log.modeldetails.log.n_subjs=length(log.modeldetails.log.subjects);  
        if size(m.xY.P,1)>log.modeldetails.log.n_subjs; error('PPI 2nd level, assumed only 1 ingoing contrast. BUT NOT TRUE. Which to load?'); end
    end
    p=load([where.data_brain filesep log.modeldetails.log.subjects{1} filesep '2 First level' filesep log.FLfolname filesep 'SPM.mat']);p=p.SPM;
    
    % Log
    diary([where.data_brain filesep 'SPM logs' filesep 'Log '  mfilename log.firstlevelmodel '-' log.secondlevelmodel ' (' date ')  '])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
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

for o1=1:1 % Template (ROI locations and names)
    
% c3 MEChoice (all)
subfol.c3allchoice='1 ME Choice\1 All ROI images';
rois.c3allchoice={
    'HPC_L';
    'HPC_R';
    'SFG_BA46';
    'MFG_BA9';
    'OFC_Mid';
    'SFG_BA11';
    'Med_FG';
    'Cingulate_Pos1';
    'Cingulate_Pos2';
    'Precuneus1';
    'Precuneus2';
    'PostcentralGyr';
    'Occipit';
    'Cereb';
    };

% c3 Frontal only
subfol.c3choicefrontal='1 ME Choice\2a Frontal only';
rois.c3choicefrontal={
    'BA46';
    'BA9';
    'OFC_mid';
    'BA11';
    'MedFG';};

% c3 MEChoice Amygdala-HPCextend
subfol.c3choiceamyg='1 ME Choice\2e Amygdala HPC-extension';
rois.c3choiceamyg={'Amyg_L_hpcextend';'Amyg_R_hpcextend';};

% c3 MEChoice HPC - Strict
subfol.c3choiceHPCstrict='3 [Strict] ME Choice masked by Task x Choice';
rois.c3choiceHPCstrict={'HPC_L_strictchoice';'HPC_R_strictchoice'};


% c3 Task x Choice HPC - Strict
subfol.c3choiceHPCstrict='4 [Strict] Task x Choice masked by Choice';
rois.c3choiceHPCstrict={'HPC_L_strictTcf';'HPC_R_strictTcf'};

% c3 StricAnat TxC clusters
subfol.c3hpcstrictanat='5 [StrictAnat] HPC clusters';
rois.c3hpcstrictanat={
'HPC_L_strictanatChoice';
'HPC_L_strictanatTcf';
'HPC_R_strictanatChoice';
'HPC_R_strictanatTcf';
};

% Battery of c3 clusters (Main ROIs for beta extraction)
rois.c3battery={
    'BA11';
    'BA46';
    'Caudate_L_anat';
    'Caudate_R_anat';
    'HPC_L_sac';
    'HPC_L_satc';
    'HPC_R_sac';
    'HPC_R_satc';
    'Putamen_L_anat';
    };
subfol.c3battery='c3 ROI battery for beta extraction';




% c3 MEChoice Striatum
subfol.c3choicestriatum='1 ME Choice\2b Striatal\Anat clusters';
rois.c3choicestriatum={
    'Caudate_L_anat';
    'Caudate_R_anat';
    'Putamen_L_anat';
    };
% subfol.c3choicestriatum='1 ME Choice\2b Striatal\Anat clusters';
% rois.c3choicestriatum={
%     'Caudate_L_raw';
%     'Caudate_R_raw';
%     };



% c7 own-clusters selected
subfol.c7ownselected='2 Selected c7 rois';
rois.c7ownselected={'BA10_c';'BA46_EmR';'SFG_med_RmE';'HPC_L_tc'};

end
% mod='c7ownselected';   % Which template to load?
 
% % ROI folder name (Betas will be saved here) ##############
% eval(['namesofrois=rois.' mod ';'])
% eval(['roi_subfol=subfol.' mod ';'];
% where.roi_foldername=['ROI' filesep roi_subfol];

% Manually define ############
where.roi_foldername=['ROI' filesep 'c3 ROI battery' filesep 'Test'];
namesofrois={'Amygdala_saC';'BA10_C'};

%% (1) Extract anatomic ROIs and volume betas 

% Load all (anatomical) ROIs (n-ary or binary masks)
disp('Loading anatomical ROIs ########################')
where.model_res=[where.experiment_folder filesep '2 Second level results' filesep log.firstlevelmodel filesep log.secondlevelthread filesep log.secondlevelmodel filesep];
if request.LoadSpecificROIs_Folder==1; where.model_resrois=[where.model_res where.roi_foldername filesep];  % Read only the ROIs in the ROI folder
else where.model_resrois=[where.model_res filesep];
end
instruc.rois_extract=cellstr(spm_select('List', where.model_resrois, '.*.nii$'));   % Change format of files here!! 
if isempty(char(instruc.rois_extract))==1;  instruc.rois_extract=[]; end
d_vol.roi_masks=cell(size(instruc.rois_extract,1),1);
for i=1:size(instruc.rois_extract,1)
    d_vol.roi_masks{i}=spm_read_vols(spm_vol([where.model_resrois instruc.rois_extract{i}]));
end

% Interface
disp('ROIs to extract:'); disp(char(instruc.rois_extract)); input('OK?   '); disp('Names:'); disp(namesofrois)

%% (2) Extract subjects' volume (contrast) betas & calculate mean

% Certain first-level contrast images only? ########### *** EDIT HERE *** ###########
disp('Available FL contrasts in this model:  ' ); disp([num2cell(1:size(p.xCon,2))' cellstr(char(p.xCon.name))])
instruc.spec_connums=[];
instruc.spec_connums=1:2; 
if isempty(instruc.spec_connums)==0; disp('Requested contrasts:  '); disp(char(p.xCon(instruc.spec_connums).name));  else disp('Extracting betas from ALL FL contrasts');  end; input('Proceed?    '); 
    
for o1=1:1
    % Load all contrast images from all subjects (d_vol.subcon, row=sub, col=con image)
    disp('Loading subjects'' contrast images ########################')
    instruc.subcon_names=cell(size(p.xCon,2),1); instruc.subcon_images=cell(size(p.xCon,2),1);
    for i=1:length(instruc.subcon_names) % Get names and available con images (sample 1st subject)
        instruc.subcon_names{i}=p.xCon(i).name; instruc.subcon_images{i}=p.xCon(i).Vcon.fname;
    end
  
    % Load first-level contrasts (see above for specification: extract from all or from some only?)
    if isempty(instruc.spec_connums)==0
        instruc.subcon_names=instruc.subcon_names(instruc.spec_connums);
        instruc.subcon_images=instruc.subcon_images(instruc.spec_connums);
        disp('Subset of contrast images only. Selected contrast images:'); disp(instruc.subcon_names); input('OK?   ');
    else disp('Extract betas from all contrast images:'); disp(instruc.subcon_names); input('OK?   ');
    end
    
    % Extract
    d_vol.subcon=cell(log.modeldetails.log.n_subjs, length(instruc.subcon_names));
    for s=1:log.modeldetails.log.n_subjs % Load volume betas for each contrast
        disp(['Subject ' num2str(s) '  -  ' log.modeldetails.log.subjects{s} ' --------------------------------'])
        ws.wheremod1st=[where.data_brain filesep log.modeldetails.log.subjects{s} filesep '2 First level' filesep log.FLfolname filesep];
        disp(ws.wheremod1st)
        for c=1:length(instruc.subcon_names)
            disp(['Contrast no. ' num2str(c) '  -  ' instruc.subcon_names{c}]);
            d_vol.subcon{s,c}=spm_read_vols(spm_vol([ws.wheremod1st instruc.subcon_images{c}]));
            %         v=spm_vol([ws.wheremod1st instruc.subcon_images{c}]);
            %         d_vol.subcon{s,c}=spm_read_vols(v);
        end
    end
    
    % Subject level (d_subroibetas, col+1=contrast image, row+1=roi)
    disp('Calculating ROI betas for each ROI x Subject x Contrast #################')
    d_subroibetas=cell(length(instruc.rois_extract)+1, length(instruc.subcon_names)+1);
    for r=1:length(d_vol.roi_masks)
        disp(['ROI no. ' num2str(r) ': ' instruc.rois_extract{r} ' -------------------------'])
        d_subroibetas{r+1,1}=instruc.rois_extract{r};
        for c=1:length(instruc.subcon_names)
            disp(['Contrast no. ' num2str(c) ': '  instruc.subcon_names{c} '----'])
            if r==1; d_subroibetas{1,c+1}=instruc.subcon_names{c}; end % Title
            
            % Read each subject's ROI betas
            wc.s=nan*zeros(log.modeldetails.log.n_subjs,1);
            for s=1:log.modeldetails.log.n_subjs
                disp(['Subject ' num2str(s) ': ' log.modeldetails.log.subjects{s}])
                wc.s(s)=mean(d_vol.subcon{s,c}(d_vol.roi_masks{r}~=0));
            end
            d_subroibetas{r+1, c+1}=wc.s;
            wc=[];
        end
    end
end

%% (4) Export to txt (for SPSS?)

request.export2text=1; request.calculateandprint_choicedifferencescores=0;
if request.calculateandprint_choicedifferencescores; disp('Choice difference scores are requested!'); end
disp('ROIs:'); disp(instruc.rois_extract); disp('All subject contrast images:'); disp(instruc.subcon_names)
input('Enter manual inputs (ROI names, Contrast names) ?')


% ########### *** EDIT HERE *** ###########
text.roi_names=namesofrois; text.roi_names=cellfun(@(x)[x log.betasuffix], text.roi_names, 'UniformOutput', 0);

% (B) Certain first-level contrasts only?
instruc.req_subcon_names=[]; if isempty(instruc.req_subcon_names);  instruc.req_subcon_names=instruc.subcon_names; end

% ########### *** EDIT HERE *** ###########



% Write to text file
if request.export2text==1
    disp('ROI names'); disp(text.roi_names); disp('Contrast names'); disp(instruc.req_subcon_names); input('ROI names & Con names ok?');
    disp('Writing to excel file for export to SPSS -------------------- ')
    d_subjectbetas=cell(length(instruc.rois_extract)*length(instruc.req_subcon_names)+1); k=1;
    
    error
    
    % Subject names 
    for s=1:log.modeldetails.log.n_subjs  
        d_subjectbetas{s+1,1}=log.modeldetails.log.subjects{s};
    end
    
    % Print subject roi x contrast beta values
    for r=1:size(instruc.rois_extract,1) 
        for n=1: length(instruc.req_subcon_names)
            c=find(strcmp(instruc.subcon_names, instruc.req_subcon_names{n})); % speccific contrasts
            if isempty(c)==1; error(['Invalid first-level name requested!   ('  instruc.req_subcon_names{n} ')']); end
            %
            d_subjectbetas{1,k+1}=[text.roi_names{r} '-' instruc.subcon_names{c}]; % Title
            disp(d_subjectbetas{1,k+1})
            for s=1:log.modeldetails.log.n_subjs
                d_subjectbetas{s+1, k+1}=d_subroibetas{r+1, c+1}(s);
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
        
        
        for r=1:size(instruc.rois_extract,1)
            
            % Locate cell contrasts
%             c.cF_Accept=find(strcmp(instruc.subcon_names, [inclusterprefix 'cF_Accept'])); if isempty(c.cF_Accept)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Reject=find(strcmp(instruc.subcon_names, [inclusterprefix 'cF_Reject'])); if isempty(c.cF_Reject)==1; error('cannot find requested contrast, difference scores'); end
            c.cF_Explore=find(strcmp(instruc.subcon_names, [inclusterprefix 'cF_Explore']));  if isempty(c.cF_Explore)==1; error('cannot find requested contrast, difference scores'); end
            %
%             c.ct_NoBomb=find(strcmp(instruc.subcon_names, [inclusterprefix 'ct_NoBomb'])); if isempty(c.ct_NoBomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Bomb=find(strcmp(instruc.subcon_names, [inclusterprefix 'ct_Bomb'])); if isempty(c.ct_Bomb)==1; error('cannot find requested contrast, difference scores'); end
            c.ct_Explore=find(strcmp(instruc.subcon_names, [inclusterprefix 'ct_Explore']));  if isempty(c.ct_Explore)==1; error('cannot find requested contrast, difference scores'); end
            
            
            
            % Apply contrasts (instructed above in diffcontrasts)
            for jj=1:size(diffcontrasts,1)
                cons=diffcontrasts{jj};
                d_subjectbetas{1,k+1}=[text.roi_names{r} '-' inclusterprefix cons{1} '_minus_' inclusterprefix cons{2}]; disp(d_subjectbetas{1,k+1})
                eval(['a{1}=c.' cons{1} ';']); eval(['a{2}=c.' cons{2} ';']);
                for s=1:log.modeldetails.log.n_subjs
                    d_subjectbetas{s+1, k+1}=d_subroibetas{r+1, a{1}+1}(s)-d_subroibetas{r+1, a{2}+1}(s);
                end
                k=k+1;
            end
            
            
            
        end
    end
    end
    
    
    
    
    
    
    
    % Print to txt
    cd(where.model_resrois); cd ..
    printok=print2txt(where.model_resrois, ['('  date ') Extracted betas' log.betasuffix], d_subjectbetas);
    disp(printok)
    
    
    openvar d_subjectbetas
end
