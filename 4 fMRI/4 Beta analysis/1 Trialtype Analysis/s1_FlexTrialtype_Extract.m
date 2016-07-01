% Extract ROI betas using SPM functions (subject-level) for FLEXIBLE models only
%       Flex models: Task x Choice x EnvThreat x NTokens i.e. empty cells for some subjects
%       Trialtype models: Task x EnvThreat x NTokens 
%       Extract betas for all available contrasts (or specified numbers). Identity of each contrast varies from subject to subject.
clear all;close all hidden; clc
% clear all;clc

% Request
log.specificsubjects={};
%
log.AntsType='_Basic';

% log.FLmodel='t1_1_Trialtype';  % ##### WHICH FL MODEL #########
log.FLmodel='t2_1_TrialtypeNc';

for o1=1:1 % General settings and specifications 
    
    % Folders
    w.w=pwd; if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SCRIPPS/5 Explore fMRI'; where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data';  where.data_brain=[where.experiment_folder filesep '1 Brain data'];  where.modellingdata=     '/Users/EleanorL/Dropbox/SCRIPPS/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs';   where.modfol='/Users/EleanorL/Dropbox/SCRIPPS/4 Explore experiment/3 Analysis/4 Fit computational models';
    else where.where='D:\Dropbox\SCRIPPS\5 Explore fMRI'; 
        where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]';  where.modfol='D:\Dropbox\SCRIPPS\4 Explore experiment\3 Analysis\4 Fit computational models';
        where.data_brain='E:\1 Explore fMRI'; addpath(where.where);  where.modellingdata='D:\Dropbox\SCRIPPS\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs';   end;  addpath(where.where); addpath(where.modfol)
    if isdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType])==0; mkdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType]); end
    log.SLfol=[where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType filesep]; if isdir(log.SLfol)==0; mkdir(log.SLfol); end
    log.FLfol=['m_' log.FLmodel log.AntsType  ' Contrasted'];
    path(pathdef), addpath(where.where)
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    request.specificsubjects=log.specificsubjects;
    
    %
    if isempty(strfind(log.FLmodel, 'Chunk'))==1; log.nlevels=6;  log.trialchar={'t';[]};
    else log.nlevels=3;  log.trialchar={'e';'n'};
    end
    if strcmp(log.FLmodel(1), 'f')==1; log.nchoices=3; else log.nchoices=1; end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end; disp(' ');
    disp(['             First level model:              ' log.FLmodel])
    disp(['             No. of EnvThreat/NTokens levels:       ' num2str(log.nlevels)]); 
    disp(['             Cells divided by how many choices:     ' num2str(log.nchoices)]); disp(' ')
    disp(['             FL folder:      ' log.FLfol])
    disp(['             SL folder:      ' log.SLfol]); disp(' ');
%     input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% (1) ROI settings + Preprocessing 
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per typical plot)

% Where are the ROI images/beta file  ############
log.roi_fol=[log.SLfol 'Anat HPC Amyg'];
% log.roi_fol=[log.SLfol 'Subfields 95'];
% log.roi_fol=[log.SLfol 'c14 battery'];
% log.roi_fol=[log.SLfol 'c13 battery'];
% log.roi_fol=[log.SLfol 'c13 amyg'];
% log.roi_fol=[log.SLfol 'v3g HPC'];
% log.roi_fol=[log.SLfol 'v28gChoice battery'];
% log.roi_fol=[log.SLfol 'v26e pvBest'];
% log.roi_fol=[log.SLfol 'v28g HPC c14 masking'];
% log.roi_fol=[log.SLfol 'c14 HPC'];
% log.roi_fol=[log.SLfol 'v25g HPC'];


% ROI details (List of ROIs + requested order)
for o=1:1  % Battery details
    for o1=1:1 % From before beh modelling done (April 2015)
        % c3 Choice ###########
        log.roilist.c3battery={
            'AmygdalaPeri_sC';  'Amygdala_saC';
            'BA10_C';  'BA46_L_C';
            'HPC_aL_CTC';  'HPC_aL_c5RTcFRejMOthers';  'HPC_aL_sC';  'HPC_aL_sTC';
            'HPC_aR_Amyg_sC';  'HPC_aR_CTC';  'HPC_aR_sTC';  'HPC_aR_saC';
            'MidFG_R_C';  'Striatum_L_C';  'Striatum_R_C';
            };
        log.requestrois.c3battery={
            'BA10_C';  'BA46_L_C';  'MidFG_R_C';  % Frontal-striatal
            'Striatum_L_C';  'Striatum_R_C';
            'HPC_aL_sTC';  'HPC_aR_sTC';  'HPC_aL_c5RTcFRejMOthers';  % HPC (clear clusters)
            'HPC_aL_sC';  'HPC_aR_saC';
            'HPC_aL_CTC'; 'HPC_aR_CTC';    % HPC complex
            'HPC_aR_Amyg_sC';  'Amygdala_saC';  'AmygdalaPeri_sC';
            };
        
        % Choice anat ###########
        log.roilist.choiceanat={'Amygdala_L';  'Amygdala_R';
            'BA10_L';  'BA10_R'; 'BA46_L';'BA46_R'
            'Caudate_L';  'Caudate_R';
            'Cingulum_Ant_L';  'Cingulum_Ant_R';
            'HPC_aL';  'HPC_aR';   'Hippocampus_L';  'Hippocampus_R';
            'Insula_L';  'Insula_R';
            'NAccCore_L';  'NAccCore_R'; 'NAccShell_L'; 'NAccShell_R'; 'NAcc_L'; 'NAcc_R';
            'Putamen_L'; 'Putamen_R';
            'SNVTA_L'; 'SNVTA_R'; 'Striatum_L'; 'Striatum_R';
            };
        
        log.requestrois.choiceanat={
            'BA10_L';  'BA10_R'; 'BA46_L';'BA46_R'; 'Cingulum_Ant_L';  'Cingulum_Ant_R';               % Cortex
            'Striatum_L'; 'Striatum_R'; 'Caudate_L';  'Caudate_R'; 'Putamen_L'; 'Putamen_R'; % Striatum
            'NAccCore_L';  'NAccCore_R'; 'NAccShell_L'; 'NAccShell_R'; 'NAcc_L'; 'NAcc_R';
            'HPC_aL';  'HPC_aR';   'Hippocampus_L';  'Hippocampus_R';    % HPC/Amygdala
            'Amygdala_L';  'Amygdala_R';
            'Insula_L';  'Insula_R'; 'SNVTA_L'; 'SNVTA_R';  % Exploratory/various
            };
        
    end
    
    %     log.roilist.c13gbattery='c13g battery';
    log.roilist.c13gbattery={'DLPFC_L';'DLPFC_R';'HPC_L_sc';'HPC_L_stc';'HPC_R_satc';'HPC_R_sc';'Parietal';'OFC_L';'Striatum_R';};
    log.roilist.anatHPCamy={'Amyg_L';'Amyg_R';'HPC_aL';'HPC_aR';};
    log.roilist.c14HPC={
        'HPC_aL_c';
        'HPC_aL_stc';
        'HPC_aL_tc';
        'HPC_aR_c';
        'HPC_aR_satc';
        'HPC_aR_stc';
        'HPC_aR_tc'; };
    log.roilist.v26e_pvBest={'Amyg_L';'CingAnt';'HPC_aL';'HPC_aR';'HPC_mL';'vmPFC'};
end
% log.batteryname='c13gbattery';
% log.batteryname='anatHPCamy';
% log.batteryname='c14HPC';
% log.batteryname='v26e_pvBest';
% eval(['log.ROInames=log.roilist.' log.batteryname ';'])
% log.ROInames={'BA10_L';'BA10_L'};
log.ROInames={};  % Empty to use roi names

%% (2) Implement extraction

% Display + Load ROIs
log.roi_files=cellstr(spm_select('List', log.roi_fol, '.nii$')); log.n_rois=length(log.roi_files);
if isfield(log, 'batteryname')==1; disp(['######    Battery name = ' log.batteryname '  ######']); end, disp('------  ROI files       ---------------      ROI short names ------'); disp([log.roi_files   log.ROInames]); input('OK?    ');
d_roimasks=cell(length(log.roi_files),2);   if isempty(log.ROInames)==1 ,  log.ROInames =cellfun(@(x)x(1:strfind(x, '.')-1),  log.roi_files, 'UniformOutput',0); end
for r=1:log.n_rois
    d_roimasks{r,1}=log.ROInames{r};
    d_roimasks{r,2}=spm_read_vols(spm_vol([log.roi_fol filesep log.roi_files{r}]));
end


% Full flexible extraction: Load all images and extract (conditionally per subject)
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per typical plot)
d_roimatrix_cF= nan(log.nlevels,log.nlevels,log.nchoices, length(log.roi_files), log.n_subjs); d_roimatrix_ct=d_roimatrix_cF; l_roimatrix_cF=d_roimatrix_cF; l_roimatrix_ct=d_roimatrix_cF; % log: each subject ok or not?
if log.nchoices==3; log.choicecF={'Accept_';'Reject_';'Explore_'};  log.choicect={'NoBomb_';'Bomb_';'Explore_'};
else log.choicecF={[]};  log.choicect={[]};
end
for s=1:log.n_subjs
    ws.c=clock;  disp(['Subject ' num2str(s) '  (' log.subjects{s} ')            [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  --------------'])
    ws.FLfol=[where.data_brain filesep log.subjects{s} filesep '2 First level s4Ants'  filesep log.FLfol  filesep];
    ws.s=load([ws.FLfol 'SPM.mat']);
    ws.mask=spm_read_vols(spm_vol([ws.FLfol 'mask.img']));
    
    % Load all contrasts for this subject
    ws.con=cell(size(ws.s.SPM.xCon,2),3);
    for c=1:size(ws.s.SPM.xCon,2)
        ws.con{c,1}=ws.s.SPM.xCon(c).name;
        ws.con{c,2}=ws.s.SPM.xCon(c).Vcon.fname;
        ws.con{c,3}=spm_read_vols(spm_vol([ws.FLfol  ws.con{c,2}]));
        %             % Fake contrast inputs;
        %             disp('Contrast extraction turned off! Fake cons!');
        %             ws.con{c,3}=nan(79,95,68);
        %             ws.con{c,3}(1:20,1:20, 1:20)=rand(20,20,20);
        %             ws.con{c,3}(:,:,:)=rand(size(ws.con{c,3}));
    end
    
    % Read betas for all available cF contrast
    disp('        Loading cF')
    for ee=1:log.nlevels
        e=log.nlevels+1-ee;
        for n=1:log.nlevels
            for ch=1:log.nchoices
                wc.conname=['cF_' log.choicecF{ch} log.trialchar{1} num2str(ee) '-' log.trialchar{2} num2str(n)];
                if isempty(find(strcmp(ws.con, wc.conname)))==0
                    wc.con=ws.con{find(strcmp(ws.con(:,1),  wc.conname)), 3};
                    for r=1: log.n_rois
                        % Read all rois for this contrast. If crashes here, check that ROI image and contrast images have the same dimensions!
                        d_roimatrix_cF(e,n,ch,r,s)= nanmean( wc.con(d_roimasks{r,2}~=0));
                        %
                        
                        l_roimatrix_cF(e,n,ch,r,s)=1;
                    end
                else
                    l_roimatrix_cF(e,n,ch,:,s)=0;
                end
                wc=[];
            end
        end
    end
    
    % Read betas for all available ct contrast
    disp('        Loading ct')
    for ee=1:log.nlevels
        e=log.nlevels+1-ee;
        for n=1:log.nlevels
            for ch=1:log.nchoices
                wc.conname=['ct_' log.choicect{ch} log.trialchar{1} num2str(ee) '-' log.trialchar{2} num2str(n)];
                if isempty(find(strcmp(ws.con, wc.conname)))==0
                    wc.con=ws.con{find(strcmp(ws.con(:,1),  wc.conname)), 3};
                    for r=1: log.n_rois % Read all rois for this contrast
                        d_roimatrix_ct(e,n,ch,r,s)=nanmean( wc.con(d_roimasks{r,2}~=0));
                        %
                        l_roimatrix_ct(e,n,ch,r,s)=1;
                    end
                else
                    l_roimatrix_ct(e,n,ch,:,s)=0;
                end
                wc=[];
            end
        end
    end
    
    %
    ws=[];
end
save([log.roi_fol filesep 'Extracted beta matrix (' date ')'], 'd_roimatrix_cF', 'd_roimatrix_ct', 'l_roimatrix_cF', 'l_roimatrix_ct', 'log')
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis script is complete: extracting betas from all fully flexible contrasts (' mfilename ')']), ' ',1);
end
disp('======================================================='); w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
disp(['             First level model:      ' log.FLmodel]); disp('=======================================================')
