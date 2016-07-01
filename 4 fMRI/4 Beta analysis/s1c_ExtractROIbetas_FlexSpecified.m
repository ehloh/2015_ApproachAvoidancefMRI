%% s1c_ExtractROIbetas_Specified
% Extract from specified FL addresses + contrasts, print to specified 2nd
%       level model. Mostly for use with flexible models. Uses SPM functions -
%       thus ROIs and subject contrasts MUST be have same voxel properties!
clear all; clc

where.where='I:\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];

% [MODEL DETAILS] ############################
log.firstlevelmodel='f1_TaskChoiceTrialType';
log.secondlevelname='AcrossDesign';

for o1=1:1 % General setup
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, [], [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    % instruc.log.subjects=log.subjects; instruc.log.n_subjs=log.n_subjs;
    
    
    % First level settings
    log.FLfol=['m_' log.firstlevelmodel(1:2) ' Contrasted   ' log.firstlevelmodel];
    p=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.FLfol filesep 'SPM.mat']);
    log.CandidateCons=cellstr(char(p.SPM.xCon.name));
    
    
    % Second level settings (folder only)
    log.secondlevelfol=[where.experiment_folder filesep '2 Second level results' filesep  log.firstlevelmodel filesep log.secondlevelname filesep];    
    if isdir([where.experiment_folder filesep '2 Second level results' filesep log.firstlevelmodel])==0; mkdir([where.data_brain log.firstlevelmodel]); end
    if isdir(log.secondlevelfol)==0; mkdir(log.secondlevelfol); end
    
end

%% Specifications

% ROIs
where.rois='C:\Users\eloh\Desktop\2 [Explore]\2 Second level results\f1_TaskChoiceTrialtype\AcrossDesign\ROIs from c3Choice';
log.rois={ 'HPC_L_strictanatChoice';'HPC_R_strictanatChoice';
                        'HPC_L_strictanatTcf'; 'HPC_R_strictanatTcf';
                        'BA11';'BA46';'Caudate_L_anat';'Caudate_R_anat'; 'Putamen_L_anat';
                        };


% Which contrast to extract from?
for o1=1:1 
    log.contrasts={'cF_e1_n1'
    'cF_e1_n2'
    'cF_e1_n3'
    'cF_e1_n4'
    'cF_e1_n5'
    'cF_e1_n6'
    'cF_e2_n1'
    'cF_e2_n2'
    'cF_e2_n3'
    'cF_e2_n4'
    'cF_e2_n5'
    'cF_e2_n6'
    'cF_e3_n1'
    'cF_e3_n2'
    'cF_e3_n3'
    'cF_e3_n4'
    'cF_e3_n5'
    'cF_e3_n6'
    'cF_e4_n1'
    'cF_e4_n2'
    'cF_e4_n3'
    'cF_e4_n4'
    'cF_e4_n5'
    'cF_e4_n6'
    'cF_e5_n1'
    'cF_e5_n2'
    'cF_e5_n3'
    'cF_e5_n4'
    'cF_e5_n5'
    'cF_e5_n6'
    'cF_e6_n1'
    'cF_e6_n2'
    'cF_e6_n3'
    'cF_e6_n4'
    'cF_e6_n5'
    'cF_e6_n6'
    'ct_e1_n1'
    'ct_e1_n2'
    'ct_e1_n3'
    'ct_e1_n4'
    'ct_e1_n5'
    'ct_e1_n6'
    'ct_e2_n1'
    'ct_e2_n2'
    'ct_e2_n3'
    'ct_e2_n4'
    'ct_e2_n5'
    'ct_e2_n6'
    'ct_e3_n1'
    'ct_e3_n2'
    'ct_e3_n3'
    'ct_e3_n4'
    'ct_e3_n5'
    'ct_e3_n6'
    'ct_e4_n1'
    'ct_e4_n2'
    'ct_e4_n3'
    'ct_e4_n4'
    'ct_e4_n5'
    'ct_e4_n6'
    'ct_e5_n1'
    'ct_e5_n2'
    'ct_e5_n3'
    'ct_e5_n4'
    'ct_e5_n5'
    'ct_e5_n6'
    'ct_e6_n1'
    'ct_e6_n2'
    'ct_e6_n3'
    'ct_e6_n4'
    'ct_e6_n5'
    'ct_e6_n6'};

log.n_cons=length(log.contrasts);
end

%% Extract betas & print to 

% Data variables (v_=volumes,  d_ = data)
v_subcon=cell(log.n_subjs, log.n_cons); % Row = subject, Col = Requested Contrast
d_betas=[[{' '}; log.rois] [cellstr(log.contrasts');  repmat({nan(log.n_subjs,1)}, log.n_rois, log.n_cons)]];  % Row =Roi, Col = Contrast; within cells, row=subject;
v_rois=cell(log.n_rois,1);

% Set up ROI masks (v_rois)
log.n_rois=length(log.rois);
for r=1:log.n_rois
    f=spm_select('List', where.rois, [log.rois{r} '.*.img']);
    if isempty(f);  f=spm_select('List', where.rois, [log.rois{r} '.*.nii']); end
    if isempty(f)==1; error(['Could not find requested ROI mask: ' log.rois{r}]); end
    %
    v_rois{r}=spm_read_vols(spm_vol([where.rois filesep f]));
end

% Extract
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) '  (' log.subjects{s} ')']);
    ws.whereFL=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.FLfol filesep];
    ws.s=load([ws.whereFL 'SPM.mat']);
    
    
    % Load all available & requested contrasts from subject
    ws.subcons=cellstr(char(ws.s.SPM.xCon(:).name));
    for i=1:length(ws.subcons)
        ws.subcons{i,2}=ws.s.SPM.xCon(i).Vcon.fname;
    end
    ws.reqcons=cell(log.n_cons,2);
    for c=1:log.n_cons
        ws.reqcons{c,1}=log.contrasts{c};
        if length(find(strcmp(ws.subcons(:,1),  log.contrasts{c})))~=1; error(['Could not find specific contrast for specific subject: '  log.contrasts{c} '  '  log.subjects{s}]); end
        ws.reqcons{c,2}=ws.subcons{find(strcmp(ws.subcons(:,1),  log.contrasts{c})),2};
        
        % Load
        v_subcon{s,c}=spm_read_vols(spm_vol([ws.whereFL ws.reqcons{c,2}]));
    end
    
    
    % Extract each ROI from each contrast
    for r=1:log.n_rois
        for c=1:log.n_cons
            disp(['Extracting roi #' num2str(r) '  -  contrast #' num2str(c) '  (' log.rois{r} '  '  log.contrasts{c}   ')'])
            d_betas{r+1, c+1}(s)=mean(v_subcon{s,c}(v_rois{r}~=0));
        end
    end
    
    %
    ws=[];
end


%% Save

details.log=log; details.where=where;
save([log.secondlevelfol 'Extracted betas (' date ')'], 'd_betas', 'details')

try % Notify researcher
    f_sendemail('kurzlich', ['Analysis script is complete (' mfilename ')'], ' ',1);
end
