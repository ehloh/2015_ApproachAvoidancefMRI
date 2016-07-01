% PPI prep: preprocessing, extract VOIs (VOI_* file)
%       VOI file is in univariate FL folder. All later variables kept in PPI FL folder 
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={};
%     log.FirstLevelThread=[]; log.func_prefix='swubf';
log.FirstLevelThread=' s4Ants';   log.func_prefix='s4ubf'; 
% log.AntsType='_Basic';
log.AntsType=[];

%
request.FLcontrast_Fomnibus=1;
request.FLcontrast_SpecificComparisons=0;
request.extract_firstlevelVOI=0;
log.regressor_famtype=1; % 1=Choice, 2= Parameter Pmods, 3=Trial event, 4=ChoicexTrial event, 5=Cluster choice

% VOI details (seed region from which to read neural signal) ####################
%       This has nothing to do with design for the PPI itself - just the seed region, which can be read from any contrast
request.extract_firstlevelVOI_ImgMask={
%     'BA10_C'; 'BA46_L_C';                       % Frontal-Striatal
%     'Striatum_L_C'; 'Striatum_R_C';
%     'HPC_aL_sTC'; 'HPC_aR_sTC';             % HPC (Task x Choice)
%     'HPC_aL_sC'; 'HPC_aR_saC';             % HPC (ME Choice)
    };
request.extract_firstlevelVOI_PeakSphere={  % Name, Centre, Radius
    'sph_HPC_aL_STC'; 'sph_HPC_aR_sTC';  % c3 TxC
    'sph_BA10_C'; 'sph_Precuneus'; 'sph_Striatum_R_C'; 'sph_SupMidFG'; % Choice Explore
    'sph_HPC_aL_sC'; 'sph_HPC_aR_saC'; % Choice HPC
}; 


% [UNIVAR MODEL DETAILS] ############################
log.firstlevelmodel='m_c3_ChoiceFull_OULPEN';  
% log.firstlevelmodel='m_c7_ChCluster6Full_OULPEN';

for o1=1:1 % General settings and specifications
    

    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==0; where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.exp_folder='C:\Users\eloh\Desktop\2 [Explore]'; where.data_brain=[where.exp_folder filesep '1 Brain data']; addpath(where.where); 
%     where.model=[where.exp_folder filesep '2 Second level results' log.FirstLevelThread filesep log.firstlevelmodel filesep log.secondlevelmodel filesep];
    else where.where='/Users/EleanorL/Dropbox/sandisk/5 Explore fMRI'; where.exp_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain=[where.exp_folder filesep '1 Brain data']; addpath(where.where); 
    end
    
   % Load subjects
   log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_ChCluster4Full_OULPEN';'m_c7_ChCluster6Full_OULPEN';'m_c8_ChCluster4MovFull_OULPEN';'m_c9_ChCluster6MovFull_OULPEN'; 'm_c10_ChCluster6FullRT_OULPEN'};
    if sum(strcmp(log.firstlevelmodel, w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep '4 Set up models' filesep 'i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
    where.subFL=cell(log.n_subjs,1);
    
    for s=1:log.n_subjs  % Locations of 1st level models  
        if strcmp(log.firstlevelmodel(1), 'f')==1;
            where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep 'm_f1 Contrasted   ' log.firstlevelmodel filesep];
            log.univmod_type='flex'; input('Folders not set up yet probably!');
        elseif strcmp(log.firstlevelmodel(1), 't')==1;
            where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep 'm_t1 Contrasted   ' log.firstlevelmodel filesep];
            log.univmod_type='flex'; input('Folders not set up yet probably!');
        else
            where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.firstlevelmodel ' Contrasted' filesep];
            log.univmod_type='par';
        end
    end
    
    % Settings that don't really change
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.firstlevelmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' '); disp('Requested analysis:'); disp(request); disp(' ')
    disp('Univariate model details ------------ ')
    disp(['Type:   ' log.univmod_type])
    disp(['Model: ' log.firstlevelmodel]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% (1) PPI Preprocessing: Change address, implement omnibus & comparison contrasts

% Instructions & Setup for PPI preprocessing
for o1=1:1  
    
    % Contrast/regressor type
    switch log.regressor_famtype
        case 1; log.reg_famname='Choices'; log.reg_famclass=1;
        case 2; log.reg_famname='ParPmods'; log.reg_famclass=2;
        case 3; log.reg_famname='TrialTypeEvent'; log.reg_famclass=1;
        case 4; log.reg_famname='FlexEvent'; log.reg_famclass=1; 
        case 5; log.reg_famname='ClusterChoice'; log.reg_famclass=1;
        otherwise; error('Invalid Contrast type selected. Choice, Parameter pmods, Trialtyle?')
    end
    
    % (A) Change address of functional scans ###########################
    % See execution module
    
    % (B) Implement Omnibus F-tests ######################################################
    %       Instruc: New contrast name, Regressor names, Regressor Types
    switch log.reg_famname
        case 'Choices'
            instruc.Fomni={'Choices'        {'cF_Accept';'cF_Reject';'cF_Explore';'ct_NoBomb';'ct_Bomb';'ct_Explore'}           1;};
            p=load([where.subFL{1} 'SPM.mat']); % Read RL variables included from 1st subject
        case 'ParPmods'
            p=load([where.subFL{1} 'SPM.mat']); 
            instruc.Fomni={'ParPmods'      cellstr(vertcat([repmat('pcF_', length(p.SPM.RLvariables),1) char(p.SPM.RLvariables)],[repmat('pct_', length(p.SPM.RLvariables),1) char(p.SPM.RLvariables)]))       2;};
        case 'TrialTypeEvent' 
            
            % Compile trial-type list of regressors to be included
            w.freg=cell(2*6*6,1); k=1;
            for t=1:2; switch t; case 1, w.pre='cF'; case 2, w.pre='ct'; end; for e=1:6; for n=1:6; w.freg{k}=[w.pre '_t' num2str(e) '-' num2str(n)]; k=k+1; end;end;end
            instruc.Fomni={'TrialTypeEvent'         w.freg          1};
            %
            p=load([where.subFL{1} 'SPM.mat']); % Read RL variables included from 1st subject
        case 'ClusterChoice'
            instruc.Fomni={'ClusterChoice'        {'in_cF_Reject';'in_cF_Explore';'in_ct_Bomb';'in_ct_Explore'}           1;};
%                                     'ClusterChoiceAll'        {'in_cF_Reject';'in_cF_Explore';'in_ct_Bomb';'in_ct_Explore';
%                                                                        'out_cF_Accept';'cF_Reject';'out_ct_NoBomb';'ct_Explore';}           1;};
            p=load([where.subFL{1} 'SPM.mat']); % Read RL variables included from 1st subject
            
        case 'FlexEvent'
            error('Need to super flexibly read the cells to adjust for, because this differs from subject to subject!!!!')
            disp('Need to flexibly read this for every point in the script where it invokes the p.SPM');
        otherwise
            error('What type of reg family for omnibus comparisons?')
    end
    
    % (C) Specific comparisons ######################################################
    %           1=Comparison name, 2=Regressor-conditions names & their weights, 3=t(1) comparison or F(2)
    
    switch log.reg_famname
        case 'Choices'
            instruc.spec_comp={
                    'cF_Rej-Explore'        {'cF_Reject' 1; 'cF_Explore' -1}                                1;      % 
                    'cF_Rej-Others'        {'cF_Reject' 2; 'cF_Accept' -1;'cF_Explore' -1}          1;
                    'cF_Explore-Others'  {'cF_Explore' 2; 'cF_Accept' -1;'cF_Reject' -1}          1;
                    'ct_Rej-Explore'        {'ct_Bomb' 1; 'ct_Explore' -1}                                1;
                    'ct_Rej-Others'        {'ct_Bomb' 2; 'ct_NoBomb' -1;'ct_Explore' -1}          1;
                    'ct_Explore-Others'  {'ct_Explore' 2; 'ct_NoBomb' -1;'ct_Bomb' -1}          1;
                    };
        case 'ClusterChoice'
            instruc.spec_comp={};
%         case 'ParPmods'
%         case 'TrialTypeEvent' 
        otherwise
            error('What type of reg family for specific comparisons?')
    end
    
end

% Execute PPI preprocessing
for o2=1:1 % Change address of fxnals, Implement omnibus F, Add specific comparisons (FL contrasts)
    % (A) Change address of functional scans if necessary
    request.change_fxnscan_address=0;
    if request.change_fxnscan_address
        request.changefxnscan_newaddress='[where.data_brain filesep log.subjects{s} filesep ''2 First level'' log.FirstLevelThread filesep ''Preproc functionals'' filesep]';
        disp('User requested that address of functional scans needs to be changed. New address:'); disp(request.changefxnscan_newaddress); input('Proceed?   ');
        for s=1:log.n_subjs
            try
                ws.FLfol=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.firstlevelmodel  ' Contrasted' filesep];
                ws.ss=load([ws.FLfol 'SPM.mat']);  % Formatted for par models only!
                ws.SPM=ws.ss.SPM; ws.SPM.old=ws.SPM;
%                 ws=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.firstlevelmodel  ' Contrasted' filesep 'SPM.mat']);  % Formatted for par models only!
%                 ws.SPM.old=ws.SPM;
                ws.preprocfolder_address=eval(request.changefxnscan_newaddress); if isdir(ws.preprocfolder_address)==0; error(['Could not find location of preproc fxnals: ' ws.preprocfolder_address]); end
                
                % Identify fxnal scans (in order), re-write to SPM.xY.P & SPM.xY.VY.fname
                ws.fxnfiles=cell(size(ws.SPM.xY.P,1),2);
                for i=1: size(ws.SPM.xY.P,1) % identify
                    if isempty(strfind(ws.SPM.xY.P(i,:), log.func_prefix))==0
                        ws.fxnfiles{i,1}=[ws.preprocfolder_address ws.SPM.xY.P(i,strfind(ws.SPM.xY.P(i,:),log.func_prefix):end)];
                        ws.fxnfiles{i,2}=[ws.preprocfolder_address ws.SPM.xY.VY(i).fname(strfind(ws.SPM.xY.VY(i).fname,log.func_prefix):end)];
                        ws.SPM.xY.VY(i).fname=ws.fxnfiles{i,2};
                    else
                        error('Changing address of func scans - could not identify func scan name')
                    end
                end
                ws.SPM.xY.P=vertcat(char(ws.fxnfiles(:,1)));
                
                % Save
                disp(['Subject ' num2str(s) ' -  ' log.subjects{s} ' :  ' ws.SPM.xY.P(1,:)])
                SPM=ws.SPM;
                save([ws.FLfol 'SPM.mat'],'SPM');
                ws=[];  SPM=[];
            catch
                errorlog{e,1}=['Couldn''t change address of fxn scans for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
            end
        end
    end
    
    % (B) Implement omnibus F-test (all involved conditions)
    if request.FLcontrast_Fomnibus
        disp('############# Implementing Omnibus contrasts ###############')
        disp(['fMRI model:   ' log.firstlevelmodel]);  disp(instruc.Fomni);  disp('Included regressors (adjusting data for effects of interest)'); for i=1:size(instruc.Fomni,1); disp(['F omni #' num2str(i)]); disp(instruc.Fomni{i,2}); end
        input('Instructions OK for this fMRI model?   ');
        for f=1:size(instruc.Fomni,1)
            disp(['F omnibus for ' instruc.Fomni{f,1} ' ######'])
            
            % Identify target regressors and generate their weights
            wf.regnums=f_GetRegressorNums(p.SPM.xX.name',[instruc.Fomni{f,2} num2cell(instruc.Fomni{f,3}*ones(size(instruc.Fomni{f,2})))]);
            
            for s=1:log.n_subjs
                try
                    disp(['Subject ' num2str(s) ' - ' log.subjects{s}])
                    
                    % Compile weights
                    ws.s=load([where.subFL{1} 'SPM.mat']);
                    ws.weights=zeros(size(wf.regnums,1), length(ws.s.SPM.xX.name));
                    for i=1:size(wf.regnums,1)
                        if length(wf.regnums{i,2})==1
                            ws.weights(i, wf.regnums{i,2})=1;
                        end
                    end
                    
                    % Batch in SPM
                    matlabbatch{1}.spm.stats.con.spmmat={[where.subFL{s} 'SPM.mat']};
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name=['Fomni_' instruc.Fomni{f,1}];
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = {ws.weights};
                    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                    matlabbatch{1}.spm.stats.con.delete = 0;
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[]; ws=[];
                catch
                    errorlog{e,1}=['Couldn''t implement Fomni contrasts (' instruc.Fomni{f,1} ') for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
                end
            end
            
            wf=[];
        end
    end

    % (C-1) [Context] Implement specific comparisons at the first-level (subject-specific contrasts)
    if request.FLcontrast_SpecificComparisons
        disp('############# Implementing Specific contrast comparisons###############')
        disp(['Regressor type:   ' log.reg_famname '  --------------'])
        
        for s=1:log.n_subjs
            disp(['Subject ' num2str(s) ' - ' log.subjects{s} ' ----------------' ])
            ws.s=load([where.subFL{s} 'SPM.mat']);
            %
            matlabbatch{1}.spm.stats.con.spmmat={[where.subFL{s} 'SPM.mat']};
            matlabbatch{1}.spm.stats.con.delete = 0;
            
            % Set up all requested comparisons
            for c=1:size(instruc.spec_comp,1)
                wc.regnums=f_GetRegressorNums(ws.s.SPM.xX.name, [instruc.spec_comp{c,2}(:,1) num2cell(repmat(log.reg_famclass,size(instruc.spec_comp{c,2},1),1))]);
                wc.instruc=[instruc.spec_comp{c,2} wc.regnums(:,2)];
                
                % Compile weights
                wc.weights=zeros(1,length(ws.s.SPM.xX.name));
                for i=1:size(wc.instruc,1)
                    if length(wc.instruc{i,3})~=1
                        error('# Identified regressors to be weighted ~=1')
                    end
                    wc.weights(wc.instruc{i,3})=wc.instruc{i,2};
                end
                
                % Set up batch
                wc.batch.name = instruc.spec_comp{c,1};
                wc.batch.convec=wc.weights;
                wc.batch.sessrep = 'none';
                switch instruc.spec_comp{c,3}
                    case 1; matlabbatch{1}.spm.stats.con.consess{c}.tcon=wc.batch;
                    case 2; matlabbatch{1}.spm.stats.con.consess{c}.fcon=wc.batch;
                end
                
                wc=[];
            end
            
            % Execute
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            matlabbatch=[]; ws=[];
        end        
    end
    
end

%% (2) Extract VOIs from all subjects (create VOI_* file)
%   Roi images are held in Experiment folder --> 3 Anatomical --> 3 VOI seeds

if request.extract_firstlevelVOI
    disp('########### Identifying VOIs ###############################')
    
    % Fetch roi images for VOIs 
%     request.VOIformat='.nii';
    log.VOIfiles=cell(size(request.extract_firstlevelVOI_ImgMask,1),1);
    for m=1:size(request.extract_firstlevelVOI_ImgMask,1)
        f=spm_select('List', [where.exp_folder filesep '3 Anatomical' filesep '3 VOI seeds' filesep], [request.extract_firstlevelVOI_ImgMask{m} '.img']);
        if isempty(f)==1; f=spm_select('List', [where.exp_folder filesep '3 Anatomical' filesep '3 VOI seeds' filesep], [request.extract_firstlevelVOI_ImgMask{m} '.nii']); end % Alternate format?
        wm.add= [where.exp_folder filesep '3 Anatomical' filesep '3 VOI seeds' filesep f];
        
        % Last check
        if isempty(f); error(['Could not find requested VOI image file (' request.extract_firstlevelVOI_ImgMask{m} '). Check it''s name/format/existence!']); end
        log.VOIfiles{m}=wm.add;
        disp(['VOI #' num2str(m) ': ' log.VOIfiles{m}]);
    end
    
    % Set up for all subjects
    disp('########### Extracting VOIs ###############################')
    for s=1:log.n_subjs
        %         try
        disp(['Subject ' num2str(s) ' (' log.subjects{s} ') ---------------'])
        ws.s=load([where.subFL{s} 'SPM.mat']);
        ws.contrastlist={ws.s.SPM.xCon(:).name}';
        
        % (A) Extract VOI on the basis of an image mask ####################
        if isempty(request.extract_firstlevelVOI_ImgMask)==0
            for m=1:size(request.extract_firstlevelVOI_ImgMask,1)
                disp(['VOI # ' num2str(m) ' - ' request.extract_firstlevelVOI_ImgMask{m} ])
                
                matlabbatch{1}.spm.util.voi.spmmat = {[where.subFL{s} 'SPM.mat']};
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famname]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
                if isempty(matlabbatch{1}.spm.util.voi.adjust)==1; error('Which contrast to adjust all data for (omnibus F contrast)?'); end
                %
                matlabbatch{1}.spm.util.voi.name = request.extract_firstlevelVOI_ImgMask{m};
                matlabbatch{1}.spm.util.voi.roi{1}.mask.image={[log.VOIfiles{m} ',1']};
                matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image ={[where.subFL{s} 'mask.img,1']};
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
                %
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[]; wr=[];
            end
        end
        
        % (B) Extract VOI from peak + sphere of x radius ####################
        if isempty(request.extract_firstlevelVOI_PeakSphere)==0
            for r=1:size(request.extract_firstlevelVOI_PeakSphere,1)
                disp(['VOI #' num2str(r) '  -  '  request.extract_firstlevelVOI_PeakSphere{r,1} ' (' num2str(request.extract_firstlevelVOI_PeakSphere{r,3}) ' mm sphere at ' num2str(request.extract_firstlevelVOI_PeakSphere{r,2}) ')'])
                %
                matlabbatch{1}.spm.util.voi.spmmat = {[where.subFL{s} 'SPM.mat']};
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famname]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
                if isempty(matlabbatch{1}.spm.util.voi.adjust)==1; error('Which contrast to adjust all data for (omnibus F contrast)?'); end
                
                matlabbatch{1}.spm.util.voi.name = ['sph_' request.extract_firstlevelVOI_PeakSphere{r,1}];
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = request.extract_firstlevelVOI_PeakSphere{r,2};
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = request.extract_firstlevelVOI_PeakSphere{r,3};
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 0.5;
                %
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image ={[where.subFL{s} 'mask.img,1']};
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
                %
                spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                matlabbatch=[]; wr=[];
            end
        end
        
        % [Disused] Conjunction of anatomical ROI + 1st level contrast result
        request.extract_firstlevelVOI_AnatInContrast={};
        if isempty(request.extract_firstlevelVOI_AnatInContrast)==0
            for o1=1:1
                % Sample request code:
                % request.extract_firstlevelVOI_AnatInContrast={
                % 'Hippocampus_L'     {'SimRxCMem_Hit'; 'SimNxCMem_Hit';  'DisRxCMem_Hit';  'DisNxCMem_Hit';};
                % 'SNVTA_R'               {'SimRxCMem_Hit'; 'SimNxCMem_Hit';  'DisRxCMem_Hit';  'DisNxCMem_Hit';};
                % 'Hippocampus_L'     {'cm_SR-SN';'cm_DR-DN'; 'cm_SimxVal'};
                % 'SNVTA_R'               {'cm_SR-SN';'cm_DR-DN'; 'cm_SimxVal'};
                % };
                %
                %     % (B) Extract VOI as conjunction of anat mask &(thresholded) 1st-level contrast result ##################
                %     for r=1:size(request.extract_firstlevelVOI_AnatInContrast,1)
                %
                %
                %
                %     matlabbatch{1}.spm.util.voi.name = [request.extract_firstlevelVOI_AnatInContrast{1} '_in_' request.extract_firstlevelVOI_AnatInContrast{2} ];
                %     matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
                %     matlabbatch{1}.spm.util.voi.adjust = find(strcmp(ws.contrastlist,['Fomni_' log.reg_famname]), 1, 'last'); % Omni-F contrast for THIS type of contrast only
                %
                %     % First-level contrast
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = matlabbatch{1}.spm.util.voi.spmmat;
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1;
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.05;
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
                %     matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
                %
                %     % Anatomical masks
                %     matlabbatch{1}.spm.util.voi.roi{2}.mask.image ={[where.exp_folder filesep '3 Anatomical' filesep '2 A priori ROIs' filesep request.extract_firstlevelVOI_AnatInContrast{1} '.img,1']};
                %     matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.1;
                %
                %     % Display
                %     disp(' '); disp(['VOI extracted:  ' request.extract_firstlevelVOI_AnatInContrast{1} '      in       ' request.extract_firstlevelVOI_AnatInContrast{2} ]);
                %     disp(['Threshold for first-level result: ' num2str(matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh)])
                %     disp(['First level SPM file:   ' matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat{1}])
                %     if s==1; input('Details ok  (Approval for 1st subject only)?    '); end
            end
            error('Requested VOI on the basis of Anatomical ROI in a 1st-level contrast. Code for this is not complete')
        end
        %         catch
        %             errorlog{e,1}=['Couldn''t extract VOI for subject  ' log.subjects{s}]; disp(errorlog{e,1}); e=e+1;
        %         end
    end
end

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
% disp(['GLM model: ' log.FirstLevelThread '      ' log.firstlevelmodel '     ' log.secondlevelmodel]); disp(' ')
disp('Analysis completed:'); disp(request);
disp(' '); disp('Errors:'); disp(errorlog)
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', ['Analysis batchscript is complete (' mfilename ')'], ' ',1);
end
