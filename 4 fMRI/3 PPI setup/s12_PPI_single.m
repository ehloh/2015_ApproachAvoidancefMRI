% PPI pipeline:
%   Prep script compiles variables that are not unique to the PPI (seed & psych comparison) -
%           i.e. VOI_*.mat file. Lives in univariate folder.
%   All other files (which are unique to the seed & psych comparison) live in the PPI first-level folder 
%           i.e. PPI_*.mat file
clear all;close all hidden; clc

% Requested analysis ###############
log.specificsubjects={};
%
request.construct_PPIterm=0;
request.PPImodel_FirstLevel=1;
    request.Spec=0; 
    request.Est=1; 
    request.Contrast=0; 
request.PPImodel_SecondLevel=0;
%
log.FirstLevelThread=[];; log.func_prefix='swubf';
%
log.PPI_modelname='ChoiceCompare';   % Choice model PPIs --------------
% log.PPI_modelname='ChoiceCells';
% log.PPI_modelname='ClusterChoiceCompare';    % Cluster-Choice model PPIs --------------

% First level PPI setup 
instruc.VOIname={
%     'HPC_aL_sTC'
%     % -------------------------------------------------------
    'BA10_C'; 'BA46_L_C';                       % Frontal-Striatal
    'Striatum_L_C';   'Striatum_R_C';
%     'HPC_aL_sTC'; 'HPC_aR_sTC';            % HPC (Task x Choice)
%     'HPC_aL_sC'; 'HPC_aR_saC';             % HPC (ME Choice)
%     'Amygdala_saC';
    };


% Second level setup ###############
request.secondlevel.model='onesampleT';
request.secondlevel.VOI='Caudate_R';
% request.secondlevel.PsychVariable='TaskxRejExp';
request.secondlevel.PsychVariable='Rej-Exp';
% request.secondlevel.PsychVariable='cF_Rej-Exp';


% [UNIVARIATE MODEL DETAILS] ############################
log.firstlevelmodel='m_c3_ChoiceFull_OULPEN';                 
% log.firstlevelmodel='m_c7_ChCluster6Full_OULPEN';
for o1=1:1 % Archive: other univariate models 
% log.firstlevelmodel='m_c1_CompeteBasic';                  % ----------- Par models ----------
% log.firstlevelmodel='m_c2_CompeteBasicVExplore';
% log.firstlevelmodel='m_c4_CompeteFull_XUPEN';
% log.firstlevelmodel='m_c5_CompeteFullRT_XUVPEN';
% log.firstlevelmodel='m_c6_Cluster4CompeteFull_XUVPEN';
% log.firstlevelmodel='f1_TaskChoiceTrialtype';     % ----------- Flex models -----------
% log.firstlevelmodel='f2_ExploreOr_AllQualCells';
% log.firstlevelmodel='f3_ExploreOr_FixWind4';
% log.firstlevelmodel='f4_ExploreOr_MovWind';
end

for o1=1:1 % General settings and specifications    
 
    % Paths: Different for the PPIs vs univariate data
    w.w=pwd; 
    if strcmp(w.w(1), '/')==0;  where.where='D:\Dropbox\SANDISK\5 Explore fMRI';  where.exp_folder='C:\Users\eloh\Desktop\2 [Explore]';  where.data_brain=[where.exp_folder filesep '1 Brain data'];
        where.datappi_brain='F:\2 Explore fMRI\1 Brain data PPI';
    else where.where='/Users/EleanorL/Dropbox/sandisk/5 Explore fMRI';  where.exp_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data';  where.data_brain=[where.exp_folder filesep '1 Brain data'];
        where.datappi_brain=[where.exp_folder filesep '1b Brain data PPI'];
    end
    path(pathdef); addpath(where.where); addpath([where.where filesep '6 PPI setup' filesep 'Second level models']);
    
   % Load subjects + selections
   w.modelsneedingsubselect={'m_c6_ChCluster4Full_OULPEN';'m_c7_ChCluster6Full_OULPEN';'m_c8_ChCluster4MovFull_OULPEN';'m_c9_ChCluster6MovFull_OULPEN'; 'm_c10_ChCluster6FullRT_OULPEN'}; 
   log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
   [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
   if sum(strcmp(log.firstlevelmodel, w.modelsneedingsubselect))==1; disp('Executing additional subject selections for this univariate model');
        [w.s w.s1 log.koshertable]=xlsread([where.where filesep '4 Set up models' filesep 'i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.firstlevelmodel);
    end
    where.subFLfol_univar=cell(log.n_subjs,1); where.subFLfol_ppi=where.subFLfol_univar;
    for s=1:log.n_subjs  % Locations of 1st level models  (univariate and ppi)
        where.subFLfol_univar{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level'  log.FirstLevelThread filesep ];
        where.subFLfol_ppi{s}=[where.datappi_brain filesep log.subjects{s} filesep '2 First level'  log.FirstLevelThread filesep ];
        if strcmp(log.firstlevelmodel(1:3), 'm_c')==0; error('Specify FL locations!'); end
%         if strcmp(log.firstlevelmodel(1), 'f')==1;
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_f1 Contrasted   ' log.firstlevelmodel filesep];
%             log.univmod_type='flex';
%         elseif strcmp(log.firstlevelmodel(1), 't')==1;
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep 'm_t1 Contrasted   ' log.firstlevelmodel filesep];
%             log.univmod_type='flex';
%         else
%             where.subFL{s}=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelmodel ' Contrasted' filesep];
%             log.univmod_type='par';
%         end
    end
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.firstlevelmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)]); if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location (brain): ' where.exp_folder]);
    disp(' '); disp('Requested analysis:'); disp(request); disp(' ');
    disp(' '); disp('             -------------------  CHECK HERE  -------------------'); disp(' ')
    disp(['Univariate model: '  log.firstlevelmodel ]);
    if request.construct_PPIterm==1 ||  request.PPImodel_FirstLevel==1
        disp('Constructing first-level models. VOIs:'); disp(instruc.VOIname);
    elseif request.PPImodel_SecondLevel==1
        disp(['VOI seed:              ' request.secondlevel.VOI]);
        disp(['Psych variable:      ' request.secondlevel.PsychVariable]);
        disp(['Second level stat:  ' request.secondlevel.model]);
    end
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
end

%% (1) Set up 

for o1=1:1 % WHICH comparisons do you want?  ###### * * * ######
    
    % PPI Model type & weights (Col 1=Condition name, Col 2=Regressor type ???ithink??, Col 3=Comparison weight
    switch log.PPI_modelname
        case 'ChoiceCompare'
        % (A) [Choice model] Specific comparisons ###############################
             instruc.PPIweights={
%                 'TaskxChoice'    {'cF_Reject'         1           1;             % (A) Task x Choice Interaction
%                                             'cF_Explore'        1          -1;
%                                             'ct_Bomb'           1          -1;
%                                             'ct_Explore'        1           1};
%                  'TaskxRejExp'  {'cF_Reject'       1            1;
%                                          'cF_Explore'      1           -1;
%                                          'ct_Bomb'         1           -1;
%                                          'ct_Explore'       1            1};
%                  'cF_Rej-Exp'     {'cF_Reject'       1           1;                 % (B) (TxC) Within-Task comparisons
%                                           'cF_Explore'     1           -1;};
%                 'cF_Rej-Oth'     {'cF_Reject'       1           2;
%                                           'cF_Accept'      1           -1;
%                                           'cF_Explore'     1           -1;};
                'cF_Exp-Oth'     {'cF_Reject'       1           -1;
                                          'cF_Accept'       1           -1;
                                          'cF_Explore'      1           2;};
%                 'ct_Rej-Exp'     {'ct_Bomb'         1           1;     
%                                          'ct_Explore'       1           -1;};
%                 'ct_Rej-Oth'     {'ct_Bomb'         1           2;
%                                          'ct_NoBomb'     1           -1;
%                                          'ct_Explore'       1           -1;};
                'ct_Exp-Oth'     {'ct_Bomb'         1           -1;
                                          'ct_NoBomb'     1           -1;
                                          'ct_Explore'       1           2;}; 
%                  'Rej_cF-ct'       {'cF_Reject'       1           1;                 % (C) (TxC) Across-Task comparisons
%                                            'ct_Bomb'       1           -1;};
%                 'Exp_cF-ct'        {'cF_Explore'     1           2;
%                                            'ct_Explore'       1           -1;};
%                 ---------------------------------------------------------------------------------------------------   
%                  'Rej-Exp'         {'cF_Reject'       1            1;                  % (D) Combined across tasks
%                                          'cF_Explore'      1           -1;
%                                          'ct_Bomb'         1            1;
%                                          'ct_Explore'       1           -1};
%                 ---------------------------------------------------------------------------------------------------   
                 'Exp-Oth'         {'cF_Accept'       1           -1;                  % (E) Exploration comparisons
                                          'cF_Reject'       1           -1;
                                          'cF_Explore'      1           2;
                                          'ct_NoBomb'     1           -1;
                                          'ct_Bomb'         1           -1;
                                         'ct_Explore'        1           2;};
%                  'Exp-Rej'         {'cF_Reject'       1           -1;
%                                           'cF_Explore'      1           1;
%                                           'ct_Bomb'         1           -1;
%                                          'ct_Explore'        1           1;};

                };
            %
            log.regressor_famtype=1;
            log.reg_famprefix='comp';
        case 'ChoiceCells'
        % (B) [Choice model] Cell means ########################################
             instruc.PPIweights={
%                 'cF_Accept'         {'cF_Accept'        1            1;         % (B) Cell means ################
%                                             'cF_Reject'         1           0;
%                                             'cF_Explore'        1          0;
%                                             'ct_NoBomb'      1           0;
%                                             'ct_Bomb'           1          0;
%                                             'ct_Explore'        1           0};
%                 'cF_Reject'         {'cF_Accept'        1            0;
%                                             'cF_Reject'         1           1;
%                                             'cF_Explore'        1          0;
%                                             'ct_NoBomb'      1           0;
%                                             'ct_Bomb'           1          0;
%                                             'ct_Explore'        1           0};
%                 'cF_Explore'         {'cF_Accept'        1            0;
%                                             'cF_Reject'         1           0;
%                                             'cF_Explore'        1          1;
%                                             'ct_NoBomb'      1           0;
%                                             'ct_Bomb'           1          0;
%                                             'ct_Explore'        1           0};
%                 'ct_NoBomb'        {'cF_Accept'        1            0;
%                                             'cF_Reject'         1           0;
%                                             'cF_Explore'        1          0;
%                                             'ct_NoBomb'      1           1;
%                                             'ct_Bomb'           1          0;
%                                             'ct_Explore'        1           0};
%                 'ct_Bomb'           {'cF_Accept'        1            0;
%                                             'cF_Reject'         1           0;
%                                             'cF_Explore'        1          0;
%                                             'ct_NoBomb'      1           0;
%                                             'ct_Bomb'           1          1;
%                                             'ct_Explore'        1           0};
%                 'ct_Explore'          {'cF_Accept'        1            0;
%                                             'cF_Reject'         1           0;
%                                             'cF_Explore'        1          0;
%                                             'ct_NoBomb'      1           0;
%                                             'ct_Bomb'           1          0;
%                                             'ct_Explore'        1           1};
                };                         
            %                            
            log.regressor_famtype=1;
            log.reg_famprefix='cell';
        case 'ClusterChoiceCompare'
            % (A) [Choice model] Specific comparisons ###############################
             instruc.PPIweights={
                'TaskxInChoice'    {'in_cF_Reject'         1           1;             % (A) Task x Choice Interaction
                                            'in_cF_Explore'        1          -1;
                                            'in_ct_Bomb'           1          -1;
                                            'in_ct_Explore'        1           1};
 %                 'cF_Rej-Exp'     {'cF_Reject'       1           1;                 % (B) Within-Task comparisons
%                                           'cF_Explore'     1           -1;};
%                 'cF_Rej-Oth'     {'cF_Reject'       1           2;
%                                           'cF_Accept'      1           -1;
%                                           'cF_Explore'     1           -1;};
%                 'cF_Exp-Oth'     {'cF_Reject'       1           -1;
%                                           'cF_Accept'       1           -1;
%                                           'cF_Explore'      1           2;};
%                 'ct_Rej-Exp'     {'ct_Bomb'         1           1;     
%                                          'ct_Explore'       1           -1;};
%                 'ct_Rej-Oth'     {'ct_Bomb'         1           2;
%                                          'ct_NoBomb'     1           -1;
%                                          'ct_Explore'       1           -1;};
%                 'ct_Exp-Oth'     {'ct_Bomb'         1           -1;
%                                           'ct_NoBomb'     1           -1;
%                                           'ct_Explore'       1           2;};}; 
                            };
            %
            log.regressor_famtype=5;
            log.reg_famprefix='comp';
        otherwise
            error('Specified PPI-model not set up yet')
    end
        
% #################################################################  
% #################################################################
% #################################################################
    
    
    % Contrast/regressor type
    switch log.regressor_famtype
        case 1;
            log.reg_famname='Choices'; 
            log.reg_famclass=1;
%         case 2; 
%             log.reg_famname='ParPmods'; 
%             log.reg_famclass=2;
%             log.reg_famprefix='ppar';
%         case 3; 
%             log.reg_famname='TrialTypeEvent'; 
%             log.reg_famclass=1;
%             log.reg_famprefix='tt';
%         case 4; 
%             log.reg_famname='FlexEvent'; 
%             log.reg_famclass=1; 
%             log.reg_famprefix='f';
        case 5;
            log.reg_famname='ClusterChoice'; 
            log.reg_famclass=1;
        otherwise; error('Invalid regressor type selected. Choice, Parameter pmods, Trialtyle?')
    end
    
    % Set up PPI names
    log.ppi_names=cell(size(instruc.VOIname,1)*size(instruc.PPIweights,1),3);
    for v=1:size(instruc.VOIname,1)
        for w=1:size(instruc.PPIweights,1);
            k=(v-1)*size(instruc.PPIweights,1)+w;
            log.ppi_names{k,1}=[log.reg_famprefix 'Fam_' instruc.VOIname{v} '_psy_' instruc.PPIweights{w,1}];
            log.ppi_names{k,2}=instruc.VOIname{v};
            log.ppi_names{k,3}=instruc.PPIweights{w,1};
        end
    end
    
    % Folders - univariate model branch at 1st level
    for s=1:log.n_subjs
        if isdir([where.subFLfol_univar{s} log.firstlevelmodel ' PPIs'])==0; mkdir([where.subFLfol_univar{s} log.firstlevelmodel ' PPIs']); end
    end
    
    % 2nd level setup
    if request.PPImodel_SecondLevel==1
        log.PPImodelname=[log.reg_famprefix 'Fam_' request.secondlevel.VOI '_psy_' request.secondlevel.PsychVariable];
        input(['Change name?  Currently =' log.PPImodelname]);
        where.secondlevelres=[where.exp_folder filesep '2 Second level results' filesep log.firstlevelmodel filesep log.PPImodelname filesep];
        if isdir(where.secondlevelres)==0; mkdir(where.secondlevelres); end
    end
  
end

% Displace interface ####################################################
if request.construct_PPIterm==1 || request.PPImodel_FirstLevel==1
    disp('---------------------------- REQUESTED PPI MODEL -------------------------------------'); disp(' ')
    disp(['PPI model name: ' log.PPI_modelname])
    disp(['Regressor type:  ' log.reg_famname]); disp(' ')
    disp('Contrast weights for all requested psychological comparisons:'); disp(' '); for i=1:size(instruc.PPIweights,1); disp([num2str(i) ')  ' instruc.PPIweights{i,1}]); disp(instruc.PPIweights{i,2});    end; disp(' ')
    disp('PPI names: '); disp(log.ppi_names(:,1));disp(' '); input('OK?    ');
end
% ###############################################################
    
%% (2) Construct PPI term (in PPI FL folder, PPI_*.mat file)

if request.construct_PPIterm
    
    for o1=1:1 % Delete existing VOI terms? 
        log.delete_Old_PPIterms=0;
        if log.delete_Old_PPIterms
            input('User requested deleting ALL existing old PPI terms from 1st level folders. Proceed?')
            disp('Deleting old PPI terms --------')
            for s=1:log.n_subjs 
                f=spm_select('List',[where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs'],'^PPI.*.mat$');
                for i=1:size(f,1)
                    delete([where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs' filesep f(i,1:strfind(f(i,:), '.mat')+3)]);
                end
            end      
        end
    end
    
    disp('########### Generating PPI terms for each VOI x Psychological comparison ###############')
  
    % Load subject univariate SPM mats (for formatting instructions within subject loop)
    log.subSPM=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        log.subSPM{s}=load([where.subFLfol_univar{s} log.firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
    end
  
    % Generate PPI terms according to instructions
    for v=1:size(instruc.VOIname,1) % For each VOI
        disp(['VOI #' num2str(v) ':  ' instruc.VOIname{v}  '  #########################'])
        for w=1:size(instruc.PPIweights,1); % For each Psychological comparison
            disp(['Psych Comparison #' num2str(w) ':  ' instruc.PPIweights{w,1}  '  --------'])
            k=(v-1)*size(instruc.PPIweights,1)+w;
            
            for s=1:log.n_subjs
                disp(['Subject no. ' num2str(s) ' (' log.subjects{s} ') ----------------'])
                ws.wheremod=[where.subFLfol_univar{s} log.firstlevelmodel ' Contrasted' filesep] ;
                ws.wherethisPPI=[where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs' filesep log.ppi_names{k,1} filesep];
                
                % Format instructions for PPI weights: Name of conditions to Condition #
                ws.conds=cell(size(log.subSPM{s}.SPM.Sess.U,2),2);
                for i=1:size(log.subSPM{s}.SPM.Sess.U,2)
                    ws.conds{i,1}=log.subSPM{s}.SPM.Sess.U(i).name{1};
                    if size(log.subSPM{s}.SPM.Sess.U(i).name,2)==2; ws.conds{i,2}=log.subSPM{s}.SPM.Sess.U(i).name{2}; end
                end
                for i=1:size(instruc.PPIweights,1)  % substitute condition # for condition names
                    for j=1:size(instruc.PPIweights{i,2},1)
                        instruc.PPIweights{i,3}(j,2)=instruc.PPIweights{i,2}{j,2};
                        instruc.PPIweights{i,3}(j,3)=instruc.PPIweights{i,2}{j,3};
                        if sum(strcmp(ws.conds(:,log.reg_famclass), instruc.PPIweights{i,2}{j,1}))==1
                            instruc.PPIweights{i,3}(j,1)=find(strcmp(ws.conds(:,log.reg_famclass), instruc.PPIweights{i,2}{j,1}));
                        else error(['Couldn''t identify condition in setting up psych comparison (generating PPI term): ' instruc.PPIweights{i,2}{j,1}]);
                        end
                    end
                end
                
                % Batch
%                 try
                    matlabbatch{1}.spm.stats.ppi.spmmat = {[ws.wheremod 'SPM.mat']};
                    matlabbatch{1}.spm.stats.ppi.disp =s==1;
                    %
                    matlabbatch{1}.spm.stats.ppi.name =log.ppi_names{k,1};
                    matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {[ws.wheremod 'VOI_' instruc.VOIname{v} '_1.mat']};
                    matlabbatch{1}.spm.stats.ppi.type.ppi.u =instruc.PPIweights{w,3};
                    %
                    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                    matlabbatch=[];
                    
                    % Append PPI term to add details
                    PPI=[]; load([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat']);
                    PPI.details.SeedVOI= instruc.VOIname{v};
                    PPI.details.PsychCond_name=instruc.PPIweights{w,1};
                    PPI.details.PsychCond_weights=instruc.PPIweights{w,3};
                    PPI.details.PsychCond_AssumedRegs=ws.conds; %
                    PPI.details.univarmodel=log.firstlevelmodel;
                    PPI.details.name=log.ppi_names(w,:);
                    PPI.details.regdetails={log.reg_famname log.reg_famclass log.reg_famprefix};                    
                    save([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat'], 'PPI');
                   
                    % Move PPI term over to its own folder
                    if isdir(ws.wherethisPPI)==0; mkdir(ws.wherethisPPI); end
                    movefile([ws.wheremod 'PPI_' log.ppi_names{k,1} '.mat'], [ws.wherethisPPI 'PPI_' log.ppi_names{k,1} '.mat']);
                    
%                 catch
%                     errorlog{e,1}=['Couldn''t generate PPI term subject  ' log.subjects{s} '  - ' instruc.VOIname{v} '  - ' instruc.PPIweights{w,1}]; disp(errorlog{e,1}); e=e+1;
%                 end
                ws=[];
            end
        end
    end
end

%% (3) Set up first-level for PPI models (Specify, Estimate, Contrast)

if request.PPImodel_FirstLevel
    
    if request.Spec
        disp('Specifying PPI models ################################')
        for o1=1:1 % Specification - Settings for first-level models
            scan=load([where.where filesep '3 Scripts - Preprocessing' filesep 'i_scanningdetails.mat']);
            settings.firstlevelmodelspec.timing.units = 'secs';
            settings.firstlevelmodelspec.timing.RT = scan.TRms/1000;
            settings.firstlevelmodelspec.timing.fmri_t = 16;
            settings.firstlevelmodelspec.timing.fmri_t0 = 1;
            settings.firstlevelmodelspec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
            settings.firstlevelmodelspec.sess.regress = struct('name', {}, 'val', {});
            %         settings.firstlevelmodelspec.sess.multi_reg = {''};
            settings.firstlevelmodelspec.sess.hpf = 128;
            settings.firstlevelmodelspec.fact = struct('name', {}, 'levels', {});
            settings.firstlevelmodelspec.bases.hrf.derivs =[1 1];
            settings.firstlevelmodelspec.volt = 1;
            settings.firstlevelmodelspec.global = 'None';
            settings.firstlevelmodelspec.mask = {''};
            settings.firstlevelmodelspec.cvi = 'AR(1)';
            %         settings.firstlevelmodelspec.cond.sess.pmod.poly = 1; % Polynomials?
        end
        for v=1: size(instruc.VOIname,1) % For each VOI
            for w=1: size(instruc.PPIweights,1); % For each Psychological comparison
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Specifying models for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
                    try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') ------------'])
                        ws.whereproc=[where.subFLfol_univar{s} 'Preproc functionals' filesep];
                        ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep];
                        ws.wherePPImodel=[where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs' filesep log.ppi_names{k,1} filesep];
                                                
                        % Typical model specification
                        matlabbatch{1}.spm.stats.fmri_spec= settings.firstlevelmodelspec;
                        matlabbatch{1}.spm.stats.fmri_spec.dir = {ws.wherePPImodel};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[ws.where log.subjects{s} '_reg_physiomovement.txt']};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
                        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
                        f=spm_select('List', ws.whereproc, ['^' log.func_prefix '.*img$']); % Functional files
                        matlabbatch{1}.spm.stats.fmri_spec.sess.scans =cellstr([repmat(ws.whereproc , size(f,1),1) f]);
                        
                        % PPI model specification
                        ws.p=load([ws.wherePPImodel 'PPI_' log.ppi_names{k,1} '.mat']);
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI';    % Interaction
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = [ws.p.PPI.ppi];
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = 'P';      % Psychological
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = [ws.p.PPI.P];
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Y';      % Extracted signal
                        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = [ws.p.PPI.Y];
                        
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[];ws=[];
                    catch
                        errorlog{e}=['ERROR: Could not specify model for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
                    end
                end
            end
        end
    end
    
    if request.Est
        disp('Estimating PPI models ################################')
        for v=1: size(instruc.VOIname,1)
            for w=1: size(instruc.PPIweights,1);
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Estimating model for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
%                     try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
                        ws.wheremodel=[where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs' filesep log.ppi_names{k,1} filesep];
                        %
                        f   = spm_select('List', ws.wheremodel, 'SPM.mat');
                        if isempty(f); error('Could not find Specified SPM.mat file to Estimate model. Has model been specified?'); end
                        matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([ws.wheremodel f]);
                        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[];wb=[];
%                     catch
%                         errorlog{e}=['ERROR: Could not estimate model for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
%                     end
                end
            end
        end
    end
    
    if request.Contrast
        disp('Running Contrasts PPI models ################################')
        for v=1:size(instruc.VOIname,1)
            for w=1:size(instruc.PPIweights,1);
                k=(v-1)*size(instruc.PPIweights,1)+w;
                disp(['Running contrasts for PPI #' num2str(k) ': ' log.ppi_names{k,1} '------------------------------------------------------'])
                
                for s=1:log.n_subjs
                    try
                        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
                        ws.wherePPImodel=[where.subFLfol_ppi{s} log.firstlevelmodel ' PPIs' filesep log.ppi_names{k,1} filesep];
                        %
                        matlabbatch{1}.spm.stats.con.spmmat = {[ws.wherePPImodel 'SPM.mat']};
                        matlabbatch{1}.spm.stats.con.delete =1;
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'PPI Pos';
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1];
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'PPI Neg';
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1];
                        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                        %
                        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
                        matlabbatch=[]; ws=[];
                    catch
                        errorlog{e}=['ERROR: Could not run contrasts for subject   ' log.subjects{s} ' - ' log.ppi_names{k,1}]; disp(errorlog{e}); e=e+1;
                    end
                end
            end
        end
    end
    
end            
            

%% (4) Second level (Specific single VOI & Psych variable only)

if request.PPImodel_SecondLevel
    disp(['Running 2nd level model: ' request.secondlevel.model ' ################################'])
    eval(['[ batch] = '  request.secondlevel.model '(where, log, ''PPI Pos'');'])
end

%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
% disp(['GLM model: ' log.FirstLevelThread '      ' log.firstlevelmodel '     ' log.secondlevelmodel]); disp(' ')
disp('Analysis completed:'); disp(request);
disp('Errors:'); disp(errorlog)
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', ['Analysis batchscript is complete (' mfilename ')'], ' ',1);
end
