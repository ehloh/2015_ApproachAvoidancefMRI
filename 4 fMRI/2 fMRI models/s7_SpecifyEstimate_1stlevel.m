% First level analysis - Specify & Estimate model
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={};
% log.specificsubjects={'p06_KB'};
log.specificsubjects={'p01_GV'};
% log.specificsubjects={'p02_YY';'p04_MP';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p35_SM';'p36_FR';'p38_MK';'p41_AC';};
    
% log.specificsubjects=  {
%     'p01_GV';'p02_YY';
%     'p06_KB';
%     'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';
%     'p34_TB';
%     'p36_FR';
%     'p38_MK';'p41_AC'};


request.specify=0;
request.RemoveEventOnsetsBeforeEstim=1; % Value, Competing & Orthog onsets families only
request.estimate=0;

% Which model? --------------
request.includederivatives=1;
log.FirstLevelThread=' s4Ants'; log.prefix='s4ubf';
% log.FirstLevelThread=' s6'; log.prefix='swubf';  % Prefix must indicate smoothing size!!

for o1=1:1 % Comleted models 
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
% log.onsetsmodel='m_o1_OrthogBasic';
% log.onsetsmodel='m_f1_ChoicexTrialtype';
% log.onsetsmodel='m_f2_ChunkChoicexTrialtype';
% log.onsetsmodel='m_f3_Event';
% log.onsetsmodel='m_t1_Trialtype';
% log.onsetsmodel='m_t2_TrialtypeNc';
% log.onsetsmodel='m_t3_ChunkTrialtype';
% log.onsetsmodel='m_t4_ChunkTrialtypeNc';
% log.onsetsmodel='m_c3_ChoiceFull_OULPEN'; 
% log.onsetsmodel='m_c4_Choice_OUPEN'; 
% log.onsetsmodel='m_c6_ChCluster4Full_OULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; 
% log.onsetsmodel='m_v1c_vChoice_bpm16bpmi11';
% log.onsetsmodel='m_v2_vBestAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v3_vChosenAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v3c_vChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v4_ChoicevChosenAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v5c_ChoiceXvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v6e_ChoicevChosenAndposneg2_b01b01';
% log.onsetsmodel='m_v7_ChClustervChosenAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v8c_ChoicesubEVposneg_bpm16bpmi11';
% log.onsetsmodel='m_v9c_vChosenAndposneg2_bpm16bpmi11';
% log.onsetsmodel='m_v10c_RejectOrvChosenAndposneg2_bpm16bpmi11';
% log.onsetsmodel='m_t1_Trialtype';
% log.onsetsmodel='m_t2_TrialtypeNc';
% log.onsetsmodel='m_v11c_vGambleOutcomeMagnitude_bpm16bpmi11'; 
% log.onsetsmodel='m_v12c_vModalchoiceOutcomeMagnitude_bpm16bpmi11'; 
% log.onsetsmodel='m_v13c_OutcomevGambleOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v14c_OutcomevModalchoiceOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoPE_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoPEsign_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoVal_bpm16bpmi11';
% log.onsetsmodel='m_v16c_vChoicePE_bpm16bpmi11';
% log.onsetsmodel='m_v17c_vChoiceOutcomeMag_bpm16bpmi11';
% log.onsetsmodel='m_v18c_vChoiceAtOutcomeMag_bpm16bpmi11';
% log.onsetsmodel='m_v19c_VExploreInfoPEOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v20c_pExplore_bpm16bpmi11';
% log.onsetsmodel='m_v21c_ExploreGamInfoOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v11e_vGambleOutcomePE_b01b01';  % OutcomePE
% log.onsetsmodel='m_v21e_ExploreGamInfoOutcomePE_b01b01';
% log.onsetsmodel='m_v22c_OutcomevChosenOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v22e_OutcomevChosenOutcomeMagnitude_b01b01';
% log.onsetsmodel='m_v22f_OutcomevChosenOutcomeMagnitude_b02b01';
% log.onsetsmodel='m_c1_Choice_ENU';
% log.onsetsmodel='m_v23c_RejectOrvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
% log.onsetsmodel='m_v23g_RejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v6g_ChoicevChosenAndposneg2_bpji08bpji11';
% log.onsetsmodel='m_v24g_predChoice_bpji08bpji11';
% log.onsetsmodel='m_v25g_RejectOrvGamble_bpji08bpji11';
% log.onsetsmodel='m_c14_Choice';
% log.onsetsmodel='m_v26e_pvBestChoice_b01b01';
% log.onsetsmodel='m_c15g_NextChoice_ULPEN';
% log.onsetsmodel='m_v4g_ChoicevChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v27g_ChoiceRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v28g_ChoicePredChoice_bpji08bpji11';
% log.onsetsmodel='m_v29g_ChoiceRejectOrvGamble_bpji08bpji11';
% log.onsetsmodel='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v31g_ChoicecFRejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_c16_ChCluster4';
% log.onsetsmodel='m_c17_ChCluster6';
% log.onsetsmodel='m_v32g_vGamble_bpji08bpji11';
% log.onsetsmodel='m_v33g_ChoicevGamble_bpji08bpji11';
% log.onsetsmodel='m_v8g_vGamblePosNeg_bpji08bpji11';
% log.onsetsmodel='m_v34g_ChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_v35g_ChoiceChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_c18_pChoice';
% log.onsetsmodel='m_c19_pChoiceFull_ULPEN';
% log.onsetsmodel='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11';
% log.onsetsmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v38g_EVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v39g_ChoiceEVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v40g_vBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v41g_ChoicevBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v42g_pLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v43g_ChoicepLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v45g_ChoicevBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11'; 
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_c6_ChCluster4Full_ULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 
% log.onsetsmodel='m_c21_RejExpFull_ULPEN'; 
% log.onsetsmodel='m_c22_ChoiceFull_HLPEN'; 
end
% log.onsetsmodel='m_c23_ChoiceFullRT_ULPEN'; 
% log.onsetsmodel='m_c24_ChCluster6RT_ULPEN'; 
log.onsetsmodel='m_c13g_ChoiceFull_ULPEN';

% log.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd;  if strcmp(w.w(1), '/')==0; 
        where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI';   
        
        where.data_brain='G:\2 [Explore]\1 Brain data';
        where.data_brainfunc=where.data_brain;
        if  isempty(strfind(log.onsetsmodel, 'Trialtype')) ==0; where.data_brain='I:\1 Explore fMRI'; end
    else; where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI';  where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data'; where.data_brainfunc =where.data_brain;
    end
    addpath(where.where);  addpath([where.where filesep '3 Scripts - Preprocessing']);
    addpath(genpath([where.where filesep '4 Set up models']));
    
   % Load subjects
   log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    % log.incluster4subs= {'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};
% log.incluster6subs={'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'}; 

    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_v7_';'m_c16';'m_c17'; 'm_c24'};
    if sum(strcmp(log.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread(['i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.onsetsmodel);
    end
    
    % Settings that don't really change
    scan=load('i_scanningdetails.mat');
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.onsetsmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' '); disp('Requested analysis:'); disp(request); disp(' ')
    disp(['Model: ' log.onsetsmodel]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% STEP 1:  Specify model

for o1=1:1 % Determine which RL variables are included (Sample 1st subject) 
    disp('Determining which RL variables are included (sampling 1st subject)  ##################')
    allvar={'EnvThreat';'NTokens'; 'pLoss'; 'ChoiceH'; 'Entropy'; 'EntropyNTok';'VExplore'; 'Conflict'; 'EV'; 'OutcomeMean'; 'OutcomeVariance'};
    w.v=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.subjects{1} '_onsets_' log.onsetsmodel '.mat']);
  
    % Names of all conditions/pmods in this onsets model
    for i=1:size(w.v.names,2) % condition-regressors
        w.var{i,1}=char(w.v.names(i));
    end
    j=size(w.v.names,2)+1;
    for i=1:size(w.v.pmod,2) % pmod
        for k=1:size(w.v.pmod(i).name,2)
            if isempty(w.v.pmod(i).name)==0
                w.var{j,1}=char(w.v.pmod(i).name{k}); j=j+1;
            end
        end
    end
    % Delete the info pmods 
    w.var(find(    cell2mat(cellfun(@(x)~isempty(strfind(x,'_ExploredBombx')), w.var, 'UniformOutput',0))))= []; 
    w.var(find(    cell2mat(cellfun(@(x)~isempty(strfind(x,'_ExploredNobombx')), w.var, 'UniformOutput',0))))= []; 
     
    % For each variable, is it included in the model?
    for i=1:length(allvar)
        %             w.varthere=cell2mat(strfind(w.var, allvar{i}));
        if isempty(cell2mat(strfind(w.var, allvar{i})))==0
            allvar{i,2}=1;
        else
            allvar{i,2}=0;
        end
    end
    
    % Compile list of RL variables in this model
    j=1; RLvariables=[];
    for i=1:size(allvar,1)
        if allvar{i,2}==1
            RLvariables{j,1}=allvar{i,1}; j=j+1;
        end
    end
    
    if strcmp(log.onsetsmodel(1:3), 'm_v')==1;
        RLvariables=[];
        
    elseif sum(strcmp(log.onsetsmodel(1:5), {'m_c22'; 'm_c23'; 'm_c24'}))==1; % ULPEN. cant be bothered to code this up properly to exclude EntropyNTok
        RLvariables={'EnvThreat'; 'NTokens';'pLoss'; 'Entropy';'EV'; };
%     elseif sum(strcmp(log.onsetsmodel(1:5), {'m_c18';'m_c19'}))==1
%         RLvariables={'pReject' 2; 'pExplore' 2};
    elseif sum(strcmp(log.onsetsmodel(1:5), {'m_c22'}))==1;  
        RLvariables={'EnvThreat'; 'NTokens';'pLoss'; 'ChoiceH'; 'EV'; };
%     elseif sum(strcmp(log.onsetsmodel(1:5), {'m_c18';'m_c19'}))==1
%         RLvariables={'pReject' 2; 'pExplore' 2};
    end
    disp('RL variables included in this model: ')
    disp(RLvariables)
    disp('###############################################')
end

if request.specify==1; % Specify: Format trial onsets & other regressors
    disp('STEP 1: Specifying model ################################')
    for o1=1:1 % Settings 
        settings.firstlevelmodelspec.timing.units = 'secs';
        settings.firstlevelmodelspec.timing.RT = scan.TRms/1000*scan.nSlicesPerVol;
        settings.firstlevelmodelspec.timing.fmri_t = 16;
        settings.firstlevelmodelspec.timing.fmri_t0 = 1;
        settings.firstlevelmodelspec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        settings.firstlevelmodelspec.sess.regress = struct('name', {}, 'val', {});
        settings.firstlevelmodelspec.sess.multi_reg = {''};
        settings.firstlevelmodelspec.sess.hpf = 128;
        settings.firstlevelmodelspec.fact = struct('name', {}, 'levels', {});
        if request.includederivatives==1; settings.firstlevelmodelspec.bases.hrf.derivs =[1 1];
        else settings.firstlevelmodelspec.bases.hrf.derivs =[0 0]; disp('Derivatives are NOT included');
        end
        settings.firstlevelmodelspec.volt = 1;
        settings.firstlevelmodelspec.global = 'None';
        settings.firstlevelmodelspec.mask = {''};
        settings.firstlevelmodelspec.cvi = 'AR(1)';
%         settings.firstlevelmodelspec.cond.sess.pmod.poly = 1; % Polynomials?
    end
    for s=1: log.n_subjs % Specify model for each subject
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        matlabbatch{1}.spm.stats.fmri_spec= settings.firstlevelmodelspec;
%         try
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep];
            wb.wherefuncs=[where.data_brainfunc filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep];
            wb.wheremodel=[wb.where log.onsetsmodel ' Estimated' filesep];
            if isdir(wb.wheremodel)==0; mkdir(wb.wheremodel); end
            % SPECIFY MODEL -------------------  
            matlabbatch{1}.spm.stats.fmri_spec.dir = {wb.wheremodel}; % Specification files
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi ={[wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat']};
            if exist([wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat'], 'file')==0; input(['ERROR: Could not find onsets for the specified model   --   ' log.onsetsmodel ]); end
            onsetsmatrix{s}=load([wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[wb.wherefuncs log.subjects{s} '_reg_physiomovement.txt']};
            % Choose functional files
            f=spm_select('List', [wb.wherefuncs 'Preproc functionals'], ['^' log.prefix '.*img$']);  % This does assume functional scans are in order of runs!
            wb.func=[repmat([wb.wherefuncs 'Preproc functionals' filesep], size(f,1),1) f];
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans =cellstr(wb.func);    
%             save('batch','matlabbatch')
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
           
            
            % Rename SPM file , append details (RL variables included)
            eval('java.io.File([wb.wheremodel ''SPM.mat'']).renameTo(java.io.File([wb.wheremodel ''SPM_'' log.onsetsmodel  ''.mat'']));')
            wb.spm=load([wb.wheremodel filesep 'SPM_' log.onsetsmodel '.mat']);
            SPM=wb.spm.SPM; SPM.RLvariables=RLvariables;
            save([wb.wheremodel filesep 'SPM_' log.onsetsmodel '.mat'],'SPM');
            matlabbatch=[];wb=[]; SPM=[];
%         catch
%             errorlog{e}=['ERROR: Could not specify model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
%         end
        ws=[];
    end
end

%% STEP 2: Remove unecessary regressors

if request.RemoveEventOnsetsBeforeEstim==1;
    disp('STEP 2: Removing regressors before model estimation ##############')
    
    % Which regressors to remove?
    switch log.onsetsmodel(1:3) 
        case 'm_o' % Orthogonalized models: remove trial onsets for each task
            reg2remove={'cF_onset' 1; 'ct_onset' 1};
        case 'm_c' % Competing models: remove trial onsets for each task x RL variable
            if sum(strcmp(log.onsetsmodel(1:5), {'m_c18'}))==1
                reg2remove={'cF_onset' 1; 'ct_onset' 1};
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_c19'; 'm_c20'}))==1
                reg2remove1=horzcat(RLvariables, num2cell(ones(size(RLvariables,1),1)));    
                reg2remove2={'cF_onset' 1; 'ct_onset' 1};
                reg2remove=[reg2remove1; reg2remove2];
            else reg2remove=horzcat(RLvariables, num2cell(ones(size(RLvariables,1),1)));    
            end
        case 'm_v'
            if strcmp(log.onsetsmodel(1:4), 'm_v5')==1
                reg2remove={'cF_Acc' 1;  'cF_Rej' 1; 'cF_Exp' 1;  'ct_NoB' 1;  'ct_Bom' 1; 'ct_Exp' 1; };
%             elseif strcmp(log.onsetsmodel(1:5), 'm_v10')==1
%                 reg2remove={'cF_onset' 1;'cF_Rej_onset' 1;'cF_Rej_onset' 1;'cF_NonRej_onset' 1;'cF_NonRej_onset' 1;'ct_onset' 1;'ct_Rej_onset' 1;'ct_Rej_onset' 1;'ct_NonRej_onset' 1;'ct_NonRej_onset' 1};
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_v13';'m_v14';'m_v22'}))==1
                reg2remove={'cF_Outcome_onset' 1; 'ct_Outcome_onset' 1};
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_v10';}))==1
                reg2remove={'cF_onset' 1; 'ct_onset' 1};
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_v15';}))==1
                if isempty(strfind(log.onsetsmodel, 'VExploreInfoVal'))==0
                    reg2remove={ 'cF_ExploreInfo_onset' 1; 'ct_ExploreInfo_onset' 1};
                else reg2remove={ 'cF_ExploreInfo_onset' 1; 'ct_ExploreInfo_onset' 1;'cF_onset' 1; 'ct_onset' 1; };
                end
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_v16';}))==1
                reg2remove={'cF_onset' 1; 'ct_onset' 1;'cF_OutcomeOnset' 1; 'ct_OutcomeOnset' 1; 'cF_ExploreInfo_onset' 1; 'ct_ExploreInfo_onset' 1};   
            elseif sum(strcmp(log.onsetsmodel(1:5), {'m_v17'; 'm_v18';'m_v19';'m_v20';'m_v21'}))==1
                reg2remove={'cF_onset' 1; 'ct_onset' 1;'cF_OutcomeOnset' 1; 'ct_OutcomeOnset' 1; 'cF_ExploreInfo_onset' 1; 'ct_ExploreInfo_onset' 1};
            else reg2remove={'cF_onset' 1; 'ct_onset' 1};
            end
        otherwise
            error('User requested that regressors to remove, but has not yet specified which')
    end
    
    % Remove regressors
    for s=1:log.n_subjs         
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        
        % Load details from the model to be altered
        wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep];
        ws=load([wb.where 'SPM_' log.onsetsmodel '.mat']);
        ws.regnames=ws.SPM.xX.name';
        
        % Compile list of regressors 
        [ ws.regdetails ] = f_GetRegressorNums(ws.regnames, reg2remove );
        
        j=1;
        for i=1:size(ws.regdetails,1) % Expand list to include derivatives
            for r=1:length(ws.regdetails{i,3})
                ws.removelist{j,1}=ws.regdetails{i,3}{r}; j=j+1;
                ws.removelist{j,1}=[ws.regdetails{i,3}{r}(1:length(ws.regdetails{i,3}{r})-2) '2)']; j=j+1; % 2nd derivative
                ws.removelist{j,1}=[ws.regdetails{i,3}{r}(1:length(ws.regdetails{i,3}{r})-2) '3)']; j=j+1; % 3rd derivative
            end
        end
        for i=1:length(ws.removelist)  % Check: Any requested removals that don't exist?
            if isempty(strfind(ws.regnames,ws.removelist{i}))
                error('Error: Regressor scheduled for removal does not exist')
            end
        end
        
        % Remove regressors
        ws.newSPM = f_RemoveRegressor(ws.SPM,ws.removelist);
        
        % Save
        SPM=ws.newSPM;
        SPM.originalSPM=ws.SPM;
        SPM.note='SPM file modified before model specification.';
        save([wb.where 'SPM_' log.onsetsmodel '.mat'], 'SPM');
        ws=[]; SPM=[];
    end
    
end

%% STEP 3: Estimate model

if request.estimate==1
    disp('STEP 3: Estimate model ################################')
    for s=1:log.n_subjs
%         try
            disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated'];
            %
            f   = spm_select('List', wb.where, ['SPM_' log.onsetsmodel '.mat']);
            matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([wb.where filesep f]);
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            matlabbatch=[];wb=[];
%         catch
%             errorlog{e}=['ERROR: Could not estimate model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
%             wb=[];
%         end
    end
end

%% END

disp('======================================================='); w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completeC:\Users\e.loh\'); disp(request)
disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' '); 
disp(log); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (s7_Firstlevel_1SpecifyEstimate  -  '  log.onsetsmodel  ')']), ' ',1);
end
