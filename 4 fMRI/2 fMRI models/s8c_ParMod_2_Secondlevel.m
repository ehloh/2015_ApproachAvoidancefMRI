% Second level analysis - Parametized (Compete/Orthog/Trialtype) models only
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={};
% log.specificsubjects={'p02_YY';'p04_MP';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p35_SM';'p36_FR';'p38_MK';'p41_AC'};

% Which model? ########################
log.FirstLevelThread=' s4Ants';  
log.AntsType='_Basic';
% log.AntsType='_LM4';
% log.AntsType='_Basic';
% log.FirstLevelThread=' s6'; log.AntsType=[];
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
% log.onsetsmodel='m_c1_Choice_ENU';
% log.onsetsmodel='m_v10c_RejectOrvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
% log.onsetsmodel='m_t1_1_Trialtype';
% log.onsetsmodel='m_v23g_RejectOrvChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_v6g_ChoicevChosenAndposneg2_bpji08bpji11';
% log.onsetsmodel='m_v25g_RejectOrvGamble_bpji08bpji11';
% log.onsetsmodel='m_c14_Choice';
% log.onsetsmodel='m_v26e_pvBestChoice_b01b01';
% log.onsetsmodel='m_v24g_predChoice_bpji08bpji11';
% log.onsetsmodel='m_v4g_ChoicevChosenAnd_bpji08bpji11';
% log.onsetsmodel='m_c15g_NextChoice_ULPEN';
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
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
% log.onsetsmodel='m_v34g_ChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_v35g_ChoiceChoiceXvGamble_bpji08bpji11';
% log.onsetsmodel='m_c18_pChoice';
% log.onsetsmodel='m_c19_pChoiceFull_ULPEN';
% log.onsetsmodel='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11';
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v38g_EVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v40g_vBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v42g_pLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v39g_ChoiceEVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v41g_ChoicevBestvWorst_bpji08bpji11'; 
end
% log.onsetsmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v45g_ChoicevBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v43g_ChoicepLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11'; 
% log.onsetsmodel='m_c6_ChCluster4Full_ULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 
% log.onsetsmodel='m_c21_RejExpFull_ULPEN';  log.secondlevelmodel='TaskVar_2x2';
% log.onsetsmodel='m_c22_ChoiceFull_HLPEN'; 
log.onsetsmodel='m_c23_ChoiceFullRT_ULPEN'; 
% log.onsetsmodel='m_c24_ChCluster6RT_ULPEN'; 


log.secondlevelmodel='choice_2x3';                 % SECOND LEVEL MODEL###############
% log.secondlevelmodel='choiceRT_2x3';
% log.secondlevelmodel='choice_cluster2x2';
% log.secondlevelmodel='choice_2x2';
% log.secondlevelmodel='choice_cluster2x2ar';
% log.secondlevelmodel='TaskVar_2x3';
% log.secondlevelmodel='TaskVar_2x2';
% log.secondlevelmodel='RejOr_vBU';
% log.secondlevelmodel='choiceRT_cluster2x2';
% log.secondlevelmodel='Identity_1samplettest';
% log.secondlevelmodel='Identity_pairedttest';
% log.secondlevelmodel='RL_2xN_fakefactorial';
% log.secondlevelmodel='ExVEx_2x2';     % Comparing RL variables
% log.secondlevelmodel='MiscSecondLevel';

% % Trial-type models - Linear FX ####################
% log.secondlevelmodel='RL_linearFX';
% log.onsetsmodel='m_t1_Trialtype';
% log.RLvariables={'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 'Conflict'; 'EV';};

for o1=1:1 % General settings and specifications
        
    % Append Ants Type?
    log.onsetsonlymodel=log.onsetsmodel; log.onsetsmodel=[log.onsetsmodel log.AntsType];
    
    % Add paths
    where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; where.experiment_folder='G:\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
    addpath(where.where);  addpath(genpath([where.where filesep '4 Set up models']));
    
    % Matching?
    if strcmp(log.secondlevelmodel, 'RL_linearFX')==1 
        log.onsetsmodel='m_t1_Trialtype'; disp('Onsets model changed to m_t1_Trialtype, to match 2nd level model requested')
    end
    if isempty(strfind( log.onsetsmodel, 'Trialtype'))==0 % Trialtype models stored elsewhere!
        where.data_brain='I:\1 Explore fMRI';
    end
    
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    
    % Apply further subject selection for some models
    log.exploreinfo_oksubs={'p02_YY';'p04_MP';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p35_SM';'p36_FR';'p38_MK';'p41_AC'};
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_c10_';'m_v7_';'m_c16';'m_c17';'m_c24'};
    if sum(strcmp(log.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread('i_Subjectdataok_SpecificModels.xlsx'); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.onsetsonlymodel);
    end
    
    % Log
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.onsetsmodel ' - ' log.secondlevelmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1; 
    log.specificRLvariables={};
%     log.specificRLvariables={'EnvThreat';'pLoss';};
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             First level model:  ' log.onsetsmodel])
    disp(['             Second level model:  ' log.secondlevelmodel]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% Which Variables are included in this model? (sample 1st subject)

for o1=1:1  % Settings for 2nd level models
    choices={'cF_Accept' [1 1];'cF_Reject' [1 2];'cF_Explore'  [1 3];'ct_NoBomb'  [2 1];'ct_Bomb'  [2 2];'ct_Explore'  [2 3]};
    if isempty(strfind(log.secondlevelmodel, 'choice')); choices=[]; end
    
    
    
%     choices={'cF_Accept' [1 1];'cF_Reject' [1 2];'cF_Explore'  [1 3];'ct_NoBomb'  [2 1];'ct_Bomb'  [2 2];'ct_Explore'  [2 3]};
    
    % Which Variables to actually include?
    if strcmp(log.onsetsonlymodel(1:3), 'm_v')==1 % Value models
        % they're not RLvariables, but specifying it this way allows for use of existing scripts
        if strcmp(log.onsetsonlymodel, 'm_v2_vBestAnd_bpmi16bpmi11')==1;  RLvariables={'vBest';'vSecondBest'};
        elseif sum(strcmp(log.onsetsonlymodel(1:5), {'m_v1c'; 'm_v1e'}))==1; RLvariables={'vExplore'}; % RLvariables={'vAccept';'vExplore'};
        elseif sum(strcmp(log.onsetsonlymodel(1:4),{'m_v3'; 'm_v4'; 'm_v7'}))==1 && isnumeric(str2num(log.onsetsonlymodel(5)))==0
            
%         elseif sum(strcmp(log.onsetsonlymodel,{'m_v3_vChosenAnd_bpmi16bpmi11'; 'm_v4_ChoicevChosenAnd_bpmi16bpmi11'; 'm_v7_ChClustervChosenAnd_bpmi16bpmi11'}))==1
            RLvariables={'vChosen';'vBestUnchosen'};
%         elseif sum(strcmp(log.onsetsonlymodel(1:4) ,{'m_v4'}))==1 % All variants of v4 !!
%             RLvariables={'vChosen';'vBestUnchosen'};
        elseif sum(strcmp(log.onsetsonlymodel,{'m_v5_ChoiceXvChosenAnd_bpmi16bpmi11';}))==1
            if strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; RLvariables={ 'Acc_vChosen';'Acc_vBestUnchosen'; 'Rej_vBestUnchosen';'Exp_vChosen';'Exp_vBestUnchosen'};
            elseif strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1; RLvariables={'Rej_vChosen';}; disp('Alter 2nd level parmod Identity_1samplettest: change prefix to ct_ only!!');
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:4),{'m_v6';}))==1
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1;RLvariables={'vBestUnchosen_neg'}; input('Identity_1samplettest script must turn off ct task. Continue if this has been done    ');
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; RLvariables={'vBestUnchosen_pos'}; %                 RLvariables={'vChosen';'vBestUnchosen_pos'};
            else RLvariables={};
            end
        elseif strcmp(log.onsetsonlymodel(1:4),'m_v8')==1
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1;
                RLvariables={'subEVneg'};  input('Identity_1samplettest script must turn off ct task. Continue if this has been done    ');
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; RLvariables={'subEVpos'};
            else RLvariables={};
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:4),{'m_v9';}))==1            
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1; RLvariables={'vBestUnchosen_neg'}; input('Identity_1samplettest script must turn off ct task. Continue if this has been done    ');
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; RLvariables={'vBestUnchosen_pos'};%                 RLvariables={'vChosen';'vBestUnchosen_pos'}; 
            else RLvariables={};
            end
        elseif strcmp(log.onsetsonlymodel(1:5),'m_v10')==1
            if strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; 
                RLvariables={'Rej_vBestUnchosen'; 'NonRej_vBestUnchosen'};    
            else input('2nd level not specified yet. Might not be valid!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v11';'m_v12'}))==1
            log.Outcometype=log.onsetsonlymodel(strfind(log.onsetsonlymodel, 'Outcome')+7:18+strfind(log.onsetsonlymodel(20:end), '_b'));
            RLvariables={'CueValue';['Outcome' log.Outcometype]};
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v13';'m_v14';'m_v22'}))==1
            log.Outcometype=log.onsetsonlymodel(14+strfind(log.onsetsonlymodel(15:end), 'Outcome')+7:18+strfind(log.onsetsonlymodel(20:end), '_b'));
            RLvariables={ 'CueValue';['Outcome' log.Outcometype]};
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v15'}))==1
            if isempty(strfind(log.onsetsmodel(1:30),'InfoVal'))==0;  RLvariables={ 'vExplore';'Info'; 'InfoPE'}; % Cue-Value models where both are
            else RLvariables={ 'vExplore';'Info'};
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v16';}))==1
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1; input('Is 2nd level (Identity_1samplettest) set to ct only?');  RLvariables={'GamblevReject'; 'RejectOutcomePE'};
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                ; RLvariables={'GamblevExplore'; 'ExploreInfotoOutcomePE'; 'ExploreInfoPE'}; %             All:    RLvariables={ 'GamblevAccept'; 'GamblevExplore'; 'AcceptOutcomePE'; 'ExploreInfotoOutcomePE'; 'ExploreInfoPE'};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v17';}))==1
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1
                input('Is 2nd level (Identity_1samplettest) set to ct only?');  RLvariables={'GamblevReject'; 'RejectOutcomePE'};
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1; 
                RLvariables={'ExploreInfoVal'};
%                     'GamblevExplore';
                    ; % All:  RLvariables={ 'GamblevAccept'; 'GamblevExplore'; 'AcceptOutcomePE'; 'ExploreInfotoOutcomePE'; 'ExploreInfoPE'};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v18';}))==1
            if strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1   % Dont need this really
                input('Is 2nd level (Identity_1samplettest) set to ct only?');  RLvariables={'GamblevReject'; 'RejectOutcomePE'};
            elseif strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'OutcomevExplore'; 'ExploreInfoVal';  'ExplorePE_GamtoOutcome'; 'ExplorePE_GamtoInfo'; 'ExplorePE_InfotoOutcome'};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v19';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'vExplore';'ExploreInfoPE';'ExploreInfotoOutcomePE';};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v20';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'predExplore';};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v21';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'ExploreGambletoOutcomePE'; 'ExploreInfoPE'; 'ExploreInfotoOutcomePE'};
            else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v1g';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'vExplore'; 'vAccept'};
            else error('Wrong second level!')
            end
            
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v26';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'pvBest'};
            else error('Wrong second level!')
            end
            
%         elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v24';}))==1
%             if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
%                 RLvariables={'predChoice'};
%             else error('Wrong second level!')
%             end

        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v32';'m_v33'}))==1
            RLvariables={'vGamble'};
%             if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                
%             else error('Wrong second level!')
%             end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v38';'m_v39';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'EVGain'};
            elseif  strcmp(log.secondlevelmodel, 'Identity_1samplettest')==1                
                RLvariables={'EVLoss'};
            else  RLvariables={};
%             else error('Wrong second level!')
            end
        elseif sum(strcmp(log.onsetsonlymodel(1:5),{'m_v42';'m_v43';}))==1
            if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1                
                RLvariables={'pLossNTok'};
            else error('Wrong second level!')
            end
    % ##### [End of value models]  #########################################################        
        else
            input('Assuming no RL variables in this model. Proceed?  ');
            RLvariables=[];
%             if  sum(strcmp(log.onsetsonlymodel(1:5),{'m_v23';'m_v25';'m_v24';'m_v27';'m_v28';'m_v29';'m_v34';'m_v35';'m_v37';'m_v3g';'m_v36';'m_v40';'m_v41';'m_v44';'m_v45';'m_v46'}))==1
%                 RLvariables=[];
%             else error('Which RL variables for this value model?');
%             end
        end
    % ##############################################################
    
    elseif strcmp(log.onsetsonlymodel(1:3), 'm_t')==0 && isempty(log.specificRLvariables)==1 % Trialtype models do not exclude RL variables
        w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']);
        RLvariables=w.spm.SPM.RLvariables;
        
        % Full OULPEN models: Paired ttests cannot be applied to EntropyNTok, EV, NTokens
        if isempty(strfind(log.onsetsonlymodel, 'OULPEN'))==0 && strcmp(log.secondlevelmodel, 'Identity_pairedttest')==1
            RLvariables(strcmp(RLvariables, 'EntropyNTok '))=[];
            RLvariables(strcmp(RLvariables, 'EV'))=[];
            RLvariables(strcmp(RLvariables, 'NTokens'))=[];
        end
        
    elseif strcmp(log.onsetsonlymodel(1:3), 'm_t')==0
        RLvariables=log.specificRLvariables;
    else
        RLvariables=[];
    end
    
    disp('#################################################')
    disp('The following conditions/RL variables are  included in this model (FL con names):')
    disp(' '); disp(RLvariables); disp(' ');
    disp('#################################################')
end

% RLvariables={'vChoMvBU'};
RLvariables={}; 

%% (2) Second-level model specification + Estimation (In function)

disp('STEP 1: Specify 2nd-level model  ####################')

% Results thread for this first level model
where.resultsfolder=[where.experiment_folder filesep '2 Second level results' log.FirstLevelThread filesep log.onsetsmodel]; %  filesep log.secondlevelmodel];
if isdir(where.resultsfolder)==0; mkdir(where.resultsfolder); end

input('USER: Continue to specify + estim model?     ');

% Specify requested model
eval([' [ matlabbatch contrastfiles] = ' log.secondlevelmodel '(log, where,  log.subjects, log.onsetsmodel, log.secondlevelmodel, RLvariables,choices);'])


%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed (s8_Secondlevel)'); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of log.specificsubjects only:'); disp(log.subjects); end
disp(' '); disp(['Data location (brain): ' where.data_brain]); disp(' ')
disp(['First level model:  ' log.onsetsmodel])
disp(['Second level model:  ' log.secondlevelmodel])
disp(' '); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s8_Secondlevel)'), ' ',1);
end
cd(where.resultsfolder)