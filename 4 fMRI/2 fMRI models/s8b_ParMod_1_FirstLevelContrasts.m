% First level contrasts - Parameterized (Compete/Orthog) & Trial type onsets models
% Standard contrasts only (read off excel)
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={};
% log.specificsubjects={'p01_GV'};
% log.specificsubjects={'p02_YY';'p04_MP';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p35_SM';'p36_FR';'p38_MK';'p41_AC';};

request.DeleteExistingFLContrasts=1;

% Which spatial? --------------
% log.FirstLevelThread=[]; 
log.FirstLevelThread=' s4Ants';
% request.AntsType_ifHalfDone='_Basic';   % If some contrasts have already been ants transformed, use this (otherwise empty)
request.AntsType_ifHalfDone=[]; 
for o=1:1  % Old models 
% log.onsetsmodel='m_c1_Choice_ENU';                        % Which model? (Compete/Orthog/Trialtype only) --------------
% log.onsetsmodel='m_c2_Choice_ENV'; 
% log.onsetsmodel='m_c3_ChoiceFull_OULPEN'; 
% log.onsetsmodel='m_c4_Choice_OUPEN'; 
% log.onsetsmodel='m_c5_ChoiceRTFull_OULPEN'; 
% log.onsetsmodel='m_c6_ChCluster4Full_OULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; 
% log.onsetsmodel='m_c10_ChCluster6FullRT_OULPEN'; 
% log.onsetsmodel='m_c11_Choice_VOUPEN'; 
% log.onsetsmodel='m_c12_Choice_VOULPEN'; 
% log.onsetsmodel='m_t1_Trialtype';           % For trialtype models, FL contrast folder has new name. Only 1 model (_1) available. 
% log.onsetsmodel='m_t2_TrialtypeNc';
% log.onsetsmodel='m_t3_ChunkTrialtype';
% log.onsetsmodel='m_t4_ChunkTrialtypeNc';
% log.onsetsmodel='m_v1c_vChoice_bpm16bpmi11';
% log.onsetsmodel='m_v2_vBestAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v3_vChosenAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v3c_vChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v4d_ChoicevChosenAnd_bpm16bpm11';
% log.onsetsmodel='m_v5c_ChoiceXvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_v6_ChoicevChosenAndposneg2_bpmi16bpmi11';
% log.onsetsmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11';
% log.onsetsmodel='m_v6e_ChoicevChosenAndposneg2_b01b01';
% log.onsetsmodel='m_v7_ChClustervChosenAnd_bpmi16bpmi11';
% log.onsetsmodel='m_v8c_ChoiceGambleposneg_bpm16bpmi11';
% log.onsetsmodel='m_v9c_vChosenAndposneg2_bpm16bpmi11';
% log.onsetsmodel='m_v10c_RejectOrvChosenAndposneg2_bpm16bpmi11';
% log.onsetsmodel='m_v11c_vGambleOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v11c_vGambleOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v12c_vModalchoiceOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v12c_vModalchoiceOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v13c_OutcomevGambleOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v14c_OutcomevModalchoiceOutcomeMagnitude_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoPE_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoPEsign_bpm16bpmi11';
% log.onsetsmodel='m_v15c_VExploreInfoVal_bpm16bpmi11';
% log.onsetsmodel='m_v16c_vChoicePE_bpm16bpmi11';
% log.onsetsmodel='m_v17c_vChoiceOutcomeMag_bpm16bpmi11';16bpmi11';
% log.onsetsmodel='m_v18c_vChoiceAtOutcomeMag_bpm16bpmi11';
% log.onsetsmodel='m_v19c_VExploreInfoPEOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v20c_pExplore_bpm16bpmi11';
% log.onsetsmodel='m_v21c_ExploreGamInfoOutcomePE_bpm16bpmi11';
% log.onsetsmodel='m_v11e_vGambleOutcomePE_b01b01';  % OutcomePE
% log.onsetsmodel='m_v1e_vChoice_b01b01';
% log.onsetsmodel='m_v22e_OutcomevChosenOutcomeMagnitude_b01b01';
% log.onsetsmodel='m_v22f_OutcomevChosenOutcomeMagnitude_b02b01';
% log.onsetsmodel='m_c1_Choice_ENU';
% log.onsetsmodel='m_v23c_RejectOrvChosenAnd_bpm16bpmi11';
% log.onsetsmodel='m_t2_TrialtypeNc';
% log.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
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
% log.onsetsmodel='m_c13_ChoiceFull_ULPEN';
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
% log.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% log.onsetsmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v38g_EVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v40g_vBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v42g_pLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_v39g_ChoiceEVGainLoss_bpji08bpji11';
% log.onsetsmodel='m_v41g_ChoicevBestvWorst_bpji08bpji11'; 
% log.onsetsmodel='m_v43g_ChoicepLossNTok_bpji08bpji11'; 
% log.onsetsmodel='m_c20g_ChoicePredChoice_bpji08bpji11';
% log.onsetsmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v45g_ChoicevBUposnegvMargChoDiff_bpji08bpji11';
% log.onsetsmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11'; 
% log.onsetsmodel='m_c21_RejExpFull_ULPEN'; 
% log.onsetsmodel='m_c22_ChoiceFull_HLPEN'; 
end
% log.onsetsmodel='m_c23_ChoiceFullRT_ULPEN'; 
log.onsetsmodel='m_c24_ChCluster6RT_ULPEN'; 


% log.onsetsmodel='m_c6_ChCluster4Full_ULPEN'; 
% log.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==0; where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; where.experiment_folder='G:\2 [Explore]'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
        if  isempty(strfind(log.onsetsmodel, 'Trialtype')) ==0; where.data_brain='I:\1 Explore fMRI'; end
    else where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI';    where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain=[where.experiment_folder filesep '1 Brain data'];
    end
    addpath(where.where); 
    addpath(genpath([where.where filesep '4 Set up models']));
    
    % Load subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
%     
% log.incluster4subs= {'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};
% log.incluster6subs={'p01_GV';'p02_YY';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p17_SJ';'p18_MS';'p21_ES';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p36_FR';'p38_MK';'p41_AC'};

    % Apply further subject selection for some models 
    log.exploreinfo_oksubs={'p02_YY';'p04_MP';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p35_SM';'p36_FR';'p38_MK';'p41_AC'};
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_c9_';'m_c10'; 'm_v7_'; 'm_c16';'m_c17'; 'm_c24'};
    if sum(strcmp(log.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 log.koshertable]=xlsread(['i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, log.koshertable,log.onsetsmodel);
    end
    
    % Settings that don't really change
    log.prefix='swubf';% scan=load([where.where filesep 'i_scanningdetails.mat']);
    if strcmp(log.onsetsmodel(1:3), 'm_v')==1; request.contrasttable='i_ParValMod_FirstlevelContrasts'; else request.contrasttable='i_ParMod_FirstlevelContrasts';end
    diary([where.data_brain filesep 'SPM logs'  filesep  'Log ' mfilename ' - ' log.onsetsmodel ' (' date ')' ])
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' ');disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp(['Model: ' log.onsetsmodel])
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% Set up Contrasts

for o1=1:1 % Load default specifications for requested model 
    
    % Which contrast table to use?
    if sum(strcmp(log.onsetsmodel,{'m_c5_ChoiceRTFull_OULPEN';}))==1 % Competing + RT model
        contrasts.contrasttable='CompeteRT';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_c6_';  'm_c7_';  'm_c8_'; 'm_c16';'m_c17';'m_c24'}))==1 % Competing + ClusterChoice model
        contrasts.contrasttable='CompeteCluster';
    elseif sum(strcmp(log.onsetsmodel,{'m_c10_ChCluster6FullRT_OULPEN'}))==1 % Competing + ClusterChoiceRT model
        contrasts.contrasttable='CompeteClusterRT';
    elseif strcmp(log.onsetsmodel(1:3), 'm_c')==1  % Competing models
            contrasts.contrasttable='Compete';
    elseif strcmp(log.onsetsmodel(1:3), 'm_o')==1 % Orthogonalized models
            contrasts.contrasttable='Orthog';
            
    % For flexible models, (f or t), this script is only good for the titular-contrasted model. For all others, see Flexmod scripts
    elseif sum(strcmp(log.onsetsmodel,{'m_t1_Trialtype'; 'm_t2_TrialtypeNc'}))==1 % TrialType models
            contrasts.contrasttable='Trialtype';
    elseif sum(strcmp(log.onsetsmodel,{'m_t3_ChunkTrialtype'; 'm_t4_ChunkTrialtypeNc'}))==1 % ChunkTrialType models
            contrasts.contrasttable='ChunkTrialtype';
%     elseif sum(strcmp(log.onsetsmodel,{'m_f1_ChoicexTrialtype';}))==1 % ChoicexTrialType models   % These don't belong here!
%         contrasts.contrasttable='ChoicexTrialtype';
%     elseif sum(strcmp(log.onsetsmodel,{'m_f2_ChunkChoicexTrialtype';}))==1 % ChunkChoicexTrialType models
%         contrasts.contrasttable='ChunkChoicexTrialtype';

    % Value models
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v1';}))==1  &&  isempty(str2num(log.onsetsmodel(5))) ==1 % vChoice
        contrasts.contrasttable='vChoice'; 
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v2_';  'm_v32';'m_v33'}))==1 % vBestAnd
        contrasts.contrasttable='vBestAnd';
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v3'; 'm_v4'}))==1  & isempty(str2num(log.onsetsmodel(5)))  % vChosenAnd
%     elseif sum(strcmp(log.onsetsmodel,{'m_v3_vChosenAnd_bpmi16bpmi11'; 'm_v4_ChoicevChosenAnd_bpmi16bpmi11';    'm_v4b_ChoicevChosenAnd_bpmi16bpmi11';   'm_v4c_ChoicevChosenAnd_bpm16bpmi11';}))==1 % vChosenAnd
        contrasts.contrasttable='vChosenAnd';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v36';'m_v37'; 'm_v46'}))==1 % 
        contrasts.contrasttable='vChosenAnd';
%     elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v4'}))==1 % vChosenAnd for all version of v4
%         contrasts.contrasttable='vChosenAnd';
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v5';}))==1 % ChoiceXvChosen
        contrasts.contrasttable='ChoiceXvChosenAnd';
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v6';}))==1 % vChosenAnd, vBestUnchosen split by pos/neg
        contrasts.contrasttable='vChosenAnd_posneg';
    elseif sum(strcmp(log.onsetsmodel,{'m_v7_ChClustervChosenAnd_bpmi16bpmi11'}))==1 % vChosenAnd + Choice cluster
        contrasts.contrasttable='vChosenAnd_ChCluster'; 
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v8';}))==1 % vGamble, pos and neg
        contrasts.contrasttable='vGamble_posneg';
    elseif sum(strcmp(log.onsetsmodel(1:4),{'m_v9';}))==1 % vGamble, pos and neg
        contrasts.contrasttable='vChosenAnd_posneg'; disp('Disable choice contrasts in excel table!')
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v10'; 'm_v23'; 'm_v27'; 'm_v30'; 'm_v31'}))==1 % vChosen, pos and neg
        contrasts.contrasttable='Rejor_vChosenAnd'; disp('Rejor_vChosenAnd FL contrasts chosen. SEE excel file to make sure correct vBUs are on!'); 
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v25';'m_v29';}))==1 % vGamble, RejOr
        contrasts.contrasttable='Rejor_vGamble';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v11';'m_v12';'m_v13';'m_v14';'m_v22'}))==1 % vGamble
        contrasts.contrasttable='vCueOutcome';  % Also use this for vModalchoice
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v15';}))==1 % Values related to exploration
        contrasts.contrasttable='vExploreVal';  % Also use this for vModalchoice
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v16';}))==1 % Full model for examining values and PEs
        contrasts.contrasttable='vChoiceExploreOutcomePE';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v17';'m_v18'}))==1 % Full model for examining values and canonical PEs
        contrasts.contrasttable='vChoiceExploreOutcomeMagnitude';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v19';'m_v21'}))==1 
        contrasts.contrasttable='vExplorePutativePEs';
        if strcmp(log.onsetsmodel(1:5),'m_v19')==1; input('Go to contrast-specifying excel sheet, change the 1st 2 columns to vExplore not ExploreGambletoOutcomePE!! Won''t work otherwise'); end
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v20'; 'm_v24';'m_v26';'m_v28'}))==1 
        contrasts.contrasttable='predChoice';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v34';'m_v35'}))==1 
        contrasts.contrasttable='ChoiceXvGamble';
    elseif sum(strcmp(log.onsetsmodel(1:5),{'m_v38';'m_v39';'m_v40'; 'm_v41';'m_v42';'m_v43';'m_v44';'m_v45';}))==1 
        contrasts.contrasttable='vMisc';
    else error('Error in Contrasts setup: Could not find requested onsets model. Which contrast table (i.e. excel sheet) to use?')
    end
    
    % Details for requested contrast table (excel file)
    switch contrasts.contrasttable
        case 'Compete'                    % Standard Choice + RL models ----------
            col.num_conds= 6+2*10 +4+4;
        case 'Orthog'
            col.num_conds=22;
        case 'CompeteRT'
            col.num_conds=28;
        case 'CompeteCluster'
            col.num_conds=30;
        case 'CompeteClusterRT'
            col.num_conds=40;
        case 'Trialtype'                    % Flex/Trialtype models ----------
            col.num_conds=78;
        case 'ChunkTrialtype'
            col.num_conds=24;
        case 'vChoice';                    % Value models ----------
            col.num_conds=6+6;
        case 'vBestAnd';
            col.num_conds=4+6+2;   disp('PredChoice #s of conds ok? If you added more, change settings');
        case 'vChosenAnd';
            col.num_conds=4+6+8;    disp('PredChoice #s of conds ok? If you added more, change settings');
        case 'vChosenAnd_ChCluster';
            col.num_conds=4+6*2;
        case 'vChosenAnd_posneg';
            col.num_conds=6+6;
        case 'ChoiceXvChosenAnd';
            col.num_conds=6+6*2;
        case 'vGamble_posneg';
            col.num_conds=4+6;
        case 'Rejor_vChosenAnd_posneg';
            col.num_conds=4+6+6;
        case 'vCueOutcome';
            col.num_conds=6+4;
        case 'vExploreVal';
            col.num_conds=6+4;
        case 'vChoiceExploreOutcomePE';
            col.num_conds=6+6+2;
        case 'vChoiceExploreOutcomeMagnitude';
            col.num_conds=6+6+2;
        case 'vExplorePutativePEs';
            col.num_conds=6;
        case 'predChoice';
            col.num_conds=6+6+2;   disp('PredChoice #s of conds ok? If you added more, change settings');
        case 'Rejor_vChosenAnd';
            col.num_conds=6+6+1;
        case 'Rejor_vGamble';
            col.num_conds=6+4;
        case 'ChoiceXvGamble';
            col.num_conds=6+6;
        case 'vMisc';
            col.num_conds=6+3+4+2+3;   disp('PredChoice #s of conds ok? If you added more, change settings');
        otherwise
            error('Error in Contrasts setup: Requested contrast table (i.e. excel sheet) not found')
    end
    col.contrastname=2;
    col.weightstart=3;
    row.conditionName=1;
    row.conditionType=2; % Condition or Pmod
    col.weightend=col.weightstart+col.num_conds-1;
    col.contrastnum=col.weightend+1;
    col.requested=col.contrastnum+1;
    
end

% input('Continue to compile requested Contrasts?    ')

for o1=1:1 % Compile requested contrasts 
    
    % Load default details regarding available contrasts
    [w.a1, w.a, w.req]=xlsread(request.contrasttable, contrasts.contrasttable);
    for i=col.contrastname+1:col.contrastnum-1 % Read details
        contrasts.condnames{i-col.contrastname}=w.req{row.conditionName,i};
        contrasts.condtypes{i-col.contrastname}=w.req{row.conditionType,i};
    end
    
    % Construct details for requested contrasts
    i=1;
    for r=row.conditionType+1:size(w.req,1);
        if isnan(w.req{r,col.requested})==0
            contrasts.contrastnames{i,1}=w.req{r,col.contrastname};
            contrasts.requested(i,1)=w.req{r,col.requested};
            
            % Load condition weights
            for c=col.contrastname+1:col.contrastnum-1 % Load all weights
                contrasts.weights(i,c-col.contrastname)=w.req{r,c};
            end
            %
            i=i+1;
        end
    end
    
    % Un-request non-applicable contrasts (e.g. for models without all RL variables)
    for o2=1:1
        if strcmp(log.onsetsmodel(1:3), 'm_v')==1
            disp('Removing choice regressors')
%             varremove={'cF_Accept';'cF_Reject'; 'cF_Explore'; 'ct_NoBomb';'ct_Bomb'; 'ct_Explore'; };
            varremove={};
            if sum(strcmp(log.onsetsmodel(1:4), {'m_v1';'m_v2';'m_v3'}))==1 % vChoice
                
                % Un-request choice regressors if these are not included in the fMRI model
                for i=1:length(varremove)
                    contrasts.requested(strcmp(contrasts.contrastnames, varremove{i}))=0;
                end
            end
        elseif strcmp(log.onsetsmodel(1:3), 'm_t')==0 
            % Trialtype models do not exclude RL variables
        
            % Which RL variables to exclude?
            allvar={'EnvThreat' 0;'NTokens' 0; 'pLoss' 0; 'Entropy' 0; 'EntropyNTok' 0; 'VExplore' 0; 'Conflict' 0; 'EV' 0; 'OutcomeMean' 0; 'OutcomeVariance' 0};
            
            if isempty(request.AntsType_ifHalfDone)
            try w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']); catch; w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']); end
            else
                try w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']); catch; w.spm=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel request.AntsType_ifHalfDone ' Contrasted' filesep 'SPM.mat']); end
            end
            log.RLvariables=w.spm.SPM.RLvariables; j=1;
            for i=1:size(allvar,1)
                if isempty(cell2mat(strfind(log.RLvariables, allvar{i})))==1
                    varremove{j,1}=allvar{i}; j=j+1;
                end
            end
            w.isthere=cell(length(varremove),1);
            for i=1:length(varremove)
                w.whichtoremove=strfind(contrasts.condnames',varremove{i});
                w.isthere{i}=zeros(size(w.whichtoremove,1),1);
                for j=1:size(w.whichtoremove,1)
                    if isempty(w.whichtoremove{j})==0
                        w.isthere{i}(j)=1;
                    end
                end
            end
            
            % Un-request
            for i=1:length(varremove)
                contrasts.requested(find(w.isthere{i}))=0;
            end
            
            disp('#################################################')
            disp('[RL VARIABLES REMOVED]       Default contrasts included all RL variables. ')
            disp('List of all RL variables:'); disp(allvar); disp(' ')
            disp('The following RL variables are not included in this model:')
            disp(' '); disp(varremove); disp(' ');
            disp('Contrasts including these variables are un-requested')
            disp('#################################################')
        end
    end
    
    % Nc= No Choice-contrast models
    if isempty(strfind(log.onsetsmodel, 'TrialtypeNc'))==0
        disp('No choice contrasts are contrasted at all. Presumably because there are no choice regressors')
        varremove={'cF_Accept';'cF_Reject'; 'cF_Explore'; 'ct_NoBomb';'ct_Bomb'; 'ct_Explore'; };
        for v=1:length(varremove)        
            contrasts.requested(find(contrasts.weights(:,find(strcmp(varremove{v}, contrasts.condnames)))))=0;
        end        
    end
    
    % Model-specific un-requests! e.g. OULPEN models, ct task cannot weight Conflict, EV & NTokens
    if isempty(strfind(log.onsetsmodel, 'OULPEN'))==0
        disp('For Full-OULPEN models only: Un-weighting ct_EntropyNTok, ct_EV, ct_NTokens')
        varremove={'ct_EntropyNTok';'ct_EV';'ct_NTokens'};
        for v=1:length(varremove)        
            contrasts.requested(find(contrasts.weights(:,find(strcmp(['p' varremove{v}], contrasts.condnames)))))=0;
        end        
    end
    
    % Value models: alter names of target regressors (to match search)
    if strcmp(log.onsetsmodel(1:3), 'm_v')==1
        if strcmp(log.onsetsmodel(1:4), 'm_v5')
            contrasts.condnames(7:end)=cellfun(@(x)[x(1:6) 'x' x], contrasts.condnames(7:end),  'UniformOutput', 0);
        elseif strcmp(contrasts.condnames(find(contrasts.weights(find(contrasts.requested==1, 1, 'first'),:))), 'cF_Accept')==1;
            % Choice regressors are NOT transformed. Assume this == first 6 regressors .
            contrasts.condnames(7:end)=cellfun(@(x)[x(1:2) '_onsetx' x],  contrasts.condnames(7:end), 'UniformOutput', 0);
        elseif strcmp(log.onsetsmodel(1:5), 'm_v7_')
%             contrasts.condnames(7:end)=cellfun(@(x)[x(1:2) '_onsetx' x],  contrasts.condnames(7:end), 'UniformOutput', 0);
        elseif strcmp(contrasts.contrasttable, 'vCueOutcome')==1
            % Specify (+ update) types of outcome and cue values
            if isempty(strfind(log.onsetsmodel(1:15),'_Outcome'))==0  % Cue-Value models where both are locked to outcome onsets
                log.ValueCueType=log.onsetsmodel(4+strfind(log.onsetsmodel(4:end), '_Outcomev')+7:14+strfind(log.onsetsmodel(15:end), 'Outcome')-1);
            else log.ValueCueType=log.onsetsmodel(4+strfind(log.onsetsmodel(4:end), '_v'):strfind(log.onsetsmodel, 'Outcome')-1);
            end
            log.OutcomePEType= log.onsetsmodel(14+strfind(log.onsetsmodel(15:end), 'Outcome'):13+strfind(log.onsetsmodel(15:end), '_b'));
            contrasts.condnames{strcmp(contrasts.condnames, 'cF_CueValue')}=['cF_' log.ValueCueType];
            contrasts.condnames{strcmp(contrasts.condnames, 'ct_CueValue')}=['ct_' log.ValueCueType];
            if strcmp(log.OutcomePEType, 'OutcomePE')==0;
                contrasts.condnames{strcmp(contrasts.condnames, 'cF_OutcomePE')}=['cF_' log.OutcomePEType];
                contrasts.condnames{strcmp(contrasts.condnames, 'ct_OutcomePE')}=['ct_' log.OutcomePEType];
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'cF_OutcomePE')}=['cF_' log.OutcomePEType];
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'ct_OutcomePE')}=['ct_' log.OutcomePEType];
            end
            % Append cue and outcome onsets 
            if isempty(strfind(log.onsetsmodel(1:15),'_Outcome'))==0 
                contrasts.condnames{strcmp(contrasts.condnames, ['cF_' log.ValueCueType])}=['cF_Outcome_onsetxcF_Outcome' log.ValueCueType];
                contrasts.condnames{strcmp(contrasts.condnames, ['ct_' log.ValueCueType])}=['ct_Outcome_onsetxct_Outcome' log.ValueCueType];
            else
                contrasts.condnames{strcmp(contrasts.condnames, ['cF_' log.ValueCueType])}=['cF_onsetxcF_' log.ValueCueType];
                contrasts.condnames{strcmp(contrasts.condnames, ['ct_' log.ValueCueType])}=['ct_onsetxct_' log.ValueCueType];
            end
            contrasts.condnames{strcmp(contrasts.condnames, ['cF_' log.OutcomePEType])}=['cF_Outcome_onsetxcF_' log.OutcomePEType];
            contrasts.condnames{strcmp(contrasts.condnames, ['ct_' log.OutcomePEType])}=['ct_Outcome_onsetxct_' log.OutcomePEType];
        elseif strcmp(contrasts.contrasttable, 'vExploreVal')==1
            % What quantities? 
            if isempty(strfind(log.onsetsmodel(1:20),'VExploreInfo'))==0  % Cue-Value models where both are locked to outcome onsets
%                 error('stopped here. figure out which VExplore type, append to Contrasts table, etc');
                log.InfoValType=log.onsetsmodel(strfind(log.onsetsmodel(1:20),'VExploreInfo')+8:strfind(log.onsetsmodel(20:end),'_b')+18);
                contrasts.condnames{strcmp(contrasts.condnames, 'cF_Info')}=['cF_vExplore' log.InfoValType];  % some m_v15s may have cF_vExplore, others cF_Explore
                contrasts.condnames{strcmp(contrasts.condnames, 'ct_Info')}=['ct_vExplore' log.InfoValType];
%                 contrasts.contrastnames{strcmp(contrasts.contrastnames, 'cF_Info')}=['cF_Explore' log.InfoValType];
%                 contrasts.contrastnames{strcmp(contrasts.contrastnames, 'ct_Info')}=['ct_Explore' log.InfoValType];
            else  error('What type of Info type?');
            end
        elseif sum(strcmp(contrasts.contrasttable, {'predChoice';'vChoiceExploreOutcomePE'; 'vChoiceExploreOutcomeMagnitude'; 'vExplorePutativePEs'}))==1  
            % No need to append anything if everything is sufficiently unique
            if isempty(strfind(log.onsetsmodel(1:5),'m_v18'))==0 
                
                % The values attached to outcomes (not gamble onsets) are weighted 
                contrasts.condnames{strcmp(contrasts.condnames, 'cF_vAccept')}='cF_OutcomeOnsetxcF_vAccept';
                contrasts.condnames{strcmp(contrasts.condnames, 'cF_vReject')}='cF_OutcomeOnsetxcF_vReject';
                contrasts.condnames{strcmp(contrasts.condnames, 'cF_vExplore')}='cF_OutcomeOnsetxcF_vExplore';
                contrasts.condnames{strcmp(contrasts.condnames, 'ct_vAccept')}='ct_OutcomeOnsetxct_vAccept';
                contrasts.condnames{strcmp(contrasts.condnames, 'ct_vReject')}='ct_OutcomeOnsetxct_vReject';
                contrasts.condnames{strcmp(contrasts.condnames, 'ct_vExplore')}='ct_OutcomeOnsetxct_vExplore';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'cF_GamblevAccept')}='cF_OutcomevAccept';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'cF_GamblevReject')}='cF_OutcomevReject';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'cF_GamblevExplore')}='cF_OutcomevExplore';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'ct_GamblevAccept')}='ct_OutcomevAccept';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'ct_GamblevReject')}='ct_OutcomevReject';
                contrasts.contrastnames{strcmp(contrasts.contrastnames, 'ct_GamblevExplore')}='ct_OutcomevExplore';
                
            end
            
            
        else contrasts.condnames=cellfun(@(x)[x(1:2) '_onsetx' x],  contrasts.condnames, 'UniformOutput', 0);
        end
    end
end

 er 

disp('Requested contrasts ##########################################')
disp(contrasts.contrastnames(find(contrasts.requested)))
input('Continue to execute Contrasts?    ')

%% Execute Contrasts

for s=1:log.n_subjs
    
    disp('Running contrasts #######################################')
    disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
    
    % Organize folder for this model
    wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep];
    if sum(strcmp(log.onsetsmodel(3), {'f';'t'}))==1 % Complex renaming for the flexible models
        wb.wheremodel=[wb.where  'm_' log.onsetsmodel(3:4) '_1_'  log.onsetsmodel(max(strfind(log.onsetsmodel, '_'))+1:end)  ' Contrasted' filesep];
        if isdir(wb.wheremodel)==0;
            copyfile([wb.where log.onsetsmodel ' Estimated'], wb.wheremodel);
        end
    else
%         wb.wheremodel=[wb.where log.onsetsmodel ' Contrasted' filesep];
        wb.wheremodel=[wb.where log.onsetsmodel request.AntsType_ifHalfDone ' Contrasted' filesep];
        if isdir(wb.wheremodel)==0;
            wb.estimfol=[wb.where log.onsetsmodel ' Estimated'];
            eval('java.io.File(wb.estimfol).renameTo(java.io.File(wb.wheremodel));')
        else; disp('Found contrasts folder. Assumed correct');
        end
    end
    
    % Grab details about this model
    f   = spm_select('List', wb.wheremodel, 'SPM.mat');
    wb.spm  =load([wb.wheremodel f]);
    wb.regnamelist=wb.spm.SPM.xX.name';
    
    % Delete existing contrasts?
    if request.DeleteExistingFLContrasts && isfield(wb.spm.SPM, 'xCon')==1 && length(wb.spm.SPM.xCon)~=0
        if s==1; input('Requested deletion of all existing FL contrasts. Continue?   '); end
        matlabbatch{1}.spm.stats.con.delete = 1;
    else; matlabbatch{1}.spm.stats.con.delete = 0;
    end
    
    
    % Specify requested contrasts
    disp('Specifying Contrasts ---------------------------');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([wb.wheremodel f]); c=1;
    for k= 1:length(contrasts.contrastnames)
        if contrasts.requested(k)==1 && mean(abs(contrasts.weights(k,:)))~=0
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = contrasts.contrastnames{k};
            
            % Which regressors to actively weight?
            wc.condnums=find(contrasts.weights(k,:)~=0);
            for i=1:length(wc.condnums)
                wc.regnames{i,1}=contrasts.condnames{wc.condnums(i)};
                wc.regnames{i,2}=contrasts.condtypes{wc.condnums(i)};
            end
            [ wc.regdetails ] = f_GetRegressorNums(wb.regnamelist, wc.regnames );
            
            
            
            % Compile weights for this contrast
            wc.weights=zeros(1,length(wb.regnamelist));
            for i=1:size(wc.regdetails ,1)
                wc.weights(wc.regdetails{i,2}(1))=contrasts.weights(k,wc.condnums(i));
            end
            matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=wc.weights;
            
            % Record
            c_weights{k}=wc.weights; c_chosenregs{k}=wc.regdetails; c_names{k}=contrasts.contrastnames{k};

            wc=[]; c=c+1;
        end
    end
    
    disp('Running contrasts -----------------')
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];wb=[];
    
end

%% END

disp('================c======================================'); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ');  disp('Analysis completed:')
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' '); disp(['Model: ' log.onsetsmodel])
disp(' '); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (s7_Firstlevel_2Contrasts  -  '  log.onsetsmodel  ')']), ' ',1);
catch
end

try
edit 'D:\Dropbox\SCRIPPS\1 Explore fMRI\9 AntsAnatomical\s1_InterfaceWithAnts.m' 
end

edit 'C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI\9 AntsAnatomical\s1_InterfaceWithAnts.m'