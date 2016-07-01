% Get onsets for set up of 1st level
% clear all; close all hidden; path(pathdef); clc;
clear all; path(pathdef); clc;

% Requested analysis
logg.specificsubjects={};
logg.specificsubjects={'p01_GV'};

% logg.specificsubjects={'p06_KB'};
% logg.specificsubjects={'p02_YY';'p04_MP';'p06_KB';'p08_SG';'p10_RC';'p13_HL';'p15_SH';'p17_SJ';'p18_MS';'p21_ES';'p23_BS';'p25_RJ';'p27_DM';'p30_KL';'p34_TB';'p35_SM';'p36_FR';'p38_MK';'p41_AC';};

logg.FLthread=' s4Ants'; logg.func_prefix='s4ubf';
logg.warpJ=0;  % Value models already have j warped - dont need to turn this on
request.loadsimplifieddata=[]; % '(10-Aug-2013) Simplified data for onsets';  % Blank to process data from scratch

% Which onsets model? ########################
for o1=1:1 % Standard models 
% logg.onsetsmodel='m_c1_Choice_ENU'; 
% logg.onsetsmodel='m_c2_Choice_ENV'; 
% logg.onsetsmodel='m_c3_ChoiceFull_OULPEN'; 
% logg.onsetsmodel='m_c4_Choice_OUPEN'; 
% logg.onsetsmodel='m_c5_ChoiceRTFull_OULPEN'; 
% logg.onsetsmodel='m_c6_ChCluster4Full_OULPEN'; 
% logg.onsetsmodel='m_c7_ChCluster6Full_OULPEN'; 
% logg.onsetsmodel='m_c8_ChCluster4MovFull_OULPEN'; 
% logg.onsetsmodel='m_c9_ChCluster6MovFull_OULPEN'; 
% logg.onsetsmodel='m_c10_ChCluster6FullRT_OULPEN'; 
% logg.onsetsmodel='m_c11_Choice_VOUPEN'; 
% logg.onsetsmodel='m_c12_Choice_VOULPEN'; 
% logg.onsetsmodel='m_o1_OrthogBasic';
end
for o1=1:1 % Flex/Trialtype models
% logg.onsetsmodel='m_f1_ChoicexTrialtype';
% logg.onsetsmodel='m_f2_ChunkChoicexTrialtype';
% logg.onsetsmodel='m_f3_Event';
% logg.onsetsmodel='m_t1_Trialtype';
% logg.onsetsmodel='m_t2_TrialtypeNc';
% logg.onsetsmodel='m_t3_ChunkTrialtype';
% logg.onsetsmodel='m_t4_ChunkTrialtypeNc';
end
for o1=1:1 % Value/behavioural modelling models
% logg.onsetsmodel='m_v1c_vChoice_bpm16bpmi11';
% logg.onsetsmodel='m_v1e_vChoice_b01b01';
% logg.onsetsmodel='m_v3c_vChosenAnd_bpm16bpmi11';
% logg.onsetsmodel='m_v2_vBestAnd_bpmi16bpmi11';
% logg.onsetsmodel='m_v3_vChosenAnd_bpmi16bpmi11';
% logg.onsetsmodel='m_v4_ChoicevChosenAnd_bpmi16bpmi11';
% logg.onsetsmodel='m_v4b_ChoicevChosenAnd_bpmi16bpmi11'; % TypicalFit
% logg.onsetsmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11';  % Mixed models, based on hierarchical fit results
% logg.onsetsmodel='m_v4d_ChoicevChosenAnd_bpm16bpm11';  % Mixed models!
% logg.onsetsmodel='m_v4c_ChoicevChosenAnd_bpm16bpmi11';  % Mixed models!
% logg.onsetsmodel='m_v5c_ChoiceXvChosenAnd_bpm16bpmi11';
% logg.onsetsmodel='m_v6_ChoicevChosenAndposneg2_bpmi16bpmi11';
% logg.onsetsmodel='m_v6c_ChoicevChosenAndposneg2_bpm16bpmi11';
% logg.onsetsmodel='m_v6e_ChoicevChosenAndposneg2_b01b01';
% logg.onsetsmodel='m_v7_ChClustervChosenAnd_bpmi16bpmi11';
% logg.onsetsmodel='m_v8c_ChoicesubEVposneg_bpm16bpmi11';
% logg.onsetsmodel='m_v9c_vChosenAndposneg2_bpm16bpmi11';
% logg.onsetsmodel='m_v10c_RejectOrvChosenAndposneg2_bpm16bpmi11';
% logg.onsetsmodel='m_v11c_vGambleOutcomePE_bpm16bpmi11';  % OutcomePE
% logg.onsetsmodel='m_v11e_vGambleOutcomePE_b01b01';  % OutcomePE
% logg.onsetsmodel='m_v11c_vGambleOutcomeMag_bpm16bpmi11';    % OutcomeMag
% logg.onsetsmodel='m_v12c_vModalchoiceOutcomePE_bpm16bpmi11';  % OutcomePE
% logg.onsetsmodel='m_v12c_vModalchoiceOutcomeMag_bpm16bpmi11';    % OutcomeMag
% logg.onsetsmodel='m_v13c_OutcomevGambleOutcomeMagnitude_bpm16bpmi11';
% logg.onsetsmodel='m_v14c_OutcomevModalchoiceOutcomeMagnitude_bpm16bpmi11';
% logg.onsetsmodel='m_v15c_VExploreInfoPE_bpm16bpmi11';
% logg.onsetsmodel='m_v15c_VExploreInfoPEsign_bpm16bpmi11';
% logg.onsetsmodel='m_v15c_VExploreInfoVal_bpm16bpmi11';
% logg.onsetsmodel='m_v16c_vChoicePE_bpm16bpmi11';
% logg.onsetsmodel='m_v17c_vChoiceOutcomeMag_bpm16bpmi11';
% logg.onsetsmodel='m_v18c_vChoiceAtOutcomeMag_bpm16bpmi11';
% logg.onsetsmodel='m_v19c_VExploreInfoPEOutcomePE_bpm16bpmi11';
% logg.onsetsmodel='m_v20c_pExplore_bpm16bpmi11';
% logg.onsetsmodel='m_v21c_ExploreGamInfoOutcomePE_bpm16bpmi11';
% logg.onsetsmodel='m_v21e_ExploreGamInfoOutcomePE_b01b01';
% logg.onsetsmodel='m_v22c_OutcomevChosenOutcomeMagnitude_bpm16bpmi11';
% logg.onsetsmodel='m_v22e_OutcomevChosenOutcomeMagnitude_b01b01';
% logg.onsetsmodel='m_v22f_OutcomevChosenOutcomeMagnitude_b02b01';
% logg.onsetsmodel='m_v23c_OutcomevChosenOutcomeMagnitude_bpm16bpmi11';
% logg.onsetsmodel='m_c1_Choice_ENU';
% logg.onsetsmodel='m_v23c_RejectOrvChosenAnd_bpm16bpmi11';
% logg.onsetsmodel='m_v1g_vChoice_bpji08bpji11';  
% logg.onsetsmodel='m_c13g_ChoiceFull_ULPEN';
% logg.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% logg.onsetsmodel='m_v9g_vChosenAndposneg_bpji08bpji11' ; 
% logg.onsetsmodel='m_v23g_RejectOrvChosenAnd_bpji08bpji11';
% logg.onsetsmodel='m_v6g_ChoicevChosenAndposneg2_bpji08bpji11';
% logg.onsetsmodel='m_v24g_predChoice_bpji08bpji11';
% logg.onsetsmodel='m_v25g_RejectOrvGamble_bpji08bpji11';
% logg.onsetsmodel='m_c14_Choice';
% logg.onsetsmodel='m_v26e_pvBestChoice_b01b01';
% logg.onsetsmodel='m_c15g_NextChoice_ULPEN';
% logg.onsetsmodel='m_v4g_ChoicevChosenAnd_bpji08bpji11';
% logg.onsetsmodel='m_v27g_ChoiceRejectOrvChosenAnd_bpji08bpji11';
% logg.onsetsmodel='m_v28g_ChoicePredChoice_bpji08bpji11';
% logg.onsetsmodel='m_v29g_ChoiceRejectOrvGamble_bpji08bpji11';
% logg.onsetsmodel='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11';
% logg.onsetsmodel='m_v31g_ChoicecFRejectOrvChosenAnd_bpji08bpji11';
% logg.onsetsmodel='m_c16_ChCluster4';
% logg.onsetsmodel='m_c17_ChCluster6';
% logg.onsetsmodel='m_v32g_vGamble_bpji08bpji11';
% logg.onsetsmodel='m_v33g_ChoicevGamble_bpji08bpji11';
% logg.onsetsmodel='m_v8g_vGamblePosNeg_bpji08bpji11';
% logg.onsetsmodel='m_v34g_ChoiceXvGamble_bpji08bpji11';
% logg.onsetsmodel='m_v35g_ChoiceChoiceXvGamble_bpji08bpji11';
% logg.onsetsmodel='m_c18_pChoice';
% logg.onsetsmodel='m_c19_pChoiceFull_ULPEN';
% logg.onsetsmodel='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11';
% logg.onsetsmodel='m_v3g_vChosenAnd_bpji08bpji11'; 
% logg.onsetsmodel='m_v36g_RejectOrvMargChoDiff_bpji08bpji11';
% logg.onsetsmodel='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11';
% logg.onsetsmodel='m_v38g_EVGainLoss_bpji08bpji11';
% logg.onsetsmodel='m_v39g_ChoiceEVGainLoss_bpji08bpji11';
% logg.onsetsmodel='m_v40g_vBestvWorst_bpji08bpji11'; 
% logg.onsetsmodel='m_v41g_ChoicevBestvWorst_bpji08bpji11'; 
% logg.onsetsmodel='m_v42g_pLossNTok_bpji08bpji11'; 
% logg.onsetsmodel='m_v43g_ChoicepLossNTok_bpji08bpji11'; 
% logg.onsetsmodel='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11';
% logg.onsetsmodel='m_v45g_ChoicevBUposnegvMargChoDiff_bpji08bpji11';
% logg.onsetsmodel='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11'; 
end
% logg.onsetsmodel='m_c6_ChCluster4Full_ULPEN'; 
% logg.onsetsmodel='m_c7_ChCluster6Full_ULPEN'; 
% logg.onsetsmodel='m_c21_RejExpFull_ULPEN'; 
% logg.onsetsmodel='m_c22_ChoiceFull_HLPEN'; 
% logg.onsetsmodel='m_c23_ChoiceFullRT_ULPEN'; 
% logg.onsetsmodel='m_c24_ChCluster6RT_ULPEN'; 
logg.onsetsmodel='m_c13g_ChoiceFull_ULPEN';

for o1=1:1 % General settings and specifications
    
    % Add paths
    w.w=pwd; if strcmp(w.w(1), '/')==0;  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; where.expt_fol='G:\2 [Explore]'; 
         
        where.data_brain='G:\2 [Explore]\1 Brain data'; where.data_beh=[where.where filesep '1 Behavioural data'];   
        where.modvals='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis\4 Fit computational models\1 Value functions';
    else where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI';  where.expt_fol='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; where.data_brain='/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data'; where.data_beh=[where.where filesep '1 Behavioural data']; where.modvals='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis/4 Fit computational models/1 Value functions';
    end
    path(pathdef); addpath(where.modvals)
    addpath(where.where);  addpath([where.where filesep '3 Scripts - Preprocessing']);
    addpath(genpath([where.where filesep '4 Set up models']));
    cd(where.where); cd ..; w.here=pwd; where.analysisfolder=[w.here filesep '2 Explore experiment' filesep '3 Analysis' filesep '4 Fit computational models']; cd(where.where)
    addpath(where.analysisfolder);
    
    % Load subjects    
    logg.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); logg.datalog=logg.w.datalog;
    [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.datalog, logg.specificsubjects, [logg.datalog vertcat('include_all', num2cell(ones(size(logg.datalog,1)-1,1)))], 'include_all');
    
    % Apply further subject selection for some models
    w.modelsneedingsubselect={'m_c6_';'m_c7_';'m_v7_'; 'm_c24';};
    if sum(strcmp(logg.onsetsmodel(1:5), w.modelsneedingsubselect))==1
        [w.s w.s1 logg.koshertable]=xlsread(['i_Subjectdataok_SpecificModels.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.datalog, logg.specificsubjects, logg.koshertable,logg.onsetsmodel);
    end
    
    % Requests & Details that don't change much
    request.saveonsets=1;
    scan=load('i_scanningdetails.mat'); scan.TRseconds=scan.TRms/1000; 
    errorlog=cell(1,1); e=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0; disp('   Subset of subjects only:'); disp(logg.subjects); end
    disp(request); disp(' ')
    disp(['Data location: ' where.data_brain])
    disp('=======================================================')
    
end

%% (1) Format data

for o1=1:1 % Data specifications (columns)
    
    col.Onset_Offer=22; % Event onsets
    col.Onset_Info=25;
    col.Onset_Outcome=28;
    col.Onset_TrialEnd=29;
    col.Onset_Motor1=23;
    col.Onset_Motor2=26;
    col.Duration_Offer=31; % Event Durations
    col.Duration_Info=32;
    col.Duration_Outcome=33;
    col.Duration_Trial2End=34;
    col.Duration_Choice=35;
    
    col.Task=5; % Event classifiers
    col.Resp1=8; 
    col.Resp2=10;
    col.TrialValid=13; 
    col.OutcomePresented=14;
    col.OutcomeMagnitude=15;
    col.RT1=9;
    col.Block=17;
    col.ShowedExploredBomb=18;

    col.EnvThreat=3; % Design columns
    col.NTokens=2;
    col.pLoss=36;
    col.Entropy=37;
    col.VExplore=38;
    col.EntropyNTok=42;
    col.EV=41;
    col.OutcomeMean=39;
    col.OutcomeVariance=40;
    col.TrialType=1;
    col.pAccept=41;
    col.pReject=42;
    col.pExplore=43;
    col.ChoiceH=structmax(col)+1; 
    
end 

if isempty(request.loadsimplifieddata)==1
    disp('Fetching data + marking variables ############')
    subjdata=cell(logg.n_subjs,2); w.blockoffsets=zeros(logg.n_subjs,6); w.blockoffsets(:,1)=0; w.offsetfromstarts=zeros(logg.n_subjs,6);
    for s=1:logg.n_subjs
        disp(['Subject ' num2str(s) '   -  ' logg.subjects{s}])
        subjdata{s,1}=logg.subjects{s};
        
        for r=1: 6 % Load data + adjust timings
            ws.r=load([where.data_beh filesep logg.subjects{s} filesep logg.subjects{s} '_file_6integratedfMRI_b' num2str(r) '.mat']);
            w.offsetfromstarts(s,r)=ws.r.integratedfMRI.times.start/1000; % Start time in seconds
            ws.r=ws.r.integratedfMRI.rdata;
            ws.r(:,col.Block)=r;
            
            % Offset from start of scans (run offsets)
            f=spm_select('List',[where.data_brain filesep logg.subjects{s} filesep '1 Preprocessed' filesep 'Func_r' num2str(r)], ['^' logg.func_prefix '*.*img$']);
            if isempty(f)==1
                f=spm_select('List',[where.data_brain filesep logg.subjects{s} filesep '2 First level' logg.FLthread filesep 'Preproc functionals'], ['^' logg.func_prefix logg.datalog{s+1,4} '-000' num2str(logg.datalog{s+1,4+r}) '.*.*img$']);
                if isempty(f)==1
                    f=spm_select('List',[where.data_brain filesep logg.subjects{s} filesep '2 First level' logg.FLthread filesep 'Preproc functionals'], ['^' logg.func_prefix logg.datalog{s+1,4} '-00' num2str(logg.datalog{s+1,4+r}) '.*.*img$']);
                end
            end
            if isempty(f)==1
                errorlog{e}=['ERROR: Could not count number of scans - Arbitrary block offsets used     (' logg.subjects{s} ')']; disp(errorlog{e});e=e+1;
                if s==1 && r==1, input('Continue with arbitrary offsets?'); end
%                 warning(errorlog{e})
                f=zeros(255,1);
            else
                w.nscans(s,r)=size(f,1);
                w.blockoffsets(s,r+1)=w.blockoffsets(s,r)+size(f,1)*scan.TRms*scan.nSlicesPerVol/1000; % accumulative offset in seconds
            end
            
            % Adjust onsets accordingly
            ws.r(:,col.Onset_Motor1)=ws.r(:,col.Onset_Motor1)/1000; % Convert keytimes from ms to s
            ws.r(:,col.Onset_Motor2)=ws.r(:,col.Onset_Motor2)/1000;
            ws.r(:,col.Onset_TrialEnd)=ws.r(:,col.Onset_TrialEnd)/1000;
            ws.r(:,[20 21 22 23 25 26 28 29])=ws.r(:,[20 21 22 23 25 26 28 29])-w.offsetfromstarts(s,r)+w.blockoffsets(s,r);

            % Mark durations
            ws.r(ws.r(:,col.OutcomePresented)==0,col.Duration_Offer)=  ws.r(ws.r(:,col.OutcomePresented)==0, col.Onset_TrialEnd)    -    ws.r(ws.r(:,col.OutcomePresented)==0, col.Onset_Offer);
            ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3,col.Duration_Offer)=ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3, col.Onset_Info)    -       ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3, col.Onset_Offer);
            ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)~=3,col.Duration_Offer)=ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)~=3,col.Onset_Outcome)     -    ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)~=3,col.Onset_Offer);
            ws.r(:,col.Duration_Info)=0;
            ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3,col.Duration_Info)=ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3,col.Onset_Outcome)    -      ws.r(ws.r(:,col.OutcomePresented)==1 & ws.r(:,col.Resp1)==3,col.Onset_Info);
            ws.r(:,col.Duration_Outcome)=0;
            ws.r(ws.r(:,col.OutcomePresented)==1,col.Duration_Outcome)=ws.r(ws.r(:,col.OutcomePresented)==1,col.Onset_TrialEnd )   -         ws.r(ws.r(:,col.OutcomePresented)==1,  col.Onset_Outcome);
            ws.r(:,col.Duration_Choice)=0;
            ws.r(ws.r(:,col.TrialValid)==1, col.Duration_Choice)=ws.r(ws.r(:,col.TrialValid)==1, col.Duration_Offer)-ws.r(ws.r(:,col.TrialValid)==1, col.RT1)/1000; % Choice duration=Time between choice & end of 1st offer
            
            %
            subjdata{s,2}.run{r}=ws.r;
            eval(['ws.r' num2str(r) '=ws.r;']);
            ws=[];
        end
        
        % Combine across runs + mark RL variables 
        subjdata{s,2}.d_conflictcontrol=vertcat(subjdata{s,2}.run{1},subjdata{s,2}.run{2},subjdata{s,2}.run{3},subjdata{s,2}.run{4},subjdata{s,2}.run{5},subjdata{s,2}.run{6});
        [ subjdata{s,2}.d_conflict] = fpar_conflict(subjdata{s,2}.d_conflictcontrol(subjdata{s,2}.d_conflictcontrol(:,col.Task)==1,:),col); 
        subjdata{s,2}.d_conflict=sortrows(subjdata{s,2}.d_conflict, col.Onset_Offer);
        [ subjdata{s,2}.d_control] = fpar_control(subjdata{s,2}.d_conflictcontrol(subjdata{s,2}.d_conflictcontrol(:,col.Task)==2,:),col); subjdata{s,2}.d_control=sortrows(subjdata{s,2}.d_control, col.Onset_Offer);
        
        % % Paranoid Control parameter specifiction: Only EV changes. Technically, fpar_control does this right already 
        % [ subjdata{s,2}.d_control] = fpar_conflict(subjdata{s,2}.d_conflictcontrol(subjdata{s,2}.d_conflictcontrol(:,col.Task)==2,:),col); subjdata{s,2}.d_control=sortrows(subjdata{s,2}.d_control, col.Onset_Offer);
        % [ subjdata{s,2}.d_control_real] = fpar_control(subjdata{s,2}.d_conflictcontrol(subjdata{s,2}.d_conflictcontrol(:,col.Task)==2,:),col); subjdata{s,2}.d_control=sortrows(subjdata{s,2}.d_control, col.Onset_Offer);
        % subjdata{s,2}.d_control(:,col.EV)=subjdata{s,2}.d_control_real(:,col.EV);
        
        % Observed p(Choice)
        for e=1:6
            for n=1:6
                wc.dcf=subjdata{s,2}.d_conflict( subjdata{s,2}.d_conflict(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_conflict(:, col.NTokens)./2==n, col.Resp1);
                subjdata{s,2}.d_conflict( subjdata{s,2}.d_conflict(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_conflict(:, col.NTokens)./2==n, col.pAccept) = mean(wc.dcf==1);
                subjdata{s,2}.d_conflict( subjdata{s,2}.d_conflict(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_conflict(:, col.NTokens)./2==n, col.pReject) = mean(wc.dcf==2);
                subjdata{s,2}.d_conflict( subjdata{s,2}.d_conflict(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_conflict(:, col.NTokens)./2==n, col.pExplore) = mean(wc.dcf==3);
                
                % Choice entropy 
                wc.pch = [0 0 0]; 
                for ch=1:3, wc.pch(ch) =  mean( wc.dcf==ch);  end
                wc.pch(wc.pch==1) = 1-eps; wc.pch(wc.pch==0) = 0+eps; 
                subjdata{s,2}.d_conflict( subjdata{s,2}.d_conflict(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_conflict(:, col.NTokens)./2==n, col.ChoiceH)= sum( (-wc.pch).*log(wc.pch)); 
                wc=[];
                %
                wc.dct=subjdata{s,2}.d_control( subjdata{s,2}.d_control(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_control(:, col.NTokens)./2==n, col.Resp1);
                subjdata{s,2}.d_control( subjdata{s,2}.d_control(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_control(:, col.NTokens)./2==n, col.pAccept) = mean(wc.dct==1);
                subjdata{s,2}.d_control( subjdata{s,2}.d_control(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_control(:, col.NTokens)./2==n, col.pReject) = mean(wc.dct==2);
                subjdata{s,2}.d_control( subjdata{s,2}.d_control(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_control(:, col.NTokens)./2==n, col.pExplore) = mean(wc.dct==3);
                
                % Choice entropy 
                wc.pch = [0 0 0]; 
                for ch=1:3, wc.pch(ch) =  mean( wc.dct==ch);  end
                wc.pch(wc.pch==1) = 1-eps; wc.pch(wc.pch==0) = 0+eps; 
                subjdata{s,2}.d_control( subjdata{s,2}.d_control(:, col.EnvThreat).*6==e &  subjdata{s,2}.d_control(:, col.NTokens)./2==n, col.ChoiceH)= sum( (-wc.pch).*log(wc.pch)); 
                wc=[];
            end
        end
        
        
        % Combine tasks 
        subjdata{s,2}.d_conflictcontrol=vertcat( subjdata{s,2}.d_conflict,  subjdata{s,2}.d_control);
        subjdata{s,2}.d_conflictcontrol=sortrows(subjdata{s,2}.d_conflictcontrol,col.Onset_Offer);
%         ws=subjdata{s,2}.d_conflictcontrol; 
%         ws= ws(:,[col.Block  col.Onset_Offer col.Onset_Outcome col.Onset_TrialEnd]);
       
        % Additional?
%         datanext{1}(:, col.Resp1) = datanext{1}(:, col.NextChoice);
%         subjdata{s,2}.d_conflictcontrol(:, [col.Task col.Block])
        
        subjdata{s,2}.d{1}=subjdata{s,2}.d_conflict; % d: holds all data (1=conflict, 2=control, 3=all)
        subjdata{s,2}.d{2}=subjdata{s,2}.d_control;
        subjdata{s,2}.d{3}=subjdata{s,2}.d_conflictcontrol;
    end
    
    
    
    
    % Save simplified data
%     save([where.data_beh filesep '(' date ') Simplified data for onsets.mat'], 'subjdata');
else load([where.data_beh filesep request.loadsimplifieddata])
end

%% (2) Add additional variables 

% Append data from behvioural modelling (if necessary)
if sum(strcmp(logg.onsetsmodel(1:5), {'m_c20'}))==1 ||    strcmp(logg.onsetsmodel(1:3), 'm_v')==1 || logg.warpJ
     request.behmod_fit=logg.onsetsmodel(strfind(logg.onsetsmodel, '_b')+1:end);
     if isempty(request.behmod_fit)==1
         
         % HARD CODE
         if isempty(strfind(logg.onsetsmodel(1:8), 'g_'))==0;
             request.behmod_fit    = 'bpji08bpji11';
         else error('Which behvioural models thread?'); 
         end
         disp(['Behavioural model thread assumed to be: '         request.behmod_fit    ]);
         input('Continue?  '); 
     end
     where.behmodfile=[where.expt_fol filesep '2 Second level results' filesep '2 Behavioural model details' filesep request.behmod_fit filesep];
     f=spm_select('List', where.behmodfile, 'Model values*');
     if size(f,1)~=1; error('Multiple model-value files! Which? Move others to archive'); end
     behmod_pars=load([where.behmodfile filesep f]);
end
if strcmp(logg.onsetsmodel(1:3), 'm_v')==1 || sum(strcmp(logg.onsetsmodel(1:5), {'m_c20'}))==1
     % New cols (check no overlap with existing)
     col.vAccept=structmax(col)+1; % v(Accept/Reject/Explore)
     col.vReject=structmax(col)+1;
     col.vExplore=structmax(col)+1;
     col.vChoice=[col.vAccept col.vReject col.vExplore];
     col.vBest=structmax(col)+1;  % v(Best option/2nd best)
     col.vSecondBest=structmax(col)+1;
     col.vBestAnd=[col.vBest col.vSecondBest];
     col.vChosen=structmax(col)+1;
     col.vBestUnchosen=structmax(col)+1;
     col.vBestUnchosen_pos=structmax(col)+1;
     col.vBestUnchosen_neg=structmax(col)+1;     
     col.vGamble=structmax(col)+1;   % Value of the gamble: What can be immediately gained, without explore. V(Accept) for cF, V(NonExplore) for ct. [cF: pLoss*SubjectiveLoss + (1-pLoss)*NTok] [ct: (1-Entropy)*NTok]
     col.vGamblePos=structmax(col)+1;  
     col.vGambleNeg=structmax(col)+1; 
     col.Rej_vBestUnchosen=structmax(col)+1;
     col.NonRej_vBestUnchosen=structmax(col)+1;
     col.PEoutcome=structmax(col)+1;
     col.vModalchoice=structmax(col)+1;
     col.vSee=structmax(col)+1;
     col.vNoSee=structmax(col)+1;
     col.vExploreInfo=structmax(col)+1;
     col.ExploreInfoPE=structmax(col)+1;
     col.ExploreInfoPEsign=structmax(col)+1;
     col.AcceptOutcomePE=structmax(col)+1;
     col.RejectOutcomePE=structmax(col)+1;
     col.ExploreInfotoOutcomePE=structmax(col)+1;
     col.AcceptOutcomeMagnitude=structmax(col)+1;
     col.RejectOutcomeMagnitude=structmax(col)+1;
     col.ExploreOutcomeMagnitude=structmax(col)+1;
     col.predAccept=structmax(col)+1;
     col.predReject=structmax(col)+1;
     col.predExplore=structmax(col)+1;
     col.predChoice=[col.predAccept col.predReject col.predExplore];
     col.ExploreGambletoOutcomePE=structmax(col)+1; 
     col.Rej_vGamble=structmax(col)+1;
     col.NonRej_vGamble=structmax(col)+1;
     col.pvBest=structmax(col)+1;
     col.Acc_vGamble=structmax(col)+1;
     col.Exp_vGamble=structmax(col)+1;
     col.vMargChosen=structmax(col)+1;   % Marginal VChosen: vCho>vBU 
     col.Rej_vMargChosen=structmax(col)+1;
     col.NonRej_vMargChosen=structmax(col)+1;
     col.EVGain=structmax(col)+1;
     col.EVLoss=structmax(col)+1;
     col.EVConflict=structmax(col)+1;
     col.vBest=structmax(col)+1;
     col.vWorst=structmax(col)+1;
     col.vBesttoWorst=structmax(col)+1;
     col.pLossNTok=structmax(col)+1;
     col.vBUpos_vMargChosen=structmax(col)+1;
     col.vBUneg_vMargChosen=structmax(col)+1;
     col.Acc_vMargChosen=structmax(col)+1; 
     col.Exp_vMargChosen=structmax(col)+1;
    
     for os=1:logg.n_subjs
         ws.d={subjdata{os,2}.d_conflict subjdata{os,2}.d_control}; s=[];
         ss=find(strcmp(behmod_pars.d_vchoice(:,1), logg.subjects{os}));  % Index of subjject in behmodpars !!!!!!!!
         
         % (a) Load values from matrix-formatted behmod pars
         %      Find right value for each trial-type, and append to data (as additional column)
         for task=1:2
             for e=1:6  % Note: e is coded as ascendingly-ordered rank, NOT in terms of imagesc
                 for n=1:6
                     
                     % V(Accept/Reject/Explore)
                     for c=1:3
                         ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vChoice(c))  =  behmod_pars.d_vchoice{ss, task+1}{c}(e,n);
                     end
                     
                     % V(Best option/2nd best)  
                     for c=1:2
                          ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vBestAnd(c))   =  behmod_pars.d_vbestchoice{ss, task+1}{c}(e,n);
                     end
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vBest)   =   behmod_pars.d_vbestchoice{ss, task+1}{1}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vWorst)   =   behmod_pars.d_vbestchoice{ss, task+1}{3}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vBesttoWorst)   =   behmod_pars.d_vbestchoice{ss, task+1}{1}(e,n) - behmod_pars.d_vbestchoice{ss, task+1}{3}(e,n);
                     
                     % Prediction p(Accept/Reject/Explore)
                     for c=1:3
                         ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.predChoice(c))  =   behmod_pars.d_predchoice{ss, task+1}{c}(e,n);
                     end
                     
                     % EV Gain/Loss
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.EVGain) =  behmod_pars.d_evgainloss{ss, task+1}{1}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.EVLoss) =  behmod_pars.d_evgainloss{ss, task+1}{2}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.EVConflict) =  behmod_pars.d_evgainloss{ss, task+1}{3}(e,n);
                     
                     % PavConflict
                     ws.pn =ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,[col.pLoss col.NTokens]); 
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,col.pLossNTok)=  ws.pn(1,1)*ws.pn(1,2); 
                     
                     % Append here to add more variables from the behavioural models.
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n,  col.vModalchoice)=   behmod_pars.d_modalchoice_andval{ss,task+1}{2}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n, col.vSee)= behmod_pars.d_vseenosee{ss,task+1}{1}(e,n);
                     ws.d{task}(ws.d{task}(:,  col.EnvThreat)*6==e & ws.d{task}(:,  col.NTokens)/2==n, col.vNoSee)= behmod_pars.d_vseenosee{ss,task+1}{2}(e,n);
                    
                 end
             end
         end
         
         % (b) Load values from trialstats-formatted behmod pars
         for t=1:2
             ws.d{t}(:, [col.vChosen col.vBestUnchosen])= behmod_pars.d_trialstats{ss,t+1}(:, [behmod_pars.scol.vChosen behmod_pars.scol.vBestUnchosen]);
             ws.d{t}(:, [col.vBestUnchosen_pos  col.vBestUnchosen_neg])=behmod_pars.d_trialstats{ss,t+1}(:, [behmod_pars.scol.vBestUnchosen_pos behmod_pars.scol.vBestUnchosen_neg]);
             ws.d{t}(:, col.vGamble)=behmod_pars.d_trialstats{ss,t+1}(:, [behmod_pars.scol.vGamble]);
             ws.d{t}(:, [col.vGamblePos col.vGambleNeg])=0;              
             ws.d{t}(ws.d{t}(:, col.vGamble)>0, col.vGamblePos ) = ws.d{t}(ws.d{t}(:, col.vGamble)>0, col.vGamble);
             ws.d{t}(ws.d{t}(:, col.vGamble)<0, col.vGambleNeg) = ws.d{t}(ws.d{t}(:, col.vGamble)<0, col.vGamble);
             ws.d{t}(:, [col.Rej_vBestUnchosen  col.NonRej_vBestUnchosen])=nan;
             ws.d{t}( ws.d{t}(:, col.Resp1)==2, col.Rej_vBestUnchosen)=ws.d{t}( ws.d{t}(:, col.Resp1)==2, col.vBestUnchosen);
             ws.d{t}( ws.d{t}(:, col.Resp1)~=2, col.NonRej_vBestUnchosen)=ws.d{t}( ws.d{t}(:, col.Resp1)~=2, col.vBestUnchosen);
             ws.d{t}(:, col.PEoutcome)=ws.d{t}(:, col.OutcomeMagnitude)-ws.d{t}(:, col.vChosen);
             ws.d{t}(:, col.AcceptOutcomePE)=nan;  ws.d{t}(:, col.AcceptOutcomeMagnitude)=nan;
             ws.d{t}( ws.d{t}(:, col.Resp1)==1,col.AcceptOutcomePE) = ws.d{t}( ws.d{t}(:, col.Resp1)==1,col.PEoutcome);
             ws.d{t}( ws.d{t}(:, col.Resp1)==1,col.AcceptOutcomeMagnitude) = ws.d{t}( ws.d{t}(:, col.Resp1)==1,col.OutcomeMagnitude);
             ws.d{t}(:, col.RejectOutcomePE)=nan; ws.d{t}(:, col.RejectOutcomeMagnitude)=nan;
             ws.d{t}( ws.d{t}(:, col.Resp1)==2,col.RejectOutcomePE) = ws.d{t}( ws.d{t}(:, col.Resp1)==2,col.PEoutcome);
             ws.d{t}( ws.d{t}(:, col.Resp1)==2,col.RejectOutcomeMagnitude) = ws.d{t}( ws.d{t}(:, col.Resp1)==2,col.OutcomeMagnitude);            
             ws.d{t}(:, [col.ExploreOutcomeMagnitude col.Rej_vGamble col.Acc_vGamble col.Exp_vGamble col.NonRej_vGamble col.Acc_vGamble col.Rej_vMargChosen col.NonRej_vMargChosen col.Exp_vGamble  col.vBUpos_vMargChosen  col.vBUneg_vMargChosen])=nan; 
             ws.d{t}( ws.d{t}(:, col.Resp1)==3, col.ExploreOutcomeMagnitude) = ws.d{t}( ws.d{t}(:, col.Resp1)==3,col.OutcomeMagnitude);
             ws.d{t}( ws.d{t}(:, col.Resp1)==1, col.Acc_vGamble)=ws.d{t}( ws.d{t}(:, col.Resp1)==1, col.vGamble);
             ws.d{t}( ws.d{t}(:, col.Resp1)==2, col.Rej_vGamble)=ws.d{t}( ws.d{t}(:, col.Resp1)==2, col.vGamble);
             ws.d{t}( ws.d{t}(:, col.Resp1)==3, col.Exp_vGamble)=ws.d{t}( ws.d{t}(:, col.Resp1)==3, col.vGamble);
             ws.d{t}( ws.d{t}(:, col.Resp1)~=2, col.NonRej_vGamble)=ws.d{t}( ws.d{t}(:, col.Resp1)~=2, col.vGamble);
             ws.v=ws.d{t}( :, [col.vAccept  col.vReject  col.vExplore]); ws.v(:,4)= max(ws.v,[], 2);
             ws.b=1; ws.d{t}(:, col.pvBest) = exp(ws.b*ws.v(:,4)) ./(exp(ws.b*ws.v(:,1))+exp(ws.b*ws.v(:,2))+exp(ws.b*ws.v(:,3)));
             ws.d{t}(:, col.vMargChosen)=ws.d{t}(:, col.vChosen) - ws.d{t}(:,col.vBestUnchosen);
             ws.d{t}(ws.d{t}(:, col.Resp1)==1, col.Acc_vMargChosen) =ws.d{t}( ws.d{t}(:, col.Resp1)==1,col.vMargChosen);
             ws.d{t}(ws.d{t}(:, col.Resp1)==2, col.Rej_vMargChosen) =ws.d{t}( ws.d{t}(:, col.Resp1)==2,col.vMargChosen);
             ws.d{t}(ws.d{t}(:, col.Resp1)==3, col.Exp_vMargChosen) =ws.d{t}( ws.d{t}(:, col.Resp1)==3,col.vMargChosen);
             ws.d{t}(ws.d{t}(:, col.Resp1)~=2, col.NonRej_vMargChosen) = ws.d{t}( ws.d{t}(:, col.Resp1)~=2,col.vMargChosen);
             ws.d{t}( ws.d{t}(:, col.vBestUnchosen)<0, col.vBUneg_vMargChosen) = ws.d{t}( ws.d{t}(:, col.vBestUnchosen)<0, col.vMargChosen);
             ws.d{t}( ws.d{t}(:, col.vBestUnchosen)>=0, col.vBUpos_vMargChosen) = ws.d{t}( ws.d{t}(:, col.vBestUnchosen)>=0, col.vMargChosen);
             
             % Explore PE
             ws.expshown_trials=ws.d{t}(:, col.Resp1)==3 & ws.d{t}(:, col.OutcomePresented)==1;
             ws.d{t}(:, [col.vExploreInfo col.ExploreInfoPE])=nan;
             ws.d{t}(:, col.ExploreInfoPE)=1:size(ws.d{t},1);
             ws.dexpinfo{t}=ws.d{t}(ws.d{t}(:, col.Resp1)==3 & ws.d{t}(:, col.OutcomePresented)==1, :);
             ws.dexpinfo_see{t}=ws.dexpinfo{t}(   ws.dexpinfo{t}(:, col.ShowedExploredBomb)==1, :);
             ws.dexpinfo_nosee{t}=ws.dexpinfo{t}(   ws.dexpinfo{t}(:, col.ShowedExploredBomb)==0, :);
             logg.n_exploreinfo(os,(t-1)*3+1:(t-1)*3+3)=[size(ws.dexpinfo{t},1) size(ws.dexpinfo_see{t},1) size(ws.dexpinfo_nosee{t},1)];
             ws.dexpinfo_see{t}(:, col.vExploreInfo)=ws.dexpinfo_see{t}(:, col.vSee);
             ws.dexpinfo_nosee{t}(:, col.vExploreInfo)=ws.dexpinfo_nosee{t}(:, col.vNoSee);
             ws.dexpinfo{t}=sortrows([ws.dexpinfo_see{t}; ws.dexpinfo_nosee{t}], col.ExploreInfoPE);
             ws.d{t}(ws.dexpinfo{t}(:, col.ExploreInfoPE), col.vExploreInfo)=ws.dexpinfo{t}(:,col.vExploreInfo);
             ws.d{t}(:, col.ExploreInfoPE)=nan;
             ws.d{t}(:, col.ExploreInfoPE)=ws.d{t}(:, col.vExploreInfo)-ws.d{t}(:, col.vExplore);
             ws.d{t}(:, col.ExploreInfoPEsign)=nan;
             ws.d{t}(ws.d{t}(:, col.ExploreInfoPE)==0, col.ExploreInfoPEsign)=0;
             ws.d{t}(ws.d{t}(:, col.ExploreInfoPE)>0, col.ExploreInfoPEsign)=1;
             ws.d{t}(ws.d{t}(:, col.ExploreInfoPE)<0, col.ExploreInfoPEsign)=-1;
             ws.d{t}(:, col.ExploreInfotoOutcomePE)=nan;
             ws.d{t}(ws.expshown_trials, col.ExploreInfotoOutcomePE)= ws.d{t}(ws.expshown_trials, col.OutcomeMagnitude) - ws.d{t}(ws.expshown_trials, col.vExploreInfo);
             ws.d{t}(:, col.ExploreGambletoOutcomePE)=nan;
             ws.d{t}( ws.d{t}(:, col.Resp1)==3 & ws.d{t}(:, col.OutcomePresented)==1, col.ExploreGambletoOutcomePE)= ws.d{t}( ws.d{t}(:, col.Resp1)==3 & ws.d{t}(:, col.OutcomePresented)==1, col.OutcomeMagnitude)-ws.d{t}( ws.d{t}(:, col.Resp1)==3 & ws.d{t}(:, col.OutcomePresented)==1, col.vExplore);             
             
             
             
             
             
             % Choice entropy 
             er 
             
             
             
             
         end
         
         % Re-save 
         subjdata{os,2}.d_conflict=ws.d{1};
         subjdata{os,2}.d_control=ws.d{2};
         subjdata{os,2}.d_conflictcontrol=sortrows([ws.d{1}; ws.d{2}],  col.Onset_Offer);
         subjdata{os,2}.d{1}=subjdata{os,2}.d_conflict; % Re-save to 'd'
         subjdata{os,2}.d{2}=subjdata{os,2}.d_control;
         subjdata{os,2}.d{3}=subjdata{os,2}.d_conflictcontrol;
         
     end
end
% disp('How many explore-shown, see and no-see trials are there?'); disp(logg.n_exploreinfo);

for o1=1:1  % Hard-coded checks for new models 
% 
% 
% 
% 
% aa=load('C:\Users\eloh\Desktop\2 [Explore]\1 Brain data\p01_GV\2 First level s4Ants\p01_GV_onsets_m_v16c_vChoicePE_bpm16bpmi11.mat');
% % aa=load('/Users/EleanorL/Desktop/2 EXPLORE fMRI data/1 Brain data/p01_GV/2 First level s4Ants/p01_GV_onsets_m_v16c_vChoicePE_bpm16bpmi11.mat');
% a=aa.names'; tt=subjdata{1,2}.d_conflictcontrol;
% tcf=tt(  tt(:, col.Task)==1, :); tct=tt(  tt(:, col.Task)==2, :);
% tcf=tcf(  tcf(:, col.TrialValid)==1, :);  tct=tct(  tct(:, col.TrialValid)==1, :); 
% tcfs= tcf(tcf(:, col.OutcomePresented)==1,:);
% tcts= tct(tct(:, col.OutcomePresented)==1,:);
% tcfu= tcf(tcf(:, col.OutcomePresented)==0,:);
% tctu= tct(tct(:, col.OutcomePresented)==0,:);
% 
% 
% sum(aa.pmod(3).param{1}-tcfu(tcfu(:,col.Resp1)==1, col.vAccept))
% sum(aa.pmod(4).param{1}-tcfu(tcfu(:,col.Resp1)==2, col.vReject))
% sum(aa.pmod(5).param{1}-tcfu(tcfu(:,col.Resp1)==3, col.vExplore))
% sum(aa.pmod(6).param{1}-tctu(tctu(:,col.Resp1)==1, col.vAccept))
% sum(aa.pmod(7).param{1}-tctu(tctu(:,col.Resp1)==2, col.vReject))
% sum(aa.pmod(8).param{1}-tctu(tctu(:,col.Resp1)==3, col.vExplore))
% %
% sum(aa.pmod(11).param{1}-tcfs(tcfs(:,col.Resp1)==1, col.AcceptOutcomePE))
% sum(aa.pmod(12).param{1}-tcfs(tcfs(:,col.Resp1)==2, col.RejectOutcomePE))
% sum(aa.pmod(13).param{1}-tcfs(tcfs(:,col.Resp1)==3, col.ExploreInfotoOutcomePE))
% sum(aa.pmod(14).param{1}-tcts(tcts(:,col.Resp1)==1, col.AcceptOutcomePE))
% sum(aa.pmod(15).param{1}-tcts(tcts(:,col.Resp1)==2, col.RejectOutcomePE))
% sum(aa.pmod(16).param{1}-tcts(tcts(:,col.Resp1)==3, col.ExploreInfotoOutcomePE))
end


% Applt j warp if necessary
if logg.warpJ
    % Values etc have already been warped in the model-export script
    d_warp.jcf=behmod_pars.d_par.cf(:,   find(strcmp(behmod_pars.details.cf_moddetails{3}, 'j'))  );
    d_warp.fcf=behmod_pars.d_par.cf(:,   find(strcmp(behmod_pars.details.cf_moddetails{3}, 'f'))  );
    d_warp.ecf=behmod_pars.d_par.cf(:,   find(strcmp(behmod_pars.details.cf_moddetails{3}, 'e'))  );
    d_warp.jct=behmod_pars.d_par.ct(:,   find(strcmp(behmod_pars.details.ct_moddetails{3}, 'j'))  );
    d_warp.ect=behmod_pars.d_par.ct(:,   find(strcmp(behmod_pars.details.ct_moddetails{3}, 'e'))  );
    
    for s=1: logg.n_subjs
        ws.cf= subjdata{s,2}.d{1};
        ws.ct= subjdata{s,2}.d{2};
        
        % cF
        ws.mv=[]; 
        try ws.mv.FixedLoss= d_warp.fcf(s); end
        try ws.mv.ExploreCost= d_warp.ecf(s); end
        ws.ov=fcf_changeEnvThreat(power(ws.cf(:, col.EnvThreat), d_warp.jcf(s)),   ws.cf(:, col.NTokens),  ws.mv); 
        ws.cf(:, col.EnvThreat)=         ws.ov.EnvThreat;
        ws.cf(:, col.NTokens)=         ws.ov.NTok;
        ws.cf(:, col.pLoss)=         ws.ov.pLoss;
        ws.cf(:, col.Entropy)=         ws.ov.Entropy;
        ws.cf(:, col.EntropyNTok)=         ws.ov.EntropyNTok;
        ws.cf(:, col.VExplore)=         ws.ov.EntropyNTok;
        ws.cf(:, col.EV)=         ws.ov.EV;
        
        % ct
        ws.mv=[];
        try ws.mv.ExploreCost= d_warp.ect(s); end
        ws.ov=fct_changeEnvThreat(power(ws.ct(:, col.EnvThreat), d_warp.jct(s)),   ws.ct(:, col.NTokens),  ws.mv); 
        ws.ct(:, col.EnvThreat)=         ws.ov.EnvThreat;
        ws.ct(:, col.NTokens)=         ws.ov.NTok;
        ws.ct(:, col.pLoss)=         ws.ov.pLoss;
        ws.ct(:, col.Entropy)=         ws.ov.Entropy;
        ws.ct(:, col.EntropyNTok)=         ws.ov.EntropyNTok;
        ws.ct(:, col.VExplore)=         ws.ov.EntropyNTok;
        ws.ct(:, col.EV)=         ws.ov.EV;
        
        % Output
        subjdata{s,2}.d_conflict=ws.cf;
        subjdata{s,2}.d_control=ws.ct;
        subjdata{s,2}.d_conflictcontrol=vertcat( subjdata{s,2}.d_conflict,  subjdata{s,2}.d_control);
        subjdata{s,2}.d_conflictcontrol=sortrows(subjdata{s,2}.d_conflictcontrol,col.Onset_Offer);
        subjdata{s,2}.d{1}=subjdata{s,2}.d_conflict; % d: holds all data (1=conflict, 2=control, 3=all)
        subjdata{s,2}.d{2}=subjdata{s,2}.d_control;
        subjdata{s,2}.d{3}=subjdata{s,2}.d_conflictcontrol;
    end
    
    % Paranoid check
    paranoidcheck=0;
    if paranoidcheck
        s=1;    request.check={'EnvThreat';'NTokens';'pLoss';'Entropy';'EntropyNTok';'EV'};
        f.subplotcols=length(request.check)  ;  f.subplot_VerHorz=[0.005 0.055];f.fig_BotTop=[0.001 0.035];  f.fig_LeftRight=[0.05 0.01];   f.figwidth= 1200;  f.figheight=600;
        f.f=figure('color','w', 'Position',[50,70,f.figwidth,f.figheight]); k=1;
        %
        ws.dcf=subjdata{s,2}.d_conflict;  ws.dct=subjdata{s,2}.d_control;
        ws.et_cf = sortrows(unique( subjdata{s,2}.d_conflict(:, col.EnvThreat)));  ws.et_ct = sortrows(unique( subjdata{s,2}.d_control(:, col.EnvThreat)));
        d_check{1}= nan(size(ws.dcf,1), length(request.check));  d_check{2}= d_check{1};
        for v=1:length(request.check)
            eval(['d_check{1}(:, v)=ws.dcf(:, col.'  request.check{v} ');'])
            eval(['d_check{2}(:, v)=ws.dct(:, col.'  request.check{v} ');'])
            
            for e=1:6
                for n=1:6
                    d_check{3}{v}(7-e, n) = unique(d_check{1}( ws.dcf(:, col.EnvThreat)==  ws.et_cf(e) & ws.dcf(:, col.NTokens)==  n*2, v));
                    d_check{4}{v}(7-e, n) = unique(d_check{2}( ws.dct(:, col.EnvThreat)==  ws.et_ct(e) & ws.dct(:, col.NTokens)==  n*2, v));
                end
            end
            
            % Plot cF
            subtightplot(2,  f.subplotcols,  v, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagesc(   d_check{3}{v} ), colorbar, axis square
            title(['[cF]  ' request.check{v}])
            
            % Plot ct
            subtightplot(2,  f.subplotcols, f.subplotcols+ v, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagesc(   d_check{4}{v} ), colorbar, axis square
            title(['[ct]  ' request.check{v}])
            wv=[];
        end
        disp('Examine d_check{3}{v} and d_check{4}{v} to check that variables are okay!');
        
        input('SAVE variable to check ET warp?');
        save([where.behmodfile 'CheckETwarp_s' num2str(s) ], 'request', 'logg', 'd_check');
    end
end


% Actually computed the added value from Exploration
if strcmp(logg.onsetsmodel(1:3), 'm_v')==1 
    for s=1:logg.n_subjs
        ws.cf= subjdata{s,2}.d{1};
        ws.ct= subjdata{s,2}.d{2};
        
        % cF: VExplore = V(Explore) -  EV(EV>0)
        ws.cf(:, col.VExplore)  = ws.cf(:, col.vExplore)  -  ws.cf(:, col.EV) .* (ws.cf(:, col.EV)>0);
        
        % ct: VExplore = V(Explore) -  EV
        ws.ct(:, col.VExplore)  = ws.ct(:, col.vExplore)  -  ws.ct(:, col.EV);
        
        % Output
        subjdata{s,2}.d_conflict=ws.cf;
        subjdata{s,2}.d_control=ws.ct;
        subjdata{s,2}.d_conflictcontrol=vertcat( subjdata{s,2}.d_conflict,  subjdata{s,2}.d_control);
        subjdata{s,2}.d_conflictcontrol=sortrows(subjdata{s,2}.d_conflictcontrol,col.Onset_Offer);
        subjdata{s,2}.d{1}=subjdata{s,2}.d_conflict; % d: holds all data (1=conflict, 2=control, 3=all)
        subjdata{s,2}.d{2}=subjdata{s,2}.d_control;
        subjdata{s,2}.d{3}=subjdata{s,2}.d_conflictcontrol;
    end
end

% PLOTS
for o=1:1
    
%     disp('If model includes j parameter, SEEEE CODDEEE line 405');
    
        % NOTE: EnvThreat at this point is already warped by j term (if
        % requested), in all the variable calculated by the modelling export scirpt.
        %         . Use EnvThreatOriginal to index. BE CAREFUL if you are
        % directly using EnvThreat at any point! - explicitly check thru!
        
    
    doplot=0;
    if doplot
        
        % Collate individual datas (in plotting format)
        domeans_set= {
%                             'expGam2Info' 'ExploreInfoPE';    
%                                 'expInfo2Out' 'ExploreInfotoOutcomePE';
%                                 'expGam2Out' 'ExploreGambletoOutcomePE'
                                'softmaxvBest'  'pvBest'
                                };  % REQUEST ###
        domeans=domeans_set(:,1);
        %
        for d=1:length(domeans)
            dd{1}={cell(6,6) cell(6,6)};  % 1=cF subject trials, 2=ct subject trials, 3=cf mean, 4=ct mean, 5=cf no. trials, 6= ct no. trials
            eval(['dd{2}=col.' domeans_set{d,2} ';'])
            
%             eval(['d_' domeans{d} '=[]'])
            for s=1:logg.n_subjs
                for t=1:2
                    ws=subjdata{s,2}.d{t};
                    ws=ws(ws(:, col.OutcomePresented)==1,:); if s==1 & t==1; input('Outcome presented only!'); end
                    
                    
                    for e=1:6
                        ee=7-e;
                        for n=1:6
                            wc=ws(ws(:, col.EnvThreat)*6==e & ws(:, col.NTokens)/2==n, :);
                            dd{1}{t}{ee,n}(s)=nanmean(wc(:, dd{2}));
                        end
                    end
                end
            end
            
            
            eval(['d_' domeans{d} '(3:4)=dd{1};'])
        end
        
        % Group-means
        minsubs=round(0.5*logg.n_subjs);
        for d=1:length(domeans)
            for t=1:2
                eval(['dd=d_' domeans{d} '{t+2};'])
                for e=1:6
                    for n=1:6
                        if sum(1-isnan(dd{e,n})    )< minsubs; disp( ['[' domeans{d} ' (' num2str(t) ')  t' num2str(e) '-' num2str(n) '] Too many nans (> ' num2str(minsubs) ' subjects)']); end
                        dn{1}(e,n)=sum(isnan(dd{e,n}));
                        %
                        dn{2}(e,n)=nanmean(dd{e,n});
                    end
                end
                eval(['d_' domeans{d} '{t}=dn{2};'])
                eval(['d_' domeans{d} '{t+4}=dn{1};']) % Record no. of trials
            end
        end
        
        % Do plot
        f.plotcols=3;  f.figwidth= 1800; f.figheight=400; f.fontsize=15; f.fontsize_title=15;
        f.subplot_VerHorz=[0.1 0.05]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05]; f.taskname={'cF';'ct'};% Loads of rois
        figure('Name', 'Group means', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
        for d=1:length(domeans)
            eval(['dd=d_' domeans{d} '(1:2);'])
            subplot(length(domeans), f.plotcols, k); axis off
            title(domeans{d}); k=k+1;
            for t=1:2
                %
                subplot(length(domeans), f.plotcols, k)
                imagescnan(dd{t}, 'NanColor', [0.9 0.9 0.9]); axis square; colorbar;
                %                 title(['[' f.taskname{t} '] ' domeans{d}], 'FontSize', f.fontsize_title)
                
%                 caxis([-8 5])
                set(gca, 'FontSize',f.fontsize,'YTick',1:6, 'YTickLabel',  flipud({'1/6';'2/6';'3/6';'4/6';'5/6';'1'}),'XTick', 1:6, 'XTickLabel',2:2:12)
                k=k+1;
            end
        end
        
        
    end
    
    
end

%% (2) Format regressors
% % subjdata: Col 1= Name, Col 2=data {.d} (1=Conflict, 2=Control, 3=all), Col 3=regressor variables

input('Format regressors?   ');

for i=1:logg.n_subjs
    s=find(strcmp(subjdata, logg.subjects{i})==1);
    disp(['Subject ' num2str(i) '  (' logg.subjects{s} ')'])
    c=1;
    
    % Create variables, according to requested model
    eval(['[ws.variables c] ='  logg.onsetsmodel   '(subjdata{s,2}.d, col, c);']);
    
    subjdata{s,3}=ws.variables;
    ws=[];
end

%% (3) Save onsets ------------------------------

% Save
if request.saveonsets==1
    disp('Saving onsets to First Level folders')
    for i=1:logg.n_subjs 
        s=find(strcmp(subjdata, logg.subjects{i})==1); % Identify correct subject
        names=subjdata{s,3}.names;
        onsets=subjdata{s,3}.onsets;
        durations=subjdata{s,3}.durations;
        if isfield(subjdata{s,3},  'pmod');  pmod=subjdata{s,3}.pmod; else pmod=[]; end
        %
        if strcmp(logg.onsetsmodel(1:3), 'm_v')==1
            save([where.data_brain filesep subjdata{s,1} filesep '2 First level' logg.FLthread filesep subjdata{s,1} '_onsets_' logg.onsetsmodel '.mat'], 'names', 'onsets', 'durations', 'pmod', 'behmod_pars');
        else
            save([where.data_brain filesep subjdata{s,1} filesep '2 First level' logg.FLthread filesep subjdata{s,1} '_onsets_' logg.onsetsmodel '.mat'], 'names', 'onsets', 'durations', 'pmod');
        end
    end
end

%% (5) Post-evaluations for specific models?


% (a) Cluster models only (c6 & c7) - Empty onsets?
if sum(strcmp(logg.onsetsmodel(1:5),{'m_c5_';'m_c6_';'m_c7_';'m_c8_';'m_c9_';'m_c10_';}))==1 
    disp('Checking if onsets variables are ok (no empty events of interest) -------------------------')
    w.n_checkregs=find(strcmp(cellstr(char(subjdata{1,3}(:).names)),'n_Error'))-1;
    RegOK=cell(logg.n_subjs+1,w.n_checkregs+2); % Col 1= Subject, Col 2= Inclusion, Col 3 onwards: n events
    
    RegOK{1,1}='Subject'; RegOK{1,2}='Inclusion';
    for i=1:w.n_checkregs % Count events 
        RegOK{1, 2+i}=subjdata{1,3}.names{i};
        for s=1:logg.n_subjs
            RegOK{s+1, 2+i}=length(subjdata{s,3}.onsets{i});
        end
    end
    for s=1:logg.n_subjs; 
        RegOK{s+1,1}=logg.subjects{s}; 
        RegOK{s+1,2}=(sum([RegOK{s+1, 3:14}]>0)==12);
    end
    
    disp(RegOK)
    disp(['N ok subjects:   ' num2str(sum(cell2mat(RegOK(2:end,2))))])
    
    for o1=1:1 % OLD check  
%     RegOks=cell(logg.n_subjs+1,12+2); % Col 1= Subj, Col 2-13=No. events in onsets file, Col 14=Include subject or not?
%     RegOks{1,1}='Subject'; RegOks{1,14}='Inclusion'; for i=1:12; RegOks{1,1+i}=subjdata{1,3}.names{i}; end
%     
%     % Log nos. for each subject
%     for s=1:logg.n_subjs
%         RegOks{s+1,1}=logg.subjects{s};
%         
%         for i=1:12
%             RegOks{s+1,1+i}=length(subjdata{s,3}.onsets{i});
%         end
%         
%         % Include subject?
% 
%          if sum(cell2mat(RegOks(s+1,[2:7]))==0)<1
%             RegOks{s+1,14}=1;
%         else
%             RegOks{s+1,14}=0;
%         end
%         
%         %  Correct no. of trials in In/Out cluster?
%         RegOks{s+1,16}=RegOks{s+1,2}+RegOks{s+1,3}+RegOks{s+1,4};
%         RegOks{s+1,17}=RegOks{s+1,5}+RegOks{s+1,6}+RegOks{s+1,7};
%         RegOks{s+1,18}=RegOks{s+1,8}+RegOks{s+1,9}+RegOks{s+1,10};
%         RegOks{s+1,19}=RegOks{s+1,11}+RegOks{s+1,12}+RegOks{s+1,13};
%         RegOks{s+1,20}=length(subjdata{s,3}.onsets{find(strcmp(subjdata{s,3}.names,'n_Error'))});
%         RegOks{s+1,21}=sum(cell2mat(RegOks(s+1,16:20)));
%     end
%     
%     disp(RegOks)
%     disp(['N good subjects: ' num2str(sum(cell2mat(RegOks(2:end, 14))))])
    end
end

%% END
disp('############################################################')
disp('END');  disp(' ')
disp('Errors:');  disp(errorlog)



%% CHECK ONSETS for nans etc
% Manually clear all and load onsets + execute following 


docheck=0;
if docheck
input('Manually check onsets etc?');


for k=1:length(durations)  % Onsets and durations 
    if sum(isnan(durations{k})) + sum(isinf(durations{k}))==0 & sum(isnan(onsets{k})) + sum(isinf(onsets{k}))==0
        disp([names{k} ':  onsets and durations ok']);
    else  disp([names{k} ':  onsets and durations BAD ']);
    end
end
sum(isnan(onsets{k}))

 

for k=1:length(pmod)  % pmod contents
    if isempty(pmod(k))==1
        for p=1:length(pmod(k))
            if sum(isnan(pmod(k).param{p})) + sum(isinf(pmod(k).param{p}))==0
                disp([names{k} '    -   ' pmod(k).name '  :  pmod values ok']);
            else  disp([names{k} '    -   ' pmod(k).name '  :  pmod values BAD ']);
            end
        end
        
        
    elseif pmod(k).name
        
        
        onsets(k)
         pmod(k)
        
        onsets(:).name
        names'
        
        
        pmod(11:20).name
        
        pmod(3)
        er 
        
    end
end


end

