% Flexible beta analysis script 
clear all; close all hidden; clc, request.rois={}; request.whichbat={}; request.roi_rename=[]; 
 
% cd('/Users/EleanorL/Desktop/2 EXPLORE fMRI data/2 Second level results Ants')
cd('G:\2 [Explore]\2 Second level results s4Ants') 
% request.where_FLrois='m_v28g_ChoicePredChoice_bpji08bpji11_Basic\choice_2x3\ROI';
% request.where_FLrois='m_c14_Choice_Basic\choice_2x3\ROI';
request.where_FLrois=['m_c13_ChoiceFull_ULPEN_Basic' filesep 'choice_2x3' filesep 'ROI'];
% request.where_FLrois='m_c7_ChCluster6Full_ULPEN_Basic\choice_2x2\ROI';
% request.where_FLrois='m_v4g_ChoicevChosenAnd_bpji08bpji11_Basic\choice_2x3\ROI';
% request.where_FLrois='m_v3g_vChosenAnd_bpji08bpji11_Basic\TaskVal\ROI';
% request.where_FLrois='m_v25g_RejectOrvGamble_bpji08bpji11_Basic\TaskbyRejOr\ROI';
% request.where_FLrois='m_v23g_RejectOrvChosenAnd_bpji08bpji11_Basic\RejOr_vBU\ROI';            % vBU Rej OR
% request.where_FLrois='m_v27g_ChoiceRejectOrvChosenAnd_bpji08bpji11_Basic\RejOr_vBU\ROI';            % vBU Rej OR
% request.where_FLrois='m_v9g_vChosenAndposneg_bpji08bpji11_Basic\vBestUnchosenPosNeg\ROI';
% request.where_FLrois='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11_Basic\vBestUnchosen\ROI'; 
% request.where_FLrois='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11_Basic\choice_2x3\ROI';
% request.where_FLrois='m_v1g_vChoice_bpji08bpji11_Basic\par vExplore\ROI';
% request.where_FLrois=['m_v36g_RejectOrvMargChoDiff_bpji08bpji11_Basic' filesep 'TaskVar_2x2' filesep 'ROI'];
% request.where_FLrois='m_v37g_ChoiceRejectOrvMargChoDiff_bpji08bpji11_Basic\TaskVar_2x2\ROI';
% request.where_FLrois='m_v38g_EVGainLoss_bpji08bpji11_Basic\EVGainLoss\ROI';
% request.where_FLrois='m_v40g_vBestvWorst_bpji08bpji11_Basic\TaskVar_2x2\ROI';
% request.where_FLrois='m_v42g_pLossNTok_bpji08bpji11_Basic\par pLossNTok\ROI';
% request.where_FLrois='m_v44g_vBUposnegvMargChoDiff_bpji08bpji11_Basic\vMargCho x vBUpn\ROI';
% request.where_FLrois='m_v46g_ChoiceXvMargChoDiff_bpji08bpji11_Basic\TaskXChoicemargvcho_2x3\ROI';
% request.where_FLrois=['m_c21_RejExpFull_ULPEN_Basic' filesep 'choice_2x2' filesep 'ROI'];
% request.where_FLrois=['m_c22_ChoiceFull_HLPEN_Basic' filesep 'choice_2x3' filesep 'ROI'];
 

% WHICH rois? #################
% request.where_rois=[request.where_FLrois filesep 'Subfield fxn masked'];   % request.whichbat='anat_hpcamy';
% request.where_rois=[request.where_FLrois filesep 'Subfields binary']; request.whichbat=[];
% request.where_rois=[request.where_FLrois filesep 'Subfields ref']; request.whichbat=[];
% request.where_rois=[request.where_FLrois filesep 'Subfield x c13 HPCunmasked tc']; request.whichbat=[]; request.rois={'HPC_aL_tc_DGref';'HPC_aR_tc_DGref'; 'HPC_aL_tc_CA3ref'; 'HPC_aR_tc_CA3ref'; 'HPC_aL_tc_CA1ref'; 'HPC_aR_tc_CA1ref'; 'HPC_aL_tc_SUBref'; 'HPC_aR_tc_SUBref'}; 
% request.where_rois=[request.where_FLrois filesep 'Subfield x c13 HPCunmasked c'];  request.whichbat=[]; request.rois={'HPC_aL_c_Amyg'; 'HPC_aL_c_DGref';'HPC_aL_c_CA3ref'; 'HPC_aR_c_Amyg';  'HPC_aL_c_CA1ref'; }; 
% request.where_rois=[request.where_FLrois filesep 'Subfields 60']; request.whichbat=[];% request.where_rois=[request.where_FLrois filesep 'c13 battery'];  request.whichbat='c13battery';
% request.where_rois=[request.where_FLrois filesep 'c13 battery'];   request.whichbat='c13battery';
% request.where_rois=[request.where_FLrois filesep 'c13 amyg'];   request.whichbat=[];
% request.where_rois=[request.where_FLrois filesep 'c13 battery and amyg'];  request.whichbat='c13batteryamyg';
% request.where_rois=[request.where_FLrois filesep 'c13 HPCunmasked'];  request.whichbat=[];  request.rois={'HPC_aL_tc'; 'HPC_aR_tc' ; 'HPC_aL_c' ; 'HPC_aR_c' }; 
% request.where_rois=[request.where_FLrois filesep 'c7 rois'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep 'v28gChoice battery'];  request.whichbat='v28gbattery';
% request.where_rois=[request.where_FLrois filesep  'v3g HPC' ];   request.rois={'HPC_aL_tv001'; 'HPC_aR_tv001'; 
%     'HPC_aL_mev05fwe'; 'HPC_aR_mev05fwe';
%     };  %     'HPC_aL_vchocfMct'; 'HPC_aR_vchocfMct'; };
% request.where_rois=[request.where_FLrois filesep 'v3g hpc amyg'];   request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep  'v3g Others' ];  request.rois=  {'HPC_mR' 'PosCing' 'Striatum_R' 'rlPFC' 'vmPFC' 'Insula_L' 'Insula_R' 'Parietal_R'  'dlPFC' 'dmPFC' 'Thalamus_R' }'; 
% request.where_rois=[request.where_FLrois filesep 'v4gChoice HPC'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep  'v36g' ];  request.rois=  {'HPC_aL'; 'HPC_aR'; 'Amyg_L'; 'Amyg_R'; 'Insula_L'; 'Parietal_L'; 'Parietal_R'; 'Temporal'};  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep 'v23g'];  request.whichbat=[]; % request.rois={'roi_tc001_Amyg_L_svc';'roi_tc001_HPC_mL2'; 'roi_tc001_HPC_mL'; 'roi_tc001_IFG'}; request.roi_rename={'Amygdala_L';'HPC_mL2';'HPC_mL';'IFG'}; 
% request.where_rois=[request.where_FLrois filesep 'v36g'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep 'HPC assorted'];  request.whichbat=[]; 
% request.rois=  {'HPC_aL_stc_c13bat';'HPC_aR_stc_c13bat';'HPC_aL_sc_c13bat'; 'HPC_aR_sc_c13bat'; 'HPC_aL_v36g'; 'HPC_aR_v36g'; 'HPC_aL_v3gtv001'; 'HPC_aR_v3gtv001';'HPC_aL_apc_v28g';'HPC_aR_apc_v28g';'HPC_aL_atc_v28g'; 'HPC_aR_atc_v28g'; 'HPC_aL_s01tc_v28g';'HPC_aR_s01tc_v28g';};
% request.where_rois=[request.where_FLrois filesep 'MEC'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois ];  request.whichbat=[]; 
%  request.where_rois=[request.where_FLrois filesep 'v3g HPC SubRefMasked'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois filesep 'c22 infHPC'];  request.whichbat=[]; 
request.where_rois=[request.where_FLrois filesep 'HPC mec tc conjunc 05unc'];  request.whichbat=[]; 
 
 

% HPC rois only !! 
request.whichbat=[];
% request.rois={'rCA1_aL' 'rCA1_aR' 'rCA3_aL' 'rCA3_aR'};     % Subfield binary/ref
% request.rois={'rDG_aL' 'rCA3_aL' 'rCA1_aL' 'rSUB_aL' 'rDG_aR' 'rCA3_aR' 'rCA1_aR' 'rSUB_aR'};  % Subfields ref circuit 
% request.rois={'HPC_aL_stc' 'HPC_aR_stc' 'HPC_aL_sc' 'HPC_aR_sc'};             % c13
request.rois={ 
%     'HPC_aL';'HPC_aR'; %
%     'HPC_aL_stc';'HPC_aR_stc'; 
%     'Amyg_L'; 'Amyg_R'

'HPC_L_mectc05uncconjunc';'HPC_R_mectc05uncconjunc';
    }; 

for o=1:1
    % request.rois={'rDG_aL_60' 'rCA3_aL_60' 'rCA1_aL_60' 'rSUB_aL_60' 'rDG_aR_60' 'rCA3_aR_60' 'rCA1_aR_60' 'rSUB_aR_60'};  % Subfields ref circuit 
    % request.rois={ 'rCA1_aL_80' 'rCA1_aR_80' 'rCA3_aL_80' 'rCA3_aR_80'};   % Subfield 80
    % request.rois={ 'rCA1_aL_60' 'rCA1_aR_60' 'rCA3_aL_60' 'rCA3_aR_60'};   % Subfield 95
% request.rois={'HPC_aL_stc' 'HPC_aR_stc' 'HPC_aL_sc' 'HPC_aR_sc' 'Amyg_L' 'Amyg_R'};             % c13 amyg
% request.rois={'HPC_aL_atc' 'HPC_aR_atc' 'HPC_aL_apc' 'HPC_aR_apc'};    % v28g
% request.rois={'CA1_L_tc'; 'CA1_R_tc'; 'CA3_L_c' 'CA1_L_c'; 'CA3_L_tc';'CA3_R_tc'};  % c13-func subfields 
% request.rois={'Parietal';'Parietal_L'; 'Striatum_R_anat'};
end

% request.nonpar_correlation=0; 
request.nonpar_correlation=1; 
request.meancentre=1;

%% Set up data files 


for o1=1:1 % ROI BATTERY requests  
    
    % [c3 battery] ##################
    request.c3batt_rois={
%         'BA46'; 'BA10'; 'Striatum_L'; 'Striatum_R'; 'SupMidFG'; 'Precuneus'
%              'HPC_aL_c'; 'HPC_aR_c'; 
% 'HPC_L_sac'; 'HPC_R_sac'; 
'HPC_L_satc'; 'HPC_R_satc';        
        };
    request.c3batt_roi_rename={
%             'BA46'; 'BA10'; 'Left striatum'; 'Right striatum'; 'Superior/Mid frontal gyrus'; 'Precuneus'; 
%         'Hippocampus Left (Choice)'; 'Hippocampus Right (Choice)'; 
'Hippocampus Left (Task x Choice)'; 'Hippocampus Right (Task x Choice)';
%         'Superior/Mid frontal gyrus'; 'Precuneus'; 'HPC_aL_c'; 'HPC_aR_c';
        };
    
     % [c13] ##################
     request.c13battery_rois={'HPC_aL_stc';'HPC_aR_stc'; 'HPC_aL_sc'; 'HPC_aR_sc'; 'Striatum_R_anat'; 'BA10';'DLPFC_BA10_L';'DLPFC_BA10_R'; };
     request.c13battery_roi_rename=request.c13battery_rois;
     request.c13batteryamyg_rois= [request.c13battery_rois; {'Amyg_L';'Amyg_R'};]; 
     request.c13batteryamyg_roi_rename=request.c13batteryamyg_rois;
     
     
    % [c14 ] ##################
    request.c14battery_rois={'HPC_aL_atc' 'HPC_aR_atc' 'HPC_aL_ac' 'HPC_aR_ac' 'DLPFC_BA10_L' 'DLPFC_BA10_R'};
    request.c14battery_roi_rename= request.c14battery_rois;
%     request.c14_selectedhpc_rois={'HPC_aL_anat';'HPC_aR_anat'; 'HPC_aL_c001'; 'HPC_aR_c001'; 'HPC_aL_tc05f';'HPC_aR_tc05f'; 'HPC_aL_tc001';  'HPC_aR_tc001'; 'HPC_aL_stc_rme'; 'HPC_aR_stc_rme';};
%     request.c14_selectedhpc_roi_rename = request.c14_selectedhpc_rois; 
%     request.c14HPC_v28gmask_rois= {
%     'HPC_aL_ac';'HPC_aR_ac';   'HPC_L_c_maskedv28gtc'; 'HPC_R_c_maskedv28gtc';  'HPC_L_atc_maskedv28gpc';  'HPC_R_atc_maskedv28gpc';};
%     request.c14HPC_v28gmask_roi_rename=request.c14HPC_v28gmask_rois;
    request.c14_HPC_rois={ 'HPC_L_s001atc'; 'HPC_R_s001atc';  'HPC_L_s001ac'; 'HPC_R_s001ac'; 
     'HPC_L_s01atc'; 'HPC_R_s01atc';  'HPC_L_s01ac'; 'HPC_R_s01ac'; 
     'HPC_L_s05atc'; 'HPC_R_s05atc';
%      'HPC_aL_satc';'HPC_aR_satc';  
     }; 
    request.c14_HPC_roi_rename =request.c14_HPC_rois;
    
    % [v28g choice] ##################
    request.v28gbattery_rois={ 'HPC_aL_atc'; 'HPC_aR_atc';
    'HPC_aL_apc'; 'HPC_aR_apc'; 
    'DLPFC_L_c'; 'DLPFC_R_c'; 'Parietal_c'; 'rlPFC_L_c'};
    request.v28gbattery_roi_rename = request.v28gbattery_rois;
    request.v28gHPC_c14gmask_rois={'HPC_aL_atc_maskedc14mec_continuous' ; 'HPC_aR_atc_maskedc14mec_continuous'
    'HPC_aL_pc_maskedc14tc'; 'HPC_aR_pc_maskedc14tc'};
    request.v28gHPC_c14gmask_roi_rename=request.v28gHPC_c14gmask_rois;
    
    
    % [HPC + amygdala ANAT] ##################
    request.anat_hpcamy_rois={ 'HPC_aL'; 'HPC_aR'; 'HPC_pL'; 'HPC_pR'; 'Amygdala_L'; 'Amygdala_R'};
    request.anat_hpcamy_roi_rename = request.v28gbattery_rois;
    
    % [v25g ] ##################
    request.v25ghpc_rois = {'HPC_aL_val_atc'; 'HPC_aR_val_atc'; 'Amyg_L_val_atc'; 'Amyg_R_val_atc'; 'HPC_aL_val_tc'; 'HPC_aR_val_tc'};
    request.v25ghpc_roi_rename = request.v25ghpc_rois;
    
end
if isempty(request.whichbat )==0, eval(['request.rois=request.' request.whichbat '_rois;  request.roi_rename=request.' request.whichbat '_roi_rename;']) , end
% request.roi_rename=[];  % Empty to omit
% request.rois={};  % Empty to use all ROIs (without re-ordering)

% ################################################

for o1=1:1   % Manual settings
    log.specificsubjects=[];
    % Get file name
    request.filename=[];
    if isempty(request.filename)
        %         f=spm_select('List',request.where_rois, 'Extracted');
        %         if size(f,1)==0; error('No beta excel files!'); end
        %         if size(f,1)~=1; error('No. of beta excel files ~=1!'); end
        %         request.filename= f(1:strfind(f, '.xls')-1);
        
        cd( request.where_rois)
        f= dir('*Extracted*');
        if size(f,1)==0; error('No beta excel files!'); end
        if size(f,1)~=1; error('No. of beta excel files ~=1!'); end
        request.filename= f.name(1:strfind(f.name, '.xls')-1) ;
    end
    
    % Other non-changing setup/details
    log.roi_con_connector='-';
    
end
for o1=1:1 % Load file + auto-read details about betas (rename beh, beta, etc)
    w=pwd;  if strcmp(w(1), '/'); where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; else where.where='D:\Dropbox\SCRIPPS\1 Explore fMRI'; end;   path(pathdef); addpath(where.where)
    % Load file + select ROIs/subjects
    [n, t r]=xlsread([request.where_rois filesep request.filename '.xlsx']);
    d_betas=[t(:,1) [t(1, 2:end); num2cell(n)]]; log.orig_d_betas = d_betas; % Original table, save!
    if isempty(request.rois)==0
        w.roicon_names=d_betas(1,2:find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,1:end), 'BEH')), 1, 'first')-1); w.okconcols=[];  % ROI selection
        for r=1:length(request.rois)
            if sum(cellfun(@(x)~isempty(x),   strfind(w.roicon_names, request.rois{r}))) ==0;   disp('Column titles:' );  disp(d_betas(1,:)');  error(['Could not find requested ROI. Check name: ' request.rois{r}]); end
            
            w.okconcols=[w.okconcols   find(cellfun(@(x)~isempty(x),   strfind(w.roicon_names, [request.rois{r} log.roi_con_connector]))) ];
        end
        d_betas=[d_betas(:, [1  w.okconcols+1]) d_betas(:,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:), 'BEH')), 1, 'first') : end)];
    end
    [log.subjects log.n_subjs d_betas] = f_selectsubjects(d_betas, log.specificsubjects,  [d_betas(:,1) [{'All'}; num2cell(ones(size(d_betas,1)-1,1))]], 'All');  % Subject selection
    
    % Rename ROIs (proper name)
    if isempty(request.roi_rename)==0
        disp('Renaming of ROIs requested. Check names:'); openvar ans; ans=[request.rois request.roi_rename]; disp([request.rois request.roi_rename]); input('Continue?  ');
        w.roicon_names=d_betas(1,:);
        for r=1: length(request.rois)
            disp(['Assumed no. contrasts per roi: '  num2str(sum(cellfun(@(x)~isempty(x),  strfind(w.roicon_names,  request.rois{r}))))])
            w.cons2rename=find(cellfun(@(x)~isempty(x),  strfind(w.roicon_names,  request.rois{r})));
            for i=1:length(w.cons2rename)
                d_betas{1, w.cons2rename(i)}=[request.roi_rename{r} d_betas{1, w.cons2rename(i)}(length(request.rois{r})+1:end)];
            end
        end
        ans=[w.roicon_names' d_betas(1,:)']; input('Check line up of rename ROIs (in var window). Proceed?  ');
        request.rois=request.roi_rename;
    end
    
    % Manual renaming of behavioural variables
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'st.State')))} = 'State anxiety scores';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'st.Trait')))} = 'Trait anxiety scores';
    d_betas(1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'B.Drive')))) = {'BIS/BAS Drive scores'};
    d_betas(1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'B.Bis')))) = {'BIS scores'};
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'info.allfMRI')))} = 'Probability of following information\n from Exploring';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'o.allfMRI')))} = 'Overall money won\n on task (points)';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'info.cF_all')))} = 'Probability of following information\n from Exploring in experimental task';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'infoexnull.cF_all')))} = 'Probability of following null information\n from Exploring in experimental task';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'rt.cF_Reject')))} ='RT for Rejecting (Exp condition)';
    d_betas{1, find(cellfun(@(x) ~isempty(x), strfind(d_betas(1,:),'per.Rej_cFmct')))} ='Difference in % Rejecting (Exp>Ctrl)';
    
    
    % Read out beta details from table
    log.data_names= d_betas(1,2:end)';
    log.roicons=log.data_names(1:find(1-cellfun(@(x)isempty(x), strfind(log.data_names, 'BEH_')),1)-1);
    
    % Auto-read FL con names and ROI names
    %     log.roicons=d_betas(1,:)'; log.roicons= log.roicons(2: find(cellfun(@(x)[1- isempty(x)],  strfind(log.roicons,'BEH_')),1,'first')-1);
    wc.rc=cellfun(@(x)fliplr(x), log.roicons, 'UniformOutput',0);
    wc.firstconcharnum=strfind(wc.rc{1}, log.roi_con_connector);
    wc.firstconname=  wc.rc{1}(1:   wc.firstconcharnum(end)  );
    wc.firstconnums_perroi= find(cellfun(@(x)[1- isempty(x)],  strfind(wc.rc,wc.firstconname)));
    wc.firstcons_perroi= log.roicons(wc.firstconnums_perroi);
    wc.charsin1stroiconname=cell2mat(strfind(wc.firstcons_perroi, fliplr(wc.firstconname)))-1;
    for i=1:size(wc.firstcons_perroi,1)
        log.rois_available{i,1}=  wc.firstcons_perroi{i}( 1: wc.charsin1stroiconname(i));
    end
    disp('Available ROIs:'), disp(log.rois_available)
    log.n_cons = wc.firstconnums_perroi(2) - wc.firstconnums_perroi(1);
    wc.roi1cons= log.roicons(1:log.n_cons );
    wc.roi1_name= wc.roi1cons{1}(1:strfind(wc.roi1cons{1}, log.roi_con_connector)-1 );
    log.cons=cellfun(@(x)x(length(wc.roi1_name) +2:end),  wc.roi1cons, 'UniformOutput',0);
    log.n_rois    = sum(1-cellfun(@(x)isempty(x), strfind(log.roicons, log.cons{1})));  % If there are repeats here, the no. of ROIs will be wrong!
    log.n_cons=length(log.cons);
    log.rois=cellfun(@(x)x(  1:  length(x)-length(log.cons{1})-length(log.roi_con_connector)   ), log.roicons(find(1-cellfun(@(x)isempty(x), strfind(log.roicons, log.cons{1})))), 'UniformOutput',0);
    %     log.rois={'HPC_L_satc'; 'HPC_R_satc'}; log.n_rois=2; disp('Manual override!!')
    %     log.rois={'Hippocampus Left (Task x Choice)'; 'Hippocampus Right (Task x Choice)';};
    %
    % Read out behaviour details from table
    log.behfams_start=find(1-cellfun(@(x)isempty(x), strfind(log.data_names, 'BEH_')))+1;  %
    log.behnames=d_betas(1,log.behfams_start(1):end)';
end

% [MANUALLY specify contrasts type]  # # # # #
log.FLtype=[]; % Empty to auto select
% log.FLtype='TaskxChoice';
% log.FLtype='predChoice';
% log.FLtype='outTaskxChoice';

for o=1:1 % Identify FL con types (determines figures, stats, correlations to perform)
    
    log.con_diffs={};
    if isempty(log.FLtype)
        if sum(strcmp(log.cons, 'cF_vChosen'))+sum(strcmp(log.cons, 'ct_vChosen'))+sum(strcmp(log.cons, 'cF_vBestUnchosen'))+ sum(strcmp(log.cons, 'ct_vBestUnchosen')) ==4
            log.FLtype='vChovBU';
        elseif sum(strcmp(log.cons, 'cF_vChosen'))+sum(strcmp(log.cons, 'ct_vChosen'))+sum(strcmp(log.cons, 'cF_vBestUnchosen_pos'))+ sum(strcmp(log.cons, 'cF_vBestUnchosen_neg'))+ sum(strcmp(log.cons, 'ct_vBestUnchosen_pos')) ==5
            log.FLtype='vChovBUposneg';
        elseif sum(strcmp(log.cons, 'cF_vChosen'))+sum(strcmp(log.cons, 'ct_vChosen'))+sum(strcmp(log.cons, 'cF_Rej_vBestUnchosen'))+ sum(strcmp(log.cons, 'cF_NonRej_vBestUnchosen'))+ sum(strcmp(log.cons, 'ct_Rej_vBestUnchosen')) + sum(strcmp(log.cons, 'ct_NonRej_vBestUnchosen')) ==6
            log.FLtype='vChovBU_RejOr';
        elseif sum(strcmp(log.cons, 'cF_vChosen'))+ sum(strcmp(log.cons, 'ct_vChosen'))+ sum(strcmp(log.cons, 'cF_Rej_vBestUnchosen'))+  sum(strcmp(log.cons, 'cF_NonRej_vBestUnchosen'))+  sum(strcmp(log.cons, 'ct_vBestUnchosen')) ==5
            log.FLtype='vChovBU_cFRejOr';
        elseif sum(strcmp(log.cons, 'cF_Rej_vGamble'))+ sum(strcmp(log.cons, 'cF_NonRej_vGamble'))+ sum(strcmp(log.cons, 'ct_Rej_vGamble')) + sum(strcmp(log.cons, 'ct_NonRej_vGamble')) ==4
            log.FLtype='vGamble_RejOr';
        elseif sum(strcmp(log.cons, 'in_cF_Reject'))+sum(strcmp(log.cons,  'in_cF_Explore'))+sum(strcmp(log.cons, 'in_ct_Bomb'))+sum(strcmp(log.cons, 'in_ct_Explore'))+sum(strcmp(log.cons,  'out_cF_Accept'))+sum(strcmp(log.cons, 'out_cF_Reject'))+sum(strcmp(log.cons,  'out_ct_NoBomb'))+sum(strcmp(log.cons,  'out_ct_Bomb')) ==8
            log.FLtype='inTaskxChoice';
        elseif sum(strcmp(log.cons, 'cF_vAccept'))+ sum(strcmp(log.cons, 'cF_vExplore'))+ sum(strcmp(log.cons, 'ct_vAccept'))+ sum(strcmp(log.cons, 'ct_vReject'))+ sum(strcmp(log.cons, 'ct_vExplore'))==5
            log.FLtype='vChoice';
        elseif  sum(strcmp(log.cons, 'cF_Rej_vMargCho'))+  sum(strcmp(log.cons, 'cF_NonRej_vMargCho'))+  sum(strcmp(log.cons, 'ct_Rej_vMargCho'))+  sum(strcmp(log.cons, 'ct_NonRej_vMargCho')) ==4
            log.FLtype='vMargCho_RejOr';
        elseif  sum(strcmp(log.cons, 'cF_vBUpos_vMargChosen'))+  sum(strcmp(log.cons, 'cF_vBUneg_vMargChosen'))+  sum(strcmp(log.cons, 'ct_vBUpos_vMargChosen'))  ==3
            log.FLtype='vMargCho_vBUpn';
        elseif  sum(strcmp(log.cons, 'cF_EVGain'))+   sum(strcmp(log.cons, 'cF_EVLoss'))+   sum(strcmp(log.cons, 'ct_EVGain'))   ==3
            log.FLtype='EVGainLoss';
        elseif  sum(strcmp(log.cons, 'cF_vBest'))+   sum(strcmp(log.cons, 'cF_vWorst'))+   sum(strcmp(log.cons, 'ct_vBest'))+   sum(strcmp(log.cons, 'ct_vWorst'))==4
            log.FLtype='vBestWorst';
        elseif  sum(strcmp(log.cons, 'cF_pLossNTok'))+   sum(strcmp(log.cons, 'ct_pLossNTok')) ==2
            log.FLtype='pLossNTok';
        elseif  sum(strcmp(log.cons,  'cF_Acc_vMargChosen')) +sum(strcmp(log.cons,  'cF_Rej_vMargChosen')) +sum(strcmp(log.cons,  'cF_Exp_vMargChosen')) +sum(strcmp(log.cons,  'ct_Acc_vMargChosen')) +sum(strcmp(log.cons,  'ct_Rej_vMargChosen')) +sum(strcmp(log.cons,  'ct_Exp_vMargChosen'))
            log.FLtype='ChoicexvMargCho';
        % This should be LAST 
        elseif sum(strcmp(log.cons, 'cF_Accept')) >0 && strcmp(log.cons{find(cellfun(@(x)1-isempty(x),  strfind(log.cons, 'cF_Accept')))+5}, 'ct_Explore')
            log.FLtype='TaskxChoice';
        elseif sum(strcmp(log.cons, 'cF_Reject'))+sum(strcmp(log.cons, 'cF_Explore'))+ sum(strcmp(log.cons, 'ct_Bomb')) +sum(strcmp(log.cons, 'ct_Explore'))  ==4 && sum(strcmp(log.cons, 'cF_Accept'))+sum(strcmp(log.cons, 'ct_NoBomb'))==0
            log.FLtype='TaskxRejOr';
        else error('FL type not specified yet!');
        end
    end

    
    
    % FL Constrast differences scores (requested)
    if sum(strcmp(log.FLtype, {'TaskxChoice'}))==1
        log.con_diffs={
%             'cF_Rej-Exp'    {'cF_Reject' '-' 'cF_Explore'}   % Selected simple effects
%             'ct_Rej-Exp'    {'ct_Bomb' '-' 'ct_Explore'}
%             'Rej_cF-ct'    {'cF_Reject' '-' 'ct_Bomb'}
            'cF_Rej-Exp'    {'cF_Reject' '-' 'cF_Explore'}  % All simple effects 
            'cF_Acc-Rej'    {'cF_Accept' '-' 'cF_Reject'}
            'cF_Acc-Exp'    {'cF_Accept' '-' 'cF_Explore'}
            'ct_Rej-Exp'    {'ct_Bomb' '-' 'ct_Explore'}    %
            'ct_Acc-Rej'    {'ct_NoBomb' '-' 'ct_Bomb'}
            'ct_Acc-Exp'    {'ct_NoBomb' '-' 'ct_Explore'}
            'Acc_cF-ct'    {'cF_Accept' '-' 'ct_NoBomb'}    %
            'Rej_cF-ct'    {'cF_Reject' '-' 'ct_Bomb'}
            'Exp_cF-ct'    {'cF_Explore' '-' 'ct_Explore'}
            };
    elseif sum(strcmp(log.FLtype, {'TaskxRejOr'}))==1
        log.con_diffs={
            'cF_Rej-Exp'    {'cF_Reject' '-' 'cF_Explore'}  % All simple effects 
            'ct_Rej-Exp'    {'ct_Bomb' '-' 'ct_Explore'}    %
            'Rej_cF-ct'    {'cF_Reject' '-' 'ct_Bomb'}
            'Exp_cF-ct'    {'cF_Explore' '-' 'ct_Explore'}
            };
    elseif sum(strcmp(log.FLtype, {'inTaskxChoice'}))==1
        log.con_diffs={
            'cF_Rej-Exp'    {'in_cF_Reject' '-' 'in_cF_Explore'}
            'ct_Rej-Exp'    {'in_ct_Bomb' '-' 'in_ct_Explore'}
            'Rej_cF-ct'    {'in_cF_Reject' '-' 'in_ct_Bomb'}
            'Exp_cF-ct'    {'in_cF_Explore' '-' 'in_ct_Explore'}
            };
    elseif sum(strcmp(log.FLtype, {'outTaskxChoice'}))==1
        log.con_diffs={
            'cF_Acc-Rej'    {'out_cF_Accept' '-' 'out_cF_Reject'}
            'ct_Acc-Rej'    {'out_ct_NoBomb' '-' 'out_ct_Bomb'}
            'Acc_cF-ct'    {'out_cF_Accept' '-' 'out_ct_NoBomb'}
            'Rej_cF-ct'    {'out_cF_Reject' '-' 'out_ct_Bomb'}
            };
%     wr.whichcons_name={ 'out_cF_Accept';  'out_cF_Reject'; 'out_ct_NoBomb'; 'out_ct_Bomb' };
    elseif sum(strcmp(log.FLtype, {'predChoice'}))==1
         log.con_diffs={
            'cF_predRejMExp'    {'cF_predReject' '-' 'cF_predExplore'}
            'ct_predRejMExp'    {'ct_predReject' '-' 'ct_predExplore'}
            'predReject_cF-ct'    {'cF_predReject' '-' 'ct_predReject'}
            'predExplore_cF-ct'    {'ct_predExplore' '-' 'ct_predExplore'}
            };
    elseif sum(strcmp(log.FLtype, {'vChovBU'}))==1
        log.con_diffs={
            'cF_vCho-vBU'    {'cF_vChosen' '-' 'cF_vBestUnchosen'}
            'ct_vCho-vBU'    {'ct_vChosen' '-' 'ct_vBestUnchosen'}
            'vChosen_cF-ct'     {'cF_vChosen' '-' 'ct_vChosen'}
            'vBestUnchosen_cF-ct'     {'cF_vBestUnchosen' '-' 'ct_vBestUnchosen'}
            };
    elseif sum(strcmp(log.FLtype, {'vGamble_RejOr'}))==1
        log.con_diffs={
            'cF_Rej-NonRej'    {'cF_Rej_vGamble' '-' 'cF_NonRej_vGamble'}
            'ct_Rej-NonRej'    {'ct_Rej_vGamble' '-' 'ct_NonRej_vGamble'}
            'Rej_cF-ct'     {'cF_Rej_vGamble' '-' 'ct_Rej_vGamble'}
            'NonRej_cF-ct'     {'cF_NonRej_vGamble' '-' 'ct_NonRej_vGamble'}
            };
    elseif sum(strcmp(log.FLtype, {'vChovBU_RejOr'}))==1
        log.con_diffs={
%             'cF_vCho-RejvBU'    {'cF_vChosen' '-' 'cF_Rej_vBestUnchosen'}
%             'ct_vCho-RejvBU'    {'ct_vChosen' '-' 'ct_Rej_vBestUnchosen'}
%             'cF_vCho-NonRejvBU'    {'cF_vChosen' '-' 'cF_NonRej_vBestUnchosen'}
%             'vChosen_cF-ct'     {'cF_vChosen' '-' 'ct_vChosen'}
            'cF_vCho-RejvBU'    {'cF_vChosen' '-' 'cF_Rej_vBestUnchosen'}
            'cF_vCho-NonRejvBU'    {'cF_vChosen' '-' 'cF_NonRej_vBestUnchosen'}
            'cF_vBU_Rej-NonRej'    {'cF_Rej_vBestUnchosen' '-' 'cF_NonRej_vBestUnchosen'}
            'ct_vCho-RejvBU'    {'ct_vChosen' '-' 'ct_Rej_vBestUnchosen'}
            'ct_vCho-NonRejvBU'    {'ct_vChosen' '-' 'ct_NonRej_vBestUnchosen'}
            'ct_vBU_Rej-NonRej'    {'ct_Rej_vBestUnchosen' '-' 'ct_NonRej_vBestUnchosen'}
            'vCho_cF-ct'     {'cF_vChosen' '-' 'ct_vChosen'}
            'Rej_vBU_cF-ct'     {'cF_Rej_vBestUnchosen' '-' 'ct_Rej_vBestUnchosen'}
            'NonRej_vBU_cF-ct'     {'cF_NonRej_vBestUnchosen' '-' 'ct_NonRej_vBestUnchosen'}
            };
    elseif sum(strcmp(log.FLtype, {'vChovBU_cFRejOr'}))==1
        log.con_diffs={
            'cF_vCho-RejvBU'    {'cF_vChosen' '-' 'cF_Rej_vBestUnchosen'}
            'cF_vCho-NonRejvBU'    {'cF_vChosen' '-' 'cF_NonRej_vBestUnchosen'}
            'cF_vBU-Rej-NonRej'    {'cF_Rej_vBestUnchosen' '-' 'cF_NonRej_vBestUnchosen'}
            'vChosen_cF-ct'     {'cF_vChosen' '-' 'ct_vChosen'}
            };
    elseif sum(strcmp(log.FLtype, {'vChovBUposneg'}))==1        
        log.con_diffs={
            'cF_vCho-vBUn'    {'cF_vChosen' '-' 'cF_vBestUnchosen_neg'} % vCho and vBU
            'cF_vCho-vBUp'    {'cF_vChosen' '-' 'cF_vBestUnchosen_pos'}
            'ct_vCho-vBU'    {'ct_vChosen' '-' 'ct_vBestUnchosen_pos'}
            'cF_vBU_pos-neg'    {'cF_vBestUnchosen_pos' '-' 'cF_vBestUnchosen_neg'}
            'vCho_cF-ct'    {'cF_vChosen' '-' 'ct_vChosen'}
            'vBUp_cF-ct'    {'cF_vBestUnchosen_pos' '-' 'ct_vBestUnchosen_pos'}
%             'cF_vBU_pos-neg'    {'cF_vBestUnchosen_pos' '-' 'cF_vBestUnchosen_neg'}   % vBU only
            };
    elseif sum(strcmp(log.FLtype, {'vMargCho_RejOr'}))==1
        log.con_diffs={
            'cF_vMargCho_Rej-NonRej'   {'cF_Rej_vMargCho' '-' 'cF_NonRej_vMargCho'}
            'ct_vMargCho_Rej-NonRej'   {'ct_Rej_vMargCho' '-' 'ct_NonRej_vMargCho'}
            'Rej_vMargCho_cF-ct'   {'cF_Rej_vMargCho' '-' 'ct_Rej_vMargCho'}
            'NonRej_vMargCho_cF-ct'   {'cF_NonRej_vMargCho' '-' 'ct_NonRej_vMargCho'}
            'cF_vMargCho_NonRej-Rej'   {'cF_NonRej_vMargCho' '-' 'cF_Rej_vMargCho'}
            };
%          log.con_diffs=cell(0,2);  disp('Simple fx turned off')
    elseif sum(strcmp(log.FLtype, {'vMargCho_vBUpn'}))==1
        log.con_diffs={
            'cF_vMargCho_vBUp-vBUn'   {'cF_vBUpos_vMargChosen' '-' 'cF_vBUneg_vMargChosen'}
            'vMargCho_vBUp_cF-ct'   {'cF_vBUpos_vMargChosen' '-' 'ct_vBUpos_vMargChosen'}
            };
    elseif sum(strcmp(log.FLtype, {'EVGainLoss'}))==1
        log.con_diffs={
            'cF_EVGain-Loss'   {'cF_EVGain' '-' 'cF_EVLoss'}
            'EVGain_cF-ct'   {'cF_EVGain' '-' 'ct_EVGain'}
            };
    elseif sum(strcmp(log.FLtype, {'vBestWorst'}))==1
        log.con_diffs={
            'cF_vBest-vWorst'   {'cF_vBest' '-' 'cF_vWorst'}
            'ct_vBest-vWorst'   {'ct_vBest' '-' 'ct_vWorst'}
            'vBest_cF-ct'   {'cF_vBest' '-' 'ct_vBest'}
            'vWorst_cF-ct'   {'cF_vWorst' '-' 'ct_vWorst'}
            };
    elseif sum(strcmp(log.FLtype, {'pLossNTok'}))==1
        log.con_diffs={
            'pLossNTok_cF-ct'   {'cF_pLossNTok' '-' 'ct_pLossNTok'}
            };
    elseif sum(strcmp(log.FLtype, {'ChoicexvMargCho'}))==1
        log.con_diffs={
            'cF_Rej-Exp'    {'cF_Rej_vMargChosen' '-' 'cF_Exp_vMargChosen'}  % All simple effects 
            'cF_Acc-Rej'    {'cF_Acc_vMargChosen' '-' 'cF_Rej_vMargChosen'}
            'cF_Acc-Exp'    {'cF_Acc_vMargChosen' '-' 'cF_Exp_vMargChosen'}
            'ct_Rej-Exp'    {'ct_Rej_vMargChosen' '-' 'ct_Exp_vMargChosen'}    %
            'ct_Acc-Rej'    {'ct_Acc_vMargChosen' '-' 'ct_Rej_vMargChosen'}
            'ct_Acc-Exp'    {'ct_Acc_vMargChosen' '-' 'ct_Exp_vMargChosen'}
            'Acc_cF-ct'    {'cF_Acc_vMargChosen' '-' 'ct_Acc_vMargChosen'}    %
            'Rej_cF-ct'    {'cF_Rej_vMargChosen' '-' 'ct_Rej_vMargChosen'}
            'Exp_cF-ct'    {'cF_Exp_vMargChosen' '-' 'ct_Exp_vMargChosen'}
            };
    else log.con_diffs=cell(0,2); disp('No simple effects')
    end
    
end
for o1=1:1  % Preprocessing of betas (mean-centre, new beta scores)
    d_newconbetas=cell(log.n_subjs+1,0); log.newcons={};
    for r=1:log.n_rois
        
        
        %         % Reverse specific betas!!!!
        %         wr=cell2mat(d_betas(2:end,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:),  log.rois{r})))));
        %         wr=-1 *wr; d_betas(2:end,  find(cellfun(@(x)~isempty(x),  strfind(d_betas(1,:),  log.rois{r}))))=num2cell(wr);
        %         disp(['Reversed-valence betas for ' log.rois{r}])
        
        if request.meancentre
                wr= cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) )    );
                wr=wr-repmat(mean(wr,2), 1, log.n_cons); wr=num2cell(wr);
                d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{1}])) : find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.cons{end}])) ) = wr;
                disp('Betas mean centred!')
        end
        
        
        % Calculate difference scores
        for d=1:size(log.con_diffs,1)
            if sum(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}])) + sum(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{end}])) ~=2; error(['Error calculating difference scores. Could not find columns headed ' log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}  '  or  ' log.rois{r} log.roi_con_connector log.con_diffs{d,2}{3} ]);  end
            
            % Find data to manipulate
            wd{1}=cell2mat(  d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{1}]))) ) ;
            wd{2}= cell2mat(    d_betas(2:end,  find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector log.con_diffs{d,2}{3}])))  );
            
            % Compute difference (assume 1 operation only)
            eval(['wd{3} = wd{1} '   log.con_diffs{d,2}{2}  'wd{2};'])
            
            % Add new column to existing data
            d_betas=[d_betas [cellstr([log.rois{r} log.roi_con_connector log.con_diffs{d}] ); num2cell(wd{3})]];
            
            wd=[];
        end
        
        
        
        
        % #######################################################
        % Manually specified new scores! ###############################
        % #######################################################
        
        
        switch  log.FLtype
%             case 'predChoice';
%                 in.newconname='cF_predRejMExp';  in.cons= {'cF_predReject'; 'cF_predExplore'; 'ct_predReject'; 'ct_predExplore'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{2}  );  % Manually compute score
%                 %
%                 in.newconname='cF_predExpMRej';  in.cons= {'cF_predReject'; 'cF_predExplore'; 'ct_predReject'; 'ct_predExplore'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{2}  -   in.b{1}  );  % Manually compute score
%                 %
%                 in.newconname='predReject_cF-ct';  in.cons= {'cF_predReject'; 'cF_predExplore'; 'ct_predReject'; 'ct_predExplore'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{3}  );  % Manually compute score
%                 
            case 'TaskxChoice';
                disp('Calculating new betas for Choice models')
                in.newconname='MeanChoice'; in.cons={'cF_Accept';'cF_Reject'; 'cF_Explore'; 'ct_NoBomb';'ct_Bomb'; 'ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname]; if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{2} in.b{3} in.b{4} in.b{5} in.b{6} ],  2)     );  % Manually compute score
                %
                in.newconname='Explore'; in.cons={'cF_Explore';'ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{2}],  2)     );  % Manually compute score
                %         %
                in.newconname='RejmExp'; in.cons={'cF_Reject'; 'ct_Bomb'; 'cF_Explore'; 'ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   mean( [in.b{1} in.b{2} ],  2)        -    mean( [in.b{3} in.b{4} ],  2)      );  % Manually compute score
                %         %
                in.newconname='TaskxRejExp'; in.cons={'cF_Reject'; 'cF_Explore'; 'ct_Bomb';  'ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   (in.b{1} -  in.b{2})  - (in.b{3} -  in.b{4}   ));  % Manually compute score
              

            case 'inTaskxChoice';
                disp('Calculating new betas for in-zone Choice models')
                in.newconname='Explore'; in.cons={'in_cF_Explore';'in_ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{2}],  2)     );  % Manually compute score
%                 %         %
%                 in.newconname='RejmExp'; in.cons={'in_cF_Reject'; 'in_ct_Bomb'; 'in_cF_Explore'; 'in_ct_Explore'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(   mean( [in.b{1} in.b{2} ],  2)        -    mean( [in.b{3} in.b{4} ],  2)      );  % Manually compute score
            
             case  'vChoice';
                disp('Calculating new betas for vChoice models')
                in.newconname='vExplore'; in.cons={'cF_vExplore'; 'ct_vExplore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{2}],  2)     );  % Manually compute score
                %
                
            case  'vChovBU';
                disp('Calculating new betas for vChovBU models')
                in.newconname='TaskxValtype'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    (in.b{1}- in.b{2}) - (in.b{3}-in.b{4})    );  % Manually compute score
                %
                in.newconname='vChosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean([in.b{1} in.b{3}],2) );  % Manually compute score
                %
                in.newconname='vBestUnchosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean([in.b{2} in.b{4}],2) );  % Manually compute score
                %
            case 'vChovBUposneg';
%                 disp('Calculating new betas for vCho + vBU_posneg models')
%                 in.newconname='vChosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{4}],  2)     );  % Manually compute score
%                 %
%                 in.newconname='vBestUnchosen_pos'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(    mean( [in.b{2} in.b{5}],  2)     );  % Manually compute score
                %
                in.newconname='vBestUnchosen_posMneg'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   mean([in.b{2}  in.b{5}], 2) -  in.b{3}    );  % Manually compute score
                %
%                 in.newconname='vBUpos_cF-ct'; in.cons={'cF_vBestUnchosen_pos'; 'ct_vBestUnchosen_pos'};
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{1} - in.b{2}      );  % Manually compute score
            case 'vChovBU_RejOr';
%                 disp('Calculating new betas for vCho + vBU (RejOr) models')
%                 in.newconname='cF_vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{2}  );  % Manually compute score
%                 %
%                 in.newconname='ct_vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(     in.b{3}  -   in.b{4}  );  % Manually compute score
%                 %
%                 in.newconname='vBestUnchosen_Rej';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(  mean([in.b{1}  in.b{3}],2)   );  % Manually compute score
%                 %
%                 in.newconname='vBestUnchosen_NonRej';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(  mean([in.b{2}  in.b{4}],2) );  % Manually compute score
%                 %
%                 in.newconname='vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(  mean([in.b{1}  in.b{3}],2)  -  mean([in.b{2}  in.b{4}],2) );  % Manually compute score
%                 %
%                 in.newconname='RejvBestUnchosen_cFMct';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%                 for c=1:length(in.cons)
%                     in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
%                 end
%                 d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
%                 d_newconbetas(2:end, end)=num2cell(   in.b{1}   -  in.b{3}  );  % Manually compute score
%                 
            case 'vMargCho_RejOr'
                    %
                    in.newconname='NonRejMRej_cF-ct';  in.cons= {'cF_Rej_vMargCho';  'cF_NonRej_vMargCho';  'ct_Rej_vMargCho';  'ct_NonRej_vMargCho'};
                    for c=1:length(in.cons)
                        in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                    end
                    d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                    d_newconbetas(2:end, end)=num2cell(     in.b{2}  -   in.b{1}  - (   in.b{4}  -   in.b{3}  ));  % Manually compute score
                    %
            otherwise %,  disp('Specific simple effects may not be set up yet!')
                % v 25
                if isempty(strfind(request.where_rois, 'm_v25'))==0
                    disp('Calculating new betas for v25')
                    %
                    in.newconname='cF_vGam_RejMNonRej';  in.cons= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                    for c=1:length(in.cons)
                        in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                    end
                    d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                    d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{2}  );  % Manually compute score
                    %
                    in.newconname='ct_vGam_RejMNonRej';  in.cons= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                    for c=1:length(in.cons)
                        in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                    end
                    d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                    d_newconbetas(2:end, end)=num2cell(     in.b{3}  -   in.b{4}  );  % Manually compute score
                    %
                    in.newconname='Rej_vGam_cFMct';  in.cons= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                    for c=1:length(in.cons)
                        in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                    end
                    d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                    d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{3}  );  % Manually compute score
                    %
                    in.newconname='cF_Rej_vGam_mOthers';  in.cons= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                    for c=1:length(in.cons)
                        in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                    end
                    d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                    d_newconbetas(2:end, end)=num2cell(     in.b{1}  - mean([in.b{2}  in.b{3} in.b{4}],2)  );  % Manually compute score
                    
                end
                
        end
        
        % #######################################################
        % #######################################################
        
        wr=[];
    end
end
d_betas=[d_betas d_newconbetas];  % Append new betas
log.cons=[log.cons; log.con_diffs(:,1); log.newcons];  log.n_cons=length(log.cons);

% % Subset of subjects ?
% whichsubs= [2;8;9;11;12;13;14;15;18;20];   % LOW  'Probability of following null information\n from Exploring in experimental task'
% whichsubs=[1;3;4;5;6;7;10;16;17;19];   % HIGH 'Probability of following null information\n from Exploring in experimental task'

% whichsubs= [1;4;10;11;12;15;17;18;20]; % HIGH Trait anxiety
% whichsubs= [2;3;5;6;7;8;9;13;14;16;19]; % LOW Trait anxiety
% 
% whichsubs= [1;4;10;11;12;13;16;17;18;20] ; % HIGH state anxiety
% whichsubs= [2;3;5;6;7;8;9;14;15;19]  ;  % LOW State anxiety 
% % 
% log.subjects= log.subjects(whichsubs); log.n_subjs=length(log.subjects);  d_betas = [d_betas(1,:);  d_betas(whichsubs+1, :)];input('Some subjects only! Proceed?') ;



%% Get all stats (r_stats)
%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dostats=1;

% Get statistics (switch things on and off here!!!)
if dostats
    r_stats=[log.rois cell(log.n_rois, 5)]; k=1;  printout={}; p=1;
    for r=1:log.n_rois
        
        % Collect betas for all contrasts
        wr.betanames=cellfun(@(x)[log.rois{r}   log.roi_con_connector  x], log.cons, 'UniformOutput',0); wr.d=[];
        for i=1:length(wr.betanames)
            wr.d=[wr.d cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:),  wr.betanames{i}))))];
        end
        
        % (i) Mean + Std error + tstat + tstate df + tstat p
        r_stats{r, 2}=[ mean(wr.d)'     (std(wr.d)/sqrt(size(wr.d,1)))'    ];
        [wr.h wr.p wr.ci wr.stats]=ttest(wr.d);
        r_stats{r, 2}=[r_stats{r, 2} wr.stats.tstat' wr.stats.df' wr.p']; k=1;

        switch log.FLtype
            case 'TaskxChoice';
                wr.d=wr.d(:, [find(strcmp(log.cons, 'cF_Accept')) find(strcmp(log.cons, 'cF_Reject')) find(strcmp(log.cons, 'cF_Explore')) find(strcmp(log.cons, 'ct_NoBomb')) find(strcmp(log.cons, 'ct_Bomb')) find(strcmp(log.cons, 'ct_Explore'))]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:6) ,[2 3], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
%                 disp('ANOVA: Reject and Explore only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[2 3 5 6]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
%                 disp('ANOVA: Accept and Reject only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[1 2 4 5]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
%                 disp('ANOVA: Accept and Explore only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[1 3 4 6]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
                
                % ----------------------------------------
                r_stats{r, 4}{k,1}='Accept_cF-ct';                            % Choice, cF vs ct
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Reject_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,5));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Explore_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,3), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Rej-Exp';                            % Within task, Rej vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Rej-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,5), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Acc-Exp';                            % Within task, Acc vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Acc-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,4), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Acc-Exp';                            % Within task, Acc vs Rej
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,2));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Acc-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,4), wr.d(:,5));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
            case 'TaskxRejOr' 
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
                % ----------------------------------------
                r_stats{r, 4}{k,1}='Reject_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Explore_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Rej-Exp';                            % Within task, Rej vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,2));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Rej-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,3), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1; 
            case 'vChovBU';
                wr.TaskValtype=[wr.d(:, find(strcmp(log.cons, 'cF_vChosen')))     wr.d(:, find(strcmp(log.cons, 'ct_vChosen')))   wr.d(:, find(strcmp(log.cons, 'cF_vBestUnchosen')))    wr.d(:, find(strcmp(log.cons, 'ct_vBestUnchosen')))];
                [wr.anova]=teg_repeated_measures_ANOVA(wr.TaskValtype,[2 2], {'Task' 'Valtype'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Valtype, TxV; Col=F, df1, df2, p
                r_stats{r, 4}{k,1}='cF_vCmvBu';                            % cF      (Simple effects)
                [wr.h wr.p wr.ci wr.stats]=ttest(wr.TaskValtype(:, [1 2]));
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                r_stats{r, 4}{k,1}='ct_vCmvBu';                            % ct
                [wr.h wr.p wr.ci wr.stats]=ttest(wr.TaskValtype(:, [3 4]));
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                r_stats{r, 4}{k,1}='vChosen_cFmct';                            % vChosen
                [wr.h wr.p wr.ci wr.stats]=ttest(wr.TaskValtype(:, [1 3]));
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                r_stats{r, 4}{k,1}='vBestUnchosen_cFmct';                            % vBestUnchosen
                [wr.h wr.p wr.ci wr.stats]=ttest(wr.TaskValtype(:, [2 4]));
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                r_stats{r, 4}{k,1}='cF_vCmvBureversed';                            % cF
                [wr.h wr.p wr.ci wr.stats]=ttest(wr.TaskValtype(:, 1)   - abs(wr.TaskValtype(:, 2)) );
                [wr.h wr.p wr.ci, wr.stats]=ttest(abs(wr.TaskValtype(:, 2))   - (wr.TaskValtype(:, 1)) );
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                if wr.p<0.051;  disp([log.rois{r} '[cF]  vC > vBur:   '  num2str(wr.stats.tstat,2) '  ,  ' num2str(wr.p,3) ' *'])
                else disp([log.rois{r} '[cF]  vC > vBur:'])
                end
                r_stats{r, 4}{k,1}='ct_vCmvBureversed';                            % ct
                [wr.h wr.p wr.ci, wr.stats]=ttest(wr.TaskValtype(:, 3)   - abs(wr.TaskValtype(:, 4)));
                %             [wr.h wr.p wr.ci wr.stats]=ttest(abs(wr.TaskValtype(:, 4))   - (wr.TaskValtype(:, 3)));
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                if wr.p<0.051;  disp([log.rois{r} '[ct]  vC > vBur:   '  num2str(wr.stats.tstat,2) '  ,  '  num2str(wr.p,3) ' *'])
                else disp([log.rois{r} '[ct]  vC > vBur:'])
                end
            case 'vChovBUposneg';
                wr.cFval=[wr.d(:, find(strcmp(log.cons, 'cF_vChosen')))      wr.d(:, find(strcmp(log.cons, 'cF_vBestUnchosen_pos')))        wr.d(:, find(strcmp(log.cons, 'cF_vBestUnchosen_neg')))];
                wr.ctval=[wr.d(:, find(strcmp(log.cons, 'ct_vChosen')))      wr.d(:, find(strcmp(log.cons, 'ct_vBestUnchosen_pos')))];
                r_stats{r, 4}{k,1}='cF_vBU_negmpos';                            % cF      (Simple effects)
                [wr.h wr.p wr.ci wr.stats]=ttest(   wr.cFval(:, 3)  - wr.cFval(:, 2)  );
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                if wr.p<0.051;  disp([log.rois{r} '  cF_vBU_negmpos:   '  num2str(wr.stats.tstat,2) '  ,  '  num2str(wr.p,3) ' *'])
                else disp([log.rois{r} '  cF_vBU_negmpos:'])
                end
                r_stats{r, 4}{k,1}='vBU_posMneg';                            % cF      (Simple effects)
                [wr.h wr.p wr.ci, wr.stats]=ttest(    mean([wr.ctval(:, 2)   wr.cFval(:, 2)],2)  -  wr.cFval(:, 3)   );
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                if wr.p<0.051;  disp([log.rois{r} '  cF_vBU_negmpos:   '  num2str(wr.stats.tstat,2) '  ,  '  num2str(wr.p,3) ' *'])
                else disp([log.rois{r} '  vBU_posMneg:'])
                end
                
                
            case 'vChovBU_RejOr'
                if r==1; input('ANOVA on vBU only. Continue?'); end; log.analysis.anovatype='vBU TxC';
                wr.d=wr.d(:, [find(strcmp(log.cons, 'cF_Rej_vBestUnchosen'))  find(strcmp(log.cons, 'cF_NonRej_vBestUnchosen'))  find(strcmp(log.cons, 'ct_Rej_vBestUnchosen'))  find(strcmp(log.cons, 'ct_NonRej_vBestUnchosen'))  ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'RejOr'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
%                 if r==1; input('ANOVA on vdiff not vBU. Continue?'); end; log.analysis.anovatype='valdiff TxC';
%                 wr.d=wr.d(:, [find(strcmp(log.cons, 'cF_vCho-RejvBU'))  find(strcmp(log.cons, 'cF_vCho-NonRejvBU'))  find(strcmp(log.cons, 'ct_vCho-RejvBU'))  find(strcmp(log.cons, 'ct_vCho-NonRejvBU'))  ]);
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'RejOr'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
                
            case 'vMargCho_RejOr'
                wr.d=wr.d(:, [find(strcmp(log.cons, 'cF_Rej_vMargCho'))  find(strcmp(log.cons, 'cF_NonRej_vMargCho'))  find(strcmp(log.cons, 'ct_Rej_vMargCho'))  find(strcmp(log.cons, 'ct_NonRej_vMargCho'))  ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'RejOr'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
                
                
                
            case 'ChoicexvMargCho';
                wr.d=wr.d(:, [ find(strcmp(log.cons, 'cF_Acc_vMargChosen'))  find(strcmp(log.cons, 'cF_Rej_vMargChosen'))  find(strcmp(log.cons, 'cF_Exp_vMargChosen'))  find(strcmp(log.cons, 'ct_Acc_vMargChosen'))  find(strcmp(log.cons, 'ct_Rej_vMargChosen'))  find(strcmp(log.cons, 'ct_Exp_vMargChosen'))  ]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:6) ,[2 3], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
% %                 
%                 disp('ANOVA: Reject and Explore only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[2 3 5 6]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
% %                 
%                 disp('ANOVA: Accept and Reject only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[1 2 4 5]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
% %                 
%                 disp('ANOVA: Accept and Explore only'); if r==1; input('Continue?'); end
%                 [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,[1 3 4 6]) ,[2 2], {'Task' 'Choice'});
%                 r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
                
                % ----------------------------------------
                r_stats{r, 4}{k,1}='Accept_ct-cF';                            % Choice, cF vs ct
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Reject_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,5));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Explore_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,3), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Rej-Exp';                            % Within task, Rej vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Rej-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,5), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Acc-Exp';                            % Within task, Acc vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Acc-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,4), wr.d(:,6));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Acc-Exp';                            % Within task, Acc vs Rej
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,2));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Acc-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,4), wr.d(:,5));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
          
                
            case 'inTaskxChoice';
                wr.d=wr.d(:, [find(strcmp(log.cons, 'in_cF_Reject')) find(strcmp(log.cons, 'in_cF_Explore'))  find(strcmp(log.cons, 'in_ct_Bomb')) find(strcmp(log.cons, 'in_ct_Explore'))]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
                % ----------------------------------------
                r_stats{r, 4}{k,1}='Reject_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Explore_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Rej-Exp';                            % Within task, Rej vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,2));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Rej-Exp';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,3), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
            case 'outTaskxChoice';
                wr.d=wr.d(:, [find(strcmp(log.cons, 'out_cF_Accept')) find(strcmp(log.cons, 'out_cF_Reject')) find(strcmp(log.cons, 'out_ct_NoBomb')) find(strcmp(log.cons, 'out_ct_Bomb'))]);
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:4) ,[2 2], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
%                 
                % ----------------------------------------
                r_stats{r, 4}{k,1}='Accept_ct-cF';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,3));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='Reject_cF-ct';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,2), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='cF_Acc-Rej';                            % Within task, Rej vs Exp
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,1), wr.d(:,2));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
                r_stats{r, 4}{k,1}='ct_Acc-Rej';
                [wr.sth wr.stp wr.stci wr.ststats]=ttest(wr.d(:,3), wr.d(:,4));
                r_stats{r, 4}{k,2}=[wr.ststats.tstat wr.ststats.df wr.stp]; k=k+1;
%             otherwise disp('Simple effects stats not set up yet!');
        end
        
        %
        wr=[];
    end
    
    % Print to display ANOVA statistics
    if strcmp(log.FLtype, 'TaskxChoice')  ;  for r=1:log.n_rois  
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME Choice: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxC: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    elseif strcmp(log.FLtype, 'TaskxRejExp')  ;      for r=1:log.n_rois  
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME Choice: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxC: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    elseif strcmp(log.FLtype, 'ChoicexvMargCho')  ;          for r=1:log.n_rois  
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME Choice: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxC: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    elseif strcmp(log.FLtype, 'inTaskxChoice') + strcmp(log.FLtype, 'outTaskxChoice') ==1;     for r=1:log.n_rois   
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME Choice: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxC: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    elseif strcmp(log.FLtype, 'vMargCho_RejOr');  
        for r=1:log.n_rois   
            printout{p,1}=log.rois{r};
            
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),4) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME RejOr: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),4) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxRejOr: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),4) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
            
        end
    elseif strcmp(log.FLtype, 'vChovBU_RejOr');    for r=1:log.n_rois   
            printout{p,1}=log.rois{r};
            
            disp([r_stats{r,1} ' ---'])
            c=1; disp( ['ME Task: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Condition';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=2; disp( ['ME RejOr: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='ME Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            c=3; disp( ['TxRejOr: F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ',   p='  num2str(r_stats{r,3}(c,4),3) ])
            printout{p,2}='Condition x Choice';
            printout{p,3}=['F(' num2str(r_stats{r,3}(c,2)) ','  num2str(r_stats{r,3}(c,3)) ')=' num2str(r_stats{r,3}(c,1),3) ', p='  num2str(r_stats{r,3}(c,4),3) ];
            p=p+1;
            
            
            p=p+1;
            disp(' ')
        end
    end
end


% error('done stats')

%% Figures
%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dofig=1;
for o1=1:1  % Figure settings
        fontsize=15;
        f.fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms
        % f.fontname='Cambria';
        % f.fontname='Arial';
end
if dofig
    
    close all hidden
    f.plotcols=4;  f.figwidth= 1800; f.figheight=400; f.fontsize=30; f.fontsize_title=30;
    f.subplot_VerHorz=[0.15 0.05]; f.fig_BotTop=[0.15 0.15]; f.fig_LeftRight=[0.05 0.05];
%     f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05]; % Loads of rois
    figure('Name', 'ROI betas', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for r=1: log.n_rois
        set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname)
        
        switch log.FLtype
            case 'TaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_Accept'; 'cF_Reject';  'cF_Explore'; 'ct_NoBomb'; 'ct_Bomb'; 'ct_Explore' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 3,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 3,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Accept', 'Reject', 'Explore'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])      
            case 'TaskxRejOr'; % [Task x RejOr]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'cF_Reject';  'cF_Explore'; 'ct_Bomb'; 'ct_Explore' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject', 'Explore'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])     
            case 'ChoicexvMargCho'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_Acc_vMargChosen'; 'cF_Rej_vMargChosen';  'cF_Exp_vMargChosen'; 'ct_Acc_vMargChosen'; 'ct_Rej_vMargChosen'; 'ct_Exp_vMargChosen' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 3,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 3,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Accept', 'Reject', 'Explore'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])     
            case 'inTaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'in_cF_Reject';  'in_cF_Explore'; 'in_ct_Bomb'; 'in_ct_Explore' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('In-zone Reject', 'In-zone Explore'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])     
            case 'outTaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'out_cF_Accept';  'out_cF_Reject'; 'out_ct_NoBomb'; 'out_ct_Bomb' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('out-zone Accept', 'out-zone Reject'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])     
            case 'outTaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'out_cF_Accept';  'out_cF_Reject'; 'out_ct_Accept'; 'out_ct_Reject' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Accept', 'Reject'), end
                % %                     xticklabel_rotate
            case 'predChoice'; % [Predicted Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'cF_predReject';  'cF_predExplore'; 'ct_predReject';'ct_predExplore'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('(Predicted) Reject', '(Predicted) Explore'), end
                % %                     xticklabel_rotate
            case 'vChovBU';  % [vChosen + vBestUnchosen]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('V(Chosen)', 'V(Best Unchosen)'), end
                xlim([0.5 2.5])     
                % %                     xticklabel_rotate
            case 'vChovBUposneg'; % [vChosen + vBestUnchosen_Pos/Neg]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
%                 wr.whichcons_name={ 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vBestUnchosen_pos';}; wr.whichcons_name_label={'Exp V(BU)+';'Exp V(BU)-';'Ctrl V(BU)+';};
                wr.whichcons_name={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};  wr.whichcons_name_label={'Exp V(Cho)';'Exp V(BU)+';'Exp V(BU)-'; 'Ctrl V(Cho)';'Ctrl V(BU)'}; 
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr(r_stats{r,2}(wr.whichcons,2)', r_stats{r,2}(wr.whichcons,1)', 'y'); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:length(wr.whichcons_name), 'XTickLabel', wr.whichcons_name_label,  'FontSize',f.fontsize, 'FontName', f.fontname);  
                % set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Conflict'; 'Control'},  'FontSize',f.fontsize, 'FontName', f.fontname);  
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 if r==log.n_rois,  legend(wr.whichcons_name), end
            case 'vChovBU_RejOr';   % [RejOr vBU]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                
                % #1 Standard ------
                wr.whichcons_name={ 'cF_vChosen'; 'cF_Rej_vBestUnchosen';'cF_NonRej_vBestUnchosen' ;  'ct_vChosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen' }; wr.legend={'V(Chosen)', 'Reject V(Best Unchosen)', 'Non-Reject V(Best Unchosen)'}; 
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0)); if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 3,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 3,2)'     ); % Mean +/- 1 Std Err

%                 % #2 v Diff
%                 wr.whichcons_name={'cF_vCho-RejvBU'; 'cF_vCho-NonRejvBU';'ct_vCho-RejvBU'; 'ct_vCho-NonRejvBU'}; wr.legend={'Reject value difference';'Non-Reject value difference'}; 
%                 wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0)); if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
%                 barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
%                 
%                 % #3 Misc ------
%                 wr.whichcons_name={'vCho_cF-ct';'Rej_vBU_cF-ct'; 'NonRej_vBU_cF-ct';'cF_vBU_Rej-NonRej'};   wr.legend=wr.whichcons_name; 
%                 wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0)); if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
%                 barwitherr( r_stats{r,2}(wr.whichcons,2) ,    r_stats{r,2}(wr.whichcons,1)); % Mean +/- 1 Std Err
%                 
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois, legend(wr.legend), end
                
                % %                     xticklabel_rotate
            case 'vGamble_RejOr'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject V(Gamble)', 'Non-Reject V(Gamble)'), end
                % %                     xticklabel_rotate
            case 'vChovBU_cFRejOr'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_vChosen'; 'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen'};
                wr.whichcons_name_label= {'cF_vCho'; 'cF_Rej_vbu'; 'cF_NonRej_vbu'; 'ct_vCho'; 'ct_vbu'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2) , r_stats{r,2}(wr.whichcons,1) ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:length(wr.whichcons_name_label), 'XTickLabel', wr.whichcons_name_label,  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject V(Gamble)', 'Non-Reject V(Gamble)'), end
                % %                     xticklabel_rotate
            case 'vChoice'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_vExplore';  'ct_vExplore';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2)' ,   r_stats{r,2}(wr.whichcons,1)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 if r==log.n_rois,  legend('Reject', 'Explore'), end
                % %                     xticklabel_rotate
            case 'vMargCho_RejOr'
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_Rej_vMargCho'; 'cF_NonRej_vMargCho'; 'ct_Rej_vMargCho';'ct_NonRej_vMargCho'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject, V(Chosen)-V(Best Unchosen)', 'Non-Reject, V(Chosen)-V(Best Unchosen)'), end
                % %                     xticklabel_rotate
                xlim([0.5 2.5])     
            case 'EVGainLoss';
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_EVGain';  'cF_EVLoss';'ct_EVGain'; };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2)' ,   r_stats{r,2}(wr.whichcons,1)'   ,'y'  ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:3, 'XTickLabel', {'Exp EV+' 'Exp EV-' 'Ctrl EV+' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                %                 if r==log.n_rois,  legend('Reject', 'Explore'), end
                % %                     xticklabel_rotate
            case 'vBestWorst';
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_vBest'; 'cF_vWorst'; 'ct_vBest';'ct_vWorst'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                                ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('V(Best)', 'V(Worst'), end, xlim([0.5 2.5])                % %                     xticklabel_rotate
                ylim('auto')
            case 'pLossNTok';
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_pLossNTok';  'ct_pLossNTok'; };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2)' ,   r_stats{r,2}(wr.whichcons,1)', 'y'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
                                if r==log.n_rois,  legend('pLoss*NTokens'), end
                % %                     xticklabel_rotate
            case 'vMargCho_vBUpn';
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'cF_vBUpos_vMargChosen';'cF_vBUneg_vMargChosen';'ct_vBUpos_vMargChosen';};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2)' ,   r_stats{r,2}(wr.whichcons,1)', 'y'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:3, 'XTickLabel', {'Exp vBU+', 'Exp vBU-', 'Ctrl vBU+'},  'FontSize',f.fontsize, 'FontName', f.fontname);  % xlim([0.3 2.7])
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                                 if r==log.n_rois,  legend('Exp vBU+', 'Exp vBU-', 'Ctrl vBU+'), end
                % %                     xticklabel_rotate
            otherwise error('Which plots? '); 
                
        end
        
        
        % Signicant?
        %         [h p ci stats]=ttest(cell2mat(d_betas(2:end, 1+[129 136])));
        %         stats
        %         p
        
        
        
        
%         ylim([-5 2])
%         xlim([0.5 2.5])
        
        
        for o=1:1  % Outcome/PE Plots
            % %         % [Outcome models PE]  ###################
            %         subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
            %         wr.whichcons_name={'cF_OutcomeMagnitude';'cF_CueValue';'ct_OutcomeMagnitude';'ct_CueValue'};
            %         wr.whichcons_namelabel={'[cF] R';'[cF] Q';'[ct] R';'[ct] Q'};
            % %         wr.whichcons_namelabel={'Exp - Outcome';'Exp - EV';'Control- Outcome';'Control- EV'};
            %         wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
            %         if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
            %         barwitherr(r_stats{r,2}(wr.whichcons,2), r_stats{r,2}(wr.whichcons,1)); % Mean +/- 1 Std Err
            % %         set(gca, 'xtick', 1:length(wr.whichcons_namelabel),'xticklabel', wr.whichcons_namelabel)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
            % %         ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
            %         title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
            %         ylabel('Parameter estimates', 'FontSize',f.fontsize, 'FontName', f.fontname)
            %         set(gca,'xticklabel', wr.whichcons_namelabel, 'FontSize',f.fontsize, 'FontName', f.fontname)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
            %         xticklabel_rotate
            %
            
            
            % %         [Explore PE ]  ###################
            %         wr.whichcons_name={'cF_vExplore';'cF_Info';'ct_vExplore';'ct_Info'}; wr.whichcons_namelabel={'[cF] vExp';'[cF] InfoVal';'[ct] vExp';'[ct] InfoVal'};   % [v15c]
            %         wr.whichcons_name={'cF_InfoPE';'ct_InfoPE'}; wr.whichcons_namelabel=wr.whichcons_name; % [v15c putative]
            %         wr.whichcons_name={'cF_Info';'ct_Info'}; wr.whichcons_namelabel=wr.whichcons_name;
            %         wr.whichcons_name={'cF_GamblevExplore';'cF_ExploreInfoVal'; 'cF_ExploreOutcomeMagnitude'; 'ct_GamblevExplore'; 'ct_ExploreInfoVal';'ct_ExploreOutcomeMagnitude';};  wr.whichcons_namelabel={'cF_Gam';'cF_Info';'cF_Out'; 'ct_Gam';'ct_Info';'ct_Out'};
            %         wr.whichcons_name={'cF_ExplorePE_GamtoOutcome'; 'cF_ExplorePE_GamtoInfo'; 'cF_ExplorePE_InfotoOutcome'; 'ct_ExplorePE_GamtoOutcome'; 'ct_ExplorePE_GamtoInfo'; 'ct_ExplorePE_InfotoOutcome'; };  wr.whichcons_namelabel={'cF_GamOut';'cF_GamInfo';'cF_InfoOut'; 'ct_GamOut';'ct_GamInfo';'ct_InfoOut';};    % v18 putative
            % wr.whichcons_name={'cF_OutcomevExplore';'cF_ExploreInfoVal'; 'cF_ExploreOutcomeMagnitude'; 'ct_OutcomevExplore'; 'ct_ExploreInfoVal';'ct_ExploreOutcomeMagnitude';}; wr.whichcons_namelabel={'Gam';'Info';'Out'; 'Gam';'Info';'Out'};
            %
            %         wr.whichcons_name={ 'cF_ExploreInfoPE'; 'cF_ExploreInfotoOutcomePE'; 'cF_ExploreGambletoOutcomePE';'ct_ExploreInfoPE'; 'ct_ExploreInfotoOutcomePE';'ct_ExploreGambletoOutcomePE'; };          wr.whichcons_namelabel={'cF_G-I';'cF_I-O'; 'cF_G-O'; 'ct_G-I';'ct_I-O';'ct_G-O';};
            % % 	wr.whichcons_namelabel={'G-I';'I-O'; 'G-O'; 'G-I';'I-O';'G-O';};
            %
            
            %         %
            %         subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
            %         wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
            %         if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
            %         barwitherr(r_stats{r,2}(wr.whichcons,2), r_stats{r,2}(wr.whichcons,1)); % Mean +/- 1 Std Err
            % %         set(gca, 'xtick', 1:length(wr.whichcons_namelabel),'xticklabel', wr.whichcons_namelabel)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
            % %         ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
            %         set(gcbo,'FaceColor','y');
            %         title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
            %         ylabel('Parameter estimates', 'FontSize',f.fontsize, 'FontName', f.fontname)
            %         set(gca,'xticklabel', wr.whichcons_namelabel)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
            % %         xticklabel_rotate
            %
        end
        
        
        
                
%                 wr.whichcons_name_label ={'+ Counterfactual value';'- Counterfactual value'};
%                 set(gca,'TickDir','out',  'XTick', 1:length(wr.whichcons_name_label), 'XTickLabel', wr.whichcons_name_label,  'FontSize',10, 'FontName', f.fontname);
%                 ylim([0 1])


% wr.whichcons_name={'cF_vCho-RejvBU'; 'cF_vCho-NonRejvBU';'ct_vCho-RejvBU'; 'ct_vCho-NonRejvBU'}; wr.legend={'Reject value difference';'Non-Reject value difference'}; 
                
    end
    
    
    for o=1:1 % AD HOC
    % [print stats] ##########
    %     whichcons=[2 1 4 3];  rr=0;  % Cue & Outcome (i.e. Gamble PEs)
    
   
    % % Ad hoc one-sample t-tests
    % [h p ci stats]=ttest(cell2mat(d_betas(2:end, 1+[129 130 137 138])));
    % stats.tstat, p p
    %
    % cell2mat(d_betas(2:end, 1+[2 3 5 6]))
    % cell2mat(d_betas(2:end, 1+[8 9 11 12]))
    % openvar ans
    %
%     rr=rr+1;  disp('########################################'); disp(log.rois{rr}), ans=[log.cons(whichcons)  num2cell([r_stats{rr,2}(whichcons, 3) r_stats{rr,2}(whichcons, 5)]) ];

    
    
%     % [Legends]  ######################
% 	wr.whichcons_name_label={'Gamble to Info';'Info to Outcome'; 'Gamble to Outcome'; 'Left ventral striatum';'Right ventral striatum';'vmPFC';};
    
%     %
%     subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); axis off;
%     offset=0.15; fontsize=20;  % Legend
%     set(gca,'FontSize',35, 'FontName', f.fontname, 'LineWidth', 0.8);
%     for l=1:length(wr.whichcons_name_label)
%         text(0,1-offset*l, wr.whichcons_name_label{l},'FontSize', fontsize, 'FontName', f.fontname, 'FontWeight','bold')
%     end
%     
    end

    for o=1:1
        do=0; 
        
        
        if do 
        % Plot one contrast for ALL rois  ##############
%         wr.con='Rej_cF-ct'; 
%         wr.con='ct_Rej-Exp'; 
        wr.con='TaskxRejExp';
        
        log.cons
        wr.roicons = cellfun(@(x)[x log.roi_con_connector wr.con], log.rois, 'UniformOutput',0); 
        wr.roicon_b = cellfun(@(x)cell2mat(d_betas(2:end, strcmp(d_betas(1,:), x))), wr.roicons, 'UniformOutput',0);   
        wr.b=nan(log.n_subjs,0);  for r=1:log.n_rois, wr.b=[wr.b wr.roicon_b{r}];  end
        
        figure('Name', 'ROI betas', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
        barwitherr(std(wr.b)./sqrt(log.n_subjs),  mean(wr.b))
        set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:log.n_rois, 'XTickLabel', log.rois,  'FontSize',f.fontsize, 'FontName', f.fontname);   
        
%         ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
        title(wr.con,  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', f.fontname); set(gca,'FontSize',f.fontsize, 'FontName', f.fontname, 'LineWidth', 0.8);
%                 if r==log.n_rois,  legend('Accept', 'Reject', 'Explore'), end
                % %                     xticklabel_rotate
%                 xlim([0.5 2.5])     
        wr=[]; 
        end
        
        
        
    end
    
    
end


% error('Done with stats and fig')

%% Battery correlations (specific)
%       r_batcorr{k}: row= beh, col= roi/con

% Certain ROIs?
request.rois=log.rois;
d_meanbetas=[d_betas(1,2:end)'  num2cell(mean(cell2mat(d_betas(2:end, 2:end)))') num2cell(std(cell2mat(d_betas(2:end, 2:end)))'./sqrt(log.n_subjs))  num2cell( ttest(cell2mat(d_betas(2:end, 2:end)))' )  ]; 
[h p ci st]=ttest(cell2mat(d_betas(2:end, 2:end))); p(p>0.1)=nan; d_meanbetas(:,4)=num2cell(p'); 
openvar d_meanbetas  % mean beta - simple effects
for i=1:size(d_meanbetas,1);
    
    if isnan(d_meanbetas{i,4})==1, 
        d_meanbetas{i,4}='nsf'; 
%         d_meanbetas{i,4}=' '; 
    else d_meanbetas{i,4}=['t('  num2str(st.df(i)) ')= '  num2str(st.tstat(i),3) ', p= '  num2str(p(i),3)];
    end
    
    
    d_meanbetas{i,4}=['t('  num2str(st.df(i)) ')= '  num2str(st.tstat(i),3) ', p= '  num2str(p(i),3)];
end
titles=d_betas(1,:)';  openvar titles

% error('done with simple fx')

% (1) Behavioural inhibition battery #########################################
request.cor_inhib.beh=d_betas(1,  find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_BehInhibition')))+1: sum(find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_'))).*(find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_')))>find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_BehInhibition')))))-1)' ; 
for o=1:1 %  Adjust behaviour for correlations
    % Remove unreportable
    request.behinhib_remove={'infoexnull.ct_all'; 'infoexnull.tcF'; 'infoexnull.tct'; 'infoexnull.t2cF'; 'infoexnull.t2ct'; 'infoexnull.cF_train2fMRI'; 'infoexnull.ct_train2fMRI'; 'per_i4.cF_Reject'; 'per_i6.cF_Reject'; 'per_i4.Reject_cF-ct'; 'per_i6.Reject_cF-ct'; 'rt_i4.cF_Reject'; 'rt_i6.cF_Reject'; 'rt.Reject_cF-ct'; 'rt_i4.Reject_cF-ct'; 'rt_i6.Reject_cF-ct'; 'BIS/BAS Drive scores'; 'B.FS'; 'B.RR'; 'db.cF_areaReject'; 'db.ct_areaReject'; 'db.areaReject_cFmct'; 'db.cF_areaAccept'; 'db.ct_areaAccept'; 'db.areaAccept_cFmct'};
    for b=1:length(request.behinhib_remove)
        request.cor_inhib.beh(strcmp(request.cor_inhib.beh, request.behinhib_remove{b}))=[];
    end
    
    % New scores 
    d_betas{1, size(d_betas,2)+1}='infoexnull.cF-ct';  request.cor_inhib.beh= [request.cor_inhib.beh; d_betas{1, end}]; 
    d_betas(2:end, end) =  num2cell(  cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'infoexnull.cF')))       -     cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'infoexnull.ct')))   );
    d_betas{1, size(d_betas,2)+1}='infoexnull.all_cF-ct'; request.cor_inhib.beh= [request.cor_inhib.beh; d_betas{1, end}]; 
    d_betas(2:end, end)=  num2cell(  cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'Probability of following null information\n from Exploring in experimental task')))       -     cell2mat(d_betas(2:end, strcmp(d_betas(1,:), 'infoexnull.ct_all')))   );
    
    
%     
    
    for o1=1:1 % New beh
        k=size(d_betas,2)+1;
        
        d_betas{1,k}='per.ct_Reject';
        d_betas(2:end, k)= num2cell( 2*cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.Reject'),1,'first'))) -  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Reject'),1,'first')))     );
        k=k+1;
        
        d_betas{1,k}='per.cF_Accept';
        d_betas(2:end, k)= num2cell( 1- (cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Reject'),1,'first'))) +  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Explore'),1,'first'))))     );
        k=k+1;
        d_betas{1,k}='per.ct_Accept';
        d_betas(2:end, k)= num2cell( 1- (cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Reject'),1,'first'))) +  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Explore'),1,'first'))))     );
        k=k+1;
        
        d_betas{1,k}='per.Accept_ct-cF';
        d_betas(2:end, k)= num2cell( cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.ct_Accept'),1,'first'))) -  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), 'per.cF_Accept'),1,'first')))     );
        k=k+1;
    end

    request.cor_inhib.beh=[request.cor_inhib.beh;  {'per.Accept_ct-cF'; 'per.cF_Accept'}]; 
    
    
    
    
end


% cell2mat(d_betas(2:end,146)) - cell2mat(d_betas(2:end,145))


request.cor_inhib.rois=request.rois;
switch log.FLtype
    case 'TaskxChoice';  
        
        request.cor_inhib.con={
%             'cF_Reject'; 
            'cF_Rej-Exp';    
%         'Rej_cF-ct'; 'TaskxRejExp'
        };
    case 'TaskxRejExplore';   request.cor_inhib.con={ 'cF_Reject';  'cF_Rej-Exp';     'Rej_cF-ct'; 'TaskxRejExp' };
    case 'inTaskxChoice';   request.cor_inhib.con={'in_cF_Reject'; 'cF_Rej-Exp';'Rej_cF-ct'};
    case 'outTaskxChoice';   request.cor_inhib.con={'in_cF_Reject'; 'Rej_cF-ct'};
    case 'predChoice';   request.cor_inhib.con={'cF_predReject'; 'cF_predExplore';'predReject_cF-ct'; 'cF_predRejMExp'};
    case 'vChovBUposneg';   request.cor_inhib.con= { 'cF_vChosen'; 'vCho_cF-ct'; 
            'cF_vBestUnchosen_pos';  'cF_vBestUnchosen_neg'; 'vBestUnchosen_posMneg';  'cF_vBU_pos-neg'; 'vBUp_cF-ct'; 
            'cF_vCho-vBUp'; 'cF_vCho-vBUn'}; 
    case 'vChovBU'; request.cor_inhib.con =  {'cF_vChosen'; 'cF_vBestUnchosen'; 'cF_vCho-vBU'; 'vChosen_cF-ct'; 'vBestUnchosen_cF-ct'; 'TaskxValtype'};
    case 'vChovBU_RejOr';  request.cor_inhib.con= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'cF_vCho-RejvBU'; 'cF_vBU_Rej-NonRej'; 'Rej_vBU_cF-ct'; 'NonRej_vBU_cF-ct'; 'vCho_cF-ct' };
    case 'vGamble_RejOr'; request.cor_inhib.con= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'cF_vGam_RejMNonRej'; 'Rej_vGam_cFMct'; 'cF_Rej_vGam_mOthers'};
    case 'vChovBU_cFRejOr'; request.cor_inhib.con =  {'cF_vCho-RejvBU'; 'cF_vCho-NonRejvBU'; 'cF_vBU-Rej-NonRej'; 'vChosen_cF-ct'};
    case 'vMargCho_RejOr' ; 
        request.cor_inhib.con =  {
%             'cF_Rej_vMargCho' ; 'Rej_vMargCho_cF-ct';
%             'cF_NonRej_vMargCho';   
            'NonRej_vMargCho_cF-ct'
%             'cF_vMargCho_Rej-NonRej'; 'NonRejMRej_cF-ct'
            };
    case 'EVGainLoss';  request.cor_inhib.con =  {'cF_EVGain';'cF_EVLoss';'cF_EVGain-Loss';'EVGain_cF-ct'}; 
    case 'vBestWorst'; request.cor_inhib.con =  {'cF_vBest-vWorst';'ct_vBest-vWorst';'vBest_cF-ct'; 'vWorst_cF-ct' };
    case 'pLossNTok'; request.cor_inhib.con =  {'cF_pLossNTok';'pLossNTok_cF-ct'};
	case 'ChoicexvMargCho'; request.cor_inhib.con =  {'cF_Acc_vMargChosen'; 'cF_Rej_vMargChosen';  'cF_Exp_vMargChosen';  'cF_Rej-Exp';'cF_Acc-Rej'; 'cF_Acc-Exp'; 'Acc_cF-ct' ; 'Rej_cF-ct'; 'Exp_cF-ct' };
    case 'vMargCho_vBUpn'; request.cor_inhib.con =  {'cF_vBUpos_vMargChosen';'cF_vBUneg_vMargChosen' ;'cF_vMargCho_vBUp-vBUn'   ; 'vMargCho_vBUp_cF-ct'};
    otherwise, request.cor_inhib.con =  {};
end



% request.cor_inhib.con= {'cF_Rej_vGamble' 'cF_vGam_RejMNonRej';  'Rej_vGam_cFMct'; 'cF_Rej_vGam_mOthers'};  % v25
% request.cor_inhib.con={    % see log.cons for list  ############
% % 'cF_vBestUnchosen_negMpos'
% 'cF_vBestUnchosen_pos';  'cF_vBestUnchosen_neg';
% 'vBestUnchosen_posMneg'; 'vBUpos_cF-ct';
%  'RejvBestUnchosen_cFMct'
% };
rc=1; request.cor_inhib.roicons={}; % Compile roicons + execute battery
for r=1: length(request.cor_inhib.rois)    
    for c=1:length(request.cor_inhib.con)
        request.cor_inhib.roicons{rc,1}=[request.cor_inhib.rois{r}  log.roi_con_connector    request.cor_inhib.con{c}];  rc=rc+1;
    end
end
k=1; r_batcorr{k,1}=[[{' ' } request.cor_inhib.roicons'];  [request.cor_inhib.beh  cell(length(request.cor_inhib.beh),  length(request.cor_inhib.roicons))]];
for b=2:size(r_batcorr{k,1},1), for rc=2:size(r_batcorr{k,1},2)
        wbc.beh=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{b,1}), 1, 'first')));
        wbc.roiconbeta= cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{1,rc}))));
        
        if request.nonpar_correlation
            [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Spearman');  disp('Spearman correlation'); %
%             [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  disp('Kendalls correlation'); %
%             wbc.stat=wbc.r;   % Print r or p?5
%             wbc.stat=wbc.p;
            wbc.stat=['ta=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
        else
            [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta);  %
            %         [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  disp('Kendalls correlation'); %
%             wbc.stat=wbc.r;   % Print r or p?5
%             wbc.stat=wbc.p;
            wbc.stat=['r=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
        end
        
        if wbc.p<0.001;  r_batcorr{k,1}{b,rc}= ['*** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
        elseif wbc.p<0.01;  r_batcorr{k,1}{b,rc}= ['** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
        elseif wbc.p<=0.05;  r_batcorr{k,1}{b,rc}= ['* ' wbc.stat ];  %[num2str(wbc.stat,3) ' *'];
        elseif wbc.p<0.1; r_batcorr{k,1}{b,rc}=wbc.stat ;
%         else r_batcorr{k,1}{b,rc}=wbc.stat ;
        end
end,end
openvar r_batcorr{1,1}  % Put probabilities in 2nd col if necessary

% error('Done with BehInhib correlation battery')






% % AD HOC PLOTS ------------------------------
% roicon= 'HPC_aL_tc-cF_Rej-Exp';
roicon='HPC_aL_tv001-NonRej_vMargCho_cF-ct'
beh='Trait anxiety scores';
% beh='State anxiety scores';
% beh='Probability of following null information\n from Exploring in experimental task';
% beh='BIS scores';
% beh='infoexnull.all';
% beh='per.cF_Accept';
% beh='infoexnull.ct'; 
% beh='RT for Rejecting (Exp condition)';


% r={'HPC_aR_sc-cF_EVGain-Loss'}
% 
% % % --- 
% % % roicon=r{1}; 
% wr.rc= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),roicon)) ));  wr.b= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),beh)) ));
% figure('color', 'w'), scatter(wr.rc, wr.b), lsline; xlabel(roicon,'FontSize',20), ylabel(sprintf(beh),'FontSize',20)
% if      request.nonpar_correlation==0, [r p]= corr(wr.rc, wr.b); title(['r='  num2str(r) ', p=' num2str(p)],'FontSize',20)
% else      [r p]= corr(wr.rc, wr.b,'type', 'Kendall'); title(['tau='  num2str(r) ', p=' num2str(p)],'FontSize',20)
% end; set(gca,'FontSize',20)
% 
% % % -------------------------------------

% 

% corrwhat={'rCA1_L_tc-cF_Rej-Exp'; 'rCA3_L_tc-cF_Rej-Exp';'Trait anxiety scores' }; % c13 masked by subfield
% corrwhat={'HPC_aL_stc-cF_Rej-Exp';  'HPC_aL_sc-cF_Rej-Exp'; 'Trait anxiety scores'};
% corrwhat={'rCA1_aL-cF_Rej-Exp'; 'rCA3_aL-cF_Rej-Exp'; 'State anxiety scores'};
% corrwhat={'rCA1_aL_60-cF_Rej-Exp'; 'rCA3_aL_60-cF_Rej-Exp'; 'State anxiety scores'};
% corrwhat={'rCA1_aL-cF_Rej-Exp'; 'rDG_aL-cF_Rej-Exp'; 'Trait anxiety scores'};
% 
% 
% d_corrwhat= cellfun(@(x)cell2mat(d_betas(2:end, strcmp(d_betas(1,:),x))), corrwhat, 'UniformOutput',0);
% d_corrwhat= [d_corrwhat{:}];
% [r p]=corr(d_corrwhat);  
% % [r p]=corr(d_corrwhat,'type', 'Kendall');  
% 
% disp([[{'r'}; corrwhat] [corrwhat'; num2cell(r)]]); disp(' ');disp(' ')
% disp([[{'p'}; corrwhat] [corrwhat'; num2cell(p)]]); disp(' ');disp(' ')

% 
% 
% 
% 
% 
% cons={'cF_Rej_vMargCho'; 'cF_NonRej_vMargCho'; 'ct_Rej_vMargCho'; 'ct_NonRej_vMargCho'};
% cons={'cF_vChosen';'cF_vBestUnchosen'; 'ct_vChosen';'ct_vBestUnchosen'}; 
% for c=1:length(cons)
%     wc.c1= cell2mat(d_betas(2:end, strcmp(d_betas(1,:), ['rCA1_aL-' cons{c}])));
%     wc.c2= cell2mat(d_betas(2:end, strcmp(d_betas(1,:), ['rCA3_aL-' cons{c}])));
%     
%     [h p]= ttest(wc.c1, wc.c2);
%     disp([cons{c} '   : p=' num2str(p)])   
% end
% 


% (2) Exploration battery  #########################################
do_explorebat=1;
if do_explorebat
    request.cor_explore.beh=d_betas(1,  find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_Explore')))+1:end)';
    request.cor_explore.con=log.cons;
    request.cor_explore.rois=request.rois;
    switch log.FLtype
%         case 'TaskxChoice';   request.cor_explore.con={'cF_Explore';'ct_Explore'; 'Explore'; 'RejmExp'};
        case 'TaskxChoice';   request.cor_explore.con={'cF_Explore';'ct_Explore'; 'Explore'; 'RejmExp'};
        case 'vChovBU';   request.cor_explore.con={'cF_vChosen';'cF_vBestUnchosen'; 'ct_vChosen';'ct_vBestUnchosen'; 'vBestUnchosen'; 'vChosen'; 'cF_vCho-vBU'; 'vChosen_cF-ct'; 'vBestUnchosen_cF-ct'; 'TaskxValtype'};
        case 'vChovBUposneg';
        case 'vChovBU_RejOr';
        case 'vChoice'; request.cor_explore.con={'cF_vExplore';'ct_vExplore';'vExplore'};
    end
    % request.cor_explore.con={
    %     'vChosen';  'vBestUnchosen_pos';   'cF_vBestUnchosen_neg';
    % 'vBestUnchosen_posMneg';   'vBestUnchosen_negMpos';
    % 'cF_vExplore';'ct_vExplore'
    % 'cF_CueValue';'cF_OutcomeMagnitude';'ct_CueValue';'ct_OutcomeMagnitude'
    %     };
    
    rc=1; request.cor_explore.roicons={}; % Compile roicons + execute battery
    for r=1: length(request.cor_explore.rois)
        for c=1:length(request.cor_explore.con)
            request.cor_explore.roicons{rc,1}=[request.cor_explore.rois{r}  log.roi_con_connector    request.cor_explore.con{c}];  rc=rc+1;
        end
    end
    k=2;  r_batcorr{k,1}=[ [{' ' } request.cor_explore.roicons'];   [request.cor_explore.beh  cell(length(request.cor_explore.beh),  length(request.cor_explore.roicons))] ];
    for b=2:size(r_batcorr{k,1},1); for rc=2:size(r_batcorr{k,1},2)
            wbc.beh=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{b,1}), 1, 'first')));
            wbc.roiconbeta= cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{1,rc}))));
            
            if request.nonpar_correlation
                [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Spearman');  disp('Spearman correlation'); %
                %             wbc.stat=wbc.r;   % Print r or p?5
                %             wbc.stat=wbc.p;
                wbc.stat=['ta=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
            else
                [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta);  %
                %         [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  disp('Kendalls correlation'); %
                %             wbc.stat=wbc.r;   % Print r or p?5
                %             wbc.stat=wbc.p;
                wbc.stat=['r=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
            end
            
            if wbc.p<0.001;  r_batcorr{k,1}{b,rc}= ['*** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
            elseif wbc.p<0.01;  r_batcorr{k,1}{b,rc}= ['** ' wbc.stat  ];;  % [num2str(wbc.stat,3) ' **'];
            elseif wbc.p<=0.05;  r_batcorr{k,1}{b,rc}= ['* ' wbc.stat ];  %[num2str(wbc.stat,3) ' *'];
                %         elseif wbc.p<0.1; r_batcorr{k,1}{b,rc}=wbc.stat ;
                %         else r_batcorr{k,1}{b,rc}=wbc.stat ;
            end
            %
            %         [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta);  %
            %         [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  %
            %         wbc.stat=wbc.r;   % Print r or p?5
            % %         wbc.stat=wbc.p;
            %         wbc.stat=['r=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];
            %
            %         if wbc.p<0.01;  r_batcorr{k,1}{b,rc}= [wbc.stat  ' **'];;  % [num2str(wbc.stat,3) ' **'];
            %         elseif wbc.p<=0.05;  r_batcorr{k,1}{b,rc}= [wbc.stat ' *'];  %[num2str(wbc.stat,3) ' *'];
            % %         elseif wbc.p<0.1; r_batcorr{k,1}{b,rc}=wbc.stat ;
            % %         else r_batcorr{k,1}{b,rc}=wbc.stat ;
            %         end
        end; end
    openvar r_batcorr{2,1}  % Put probabilities in 2nd col if necessary
    error('Done with Explore correlation battery')
end


for o1=1:1  % Request correlations (not recorded)
%     wr.corr={'per.cF_Reject';  'Trait anxiety scores'};
%     [wr.r wr.p]=corr(cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), wr.corr{1})))), cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), wr.corr{2})))));
%     %
%     wr.partialcorr={'Hippocampus Left (Task x Choice)-cF_Rej-Exp';  'per.cF_Reject';  'Trait anxiety scores'};
%     wr.partialcorr={'Hippocampus Left (Task x Choice)-cF_Rej-Exp';  'per.cF_Reject';  'State anxiety scores'};
%     [wr.r wr.p]=partialcorr([cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), wr.partialcorr{1})))) cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), wr.partialcorr{2}))))  cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), wr.partialcorr{3}))))]   );
%     disp(wr.r);  disp(wr.p)
end

%% Adhoc random plots/analysis etc

% % Plot imagesc for coactivation
% input('Select rois for imagesc-ing coactivation?'); log.rois=log.rois(1:end-2); log.n_rois=length(log.rois);
% log.rois_names={'BA46';'BA10';'Right striatum';'Superior MFG';'Precuneus'};
% close all hidden; figure; set(gcf,'color','w'); a=cell2mat(  r_batcorr{2,1}(2:end, 2:end));
% imagescnan(a, 'NanColor', [0.9 0.9 0.9]); axis square; colorbar;
% title('Parameter estimates for Explore choices', 'FontSize', 25)
% set(gca, 'FontSize', 15, 'YTick',1:log.n_rois, 'YTickLabel',log.rois_names,'XTick',1:log.n_rois, 'XTickLabel',log.rois_names)
% 
% % Correlations?
% d_exp=d_betas(:,82:end);
% corrcols=[51 81];
% [r p]=corr(cell2mat(d_betas(2:end, corrcols(1))), cell2mat(d_betas(2:end, corrcols(2)))); r,p

%% Targetted/requested correlations + plots (r_corr)

r_corr=cell(0, 5); k=1;% (1) Beta name, (2) Beh name, (3) corr r, (4) corr p, (5) data, (6) corr sig?

request.plotcorr={'o1'}; 



for o1=1:1 % [Choice & c13 rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o1', request.plotcorr))==1
        
%         d_betas= d_betas([1 2 4:9 11:end],:);    disp('Excluding 2 subjects who show opposite pattern of betas!'); input('continue?'); 
        % LEFT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_stc'; wc.roiconname='Parameter estimates from Left Inferior\n Hippocampus (Ap/Av, Reject > Explore)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
        
        
        
%         
%         %
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State'; wc.behname='State anxiety scores'; 
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         % RIGHT HPC, cF Rej-Exp ##############################
%         wc.roiname='HPC_aR_stc'; wc.roiconname='Parameter estimates from Right Inferior\n Hippocampus(Exp, Reject > Explore)';
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.State';   wc.behname='State anxiety scores'; 
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         
%         
%         % Superior HPC rois 
%         
%         % LEFT HPC, cF Rej-Exp ##############################
%         wc.roiname='HPC_aL_sc'; wc.roiconname='Parameter estimates from Left Superior\n Hippocampus(Exp, Reject > Explore)';
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State'; wc.behname='State anxiety scores'; 
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         % RIGHT HPC, cF Rej-Exp ##############################
%         wc.roiname='HPC_aR_sc'; wc.roiconname='Parameter estimates from Right Superior\n Hippocampus(Exp, Reject > Explore)';
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.behname= 'st.State';   wc.behname='State anxiety scores'; 
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         
%     
    
    end
end
for o2=1:1 % [Choice & Subfield ROIs: HPC TxC Beh inhibition]
    
    if sum(strcmp('o2', request.plotcorr))==1
        sfx='_60';
        sfx=[]; 
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';
        wc.roiname=['rCA1_aL' sfx];  wc.betaname=[wc.roiname '-cF_Rej-Exp'];   wc.roiconname='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA1_aR' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA3_aL' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA3_aR' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; wc=[];
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.roiname=['rCA1_aL' sfx];  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA1_aR' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA3_aL' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname=['rCA3_aR' sfx]; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; wc=[];
      
        
    end
end
for o3=1:1 % [v28g FL & v28g rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o3', request.plotcorr))==1
        
        % Inferior HPC, cF Rej-Exp ##############################
%         wc.roiname='HPC_aL_atc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.roiname='HPC_aL_s01tc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
        wc.roiconname='Parameter estimates from Left Inferior \n Hippocampus (Exp, Reject > Explore)'; 
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        % % % ----------------
%         wc.roiname='HPC_aR_atc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.roiname='HPC_aR_s01tc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
        wc.roiconname='Parameter estimates from Right Inferior \n Hippocampus (Exp, Reject > Explore)'; 
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
%         
%         % Superior HPC, cF Rej-Exp ##############################
%         wc.roiname='HPC_aL_apc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Left Superior \n Hippocampus (Exp, Reject > Explore)'; 
%         %
%         wc.behname= 'st.State';   wc.behname='State anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.behname= 'B.Bis'; wc.behname='BIS scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         % % % --------------------------------------
%         wc.roiname='HPC_aR_apc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Right Superior \n Hippocampus (Exp, Reject > Explore)'; 
%         %
%         wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.behname= 'B.Bis'; wc.behname='BIS scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         
%         
    
    
    end
end
for o4=1:1 % [Choice & c13-Subfield ROIs: HPC TxC Beh inhibition]
    
    if sum(strcmp('o4', request.plotcorr))==1
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';
        wc.roiname='CA1_L_tc';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];   wc.roiconname='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname='CA1_R_tc';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname='CA3_L_c';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; wc=[];       
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.State';   wc.behname='State anxiety scores'; 
        wc.roiname='CA1_L_tc';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];   wc.roiconname='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname='CA1_R_tc';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Right \n CA1 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        %
        wc.roiname='CA3_L_c';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  wc.roiconname='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1; 
        
        
    end
end
for o5=1:1 % [v36g FL & v36g rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o5', request.plotcorr))==1
        
        % Non-Reject value tracking ##############################
        % Personality correlations  -----
        wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value difference tracking on Non-Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-NonRej_vMargCho_cF-ct']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r, wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        % Exploration  -----
        wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
        wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'm.cF_i';   wc.behname='BIS scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        
%         % Reject trials, value tracking ##############################
%         % Personality correlations  -----
%         wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
%         wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.roiname='HPC_aR'; wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
%         wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
        
        %   ===
        wc.roiname='HPC_aL';  wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct'];
        wc.behname= 'b.BIS';   wc.behname='BIS scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.roiname='HPC_aR';  wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
%         wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct'];
%         wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
        
        
        
        % Behaviour -----
        wc.roiname='HPC_aR'; wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct']; wc.behname='RT for Rejecting (Exp condition)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname='Difference in % Rejecting (Exp>Ctrl)';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
    
    end
end
for o6=1:1 % [v36g FL & v36g rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o6', request.plotcorr))==1
        
        % Non-Reject value tracking ##############################
        % Personality correlations  -----
        wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value difference tracking on Non-Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-NonRej_vMargCho_cF-ct']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r, wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        % Exploration  -----
        wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
        wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'm.cF_i';   wc.behname='BIS scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        
%         % Reject trials, value tracking ##############################
%         % Personality correlations  -----
%         wc.roiname='HPC_aL'; wc.roiconname='Parameter estimates from Left Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
%         wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.roiname='HPC_aR'; wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp condition)';
%         wc.betaname=[wc.roiname '-cF_Rej_vMargCho']; wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
        
        %   ===
        wc.roiname='HPC_aL';  wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct'];
        wc.behname= 'b.BIS';   wc.behname='BIS scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
%         %
%         wc.roiname='HPC_aR';  wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
%         wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct'];
%         wc.behname= 'b.BIS';   wc.behname='BIS scores';  
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}=wc.roiconname; k=k+1;
%         %
        
        
        
        % Behaviour -----
        wc.roiname='HPC_aR'; wc.roiconname='Parameter estimates from Right Anterior Hippocampus\n(Value tracking on Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-Rej_vMargCho_cF-ct']; wc.behname='RT for Rejecting (Exp condition)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname='Difference in % Rejecting (Exp>Ctrl)';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
    
    end
end
for o7=1:1 % [v23g FL & v23g rois: HPC amyg Beh inhibition]
    
    if sum(strcmp('o7', request.plotcorr))==1
        % Summary: The direction of these Correlations doesnt really make
        % sense ot me right now. 
        
        % Non-Reject value tracking ##############################
        % Personality correlations  -----
        wc.roiname='roi_tc001_HPC_mL2'; wc.roiconname='Parameter estimates from Left Mid Hippocampus\n(V(Best Unchosen) tracking on Non-Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-NonRej_vBU_cF-ct']; wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
        % Choice correlations  -----
        wc.roiname='roi_tc001_HPC_mL2'; wc.roiconname='Parameter estimates from Left Mid Hippocampus\n(V(Best Unchosen) tracking on Non-Reject trials, Exp>Ctrl)';
        wc.betaname=[wc.roiname '-NonRej_vBU_cF-ct']; wc.behname= 'per.cF_NotAccept';   
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.roiname='roi_tc001_HPC_mL2'; wc.roiconname='Parameter estimates from Left Mid Hippocampus\n(V(Best Unchosen) tracking on Non-Reject trials, Exp condition)';
        wc.betaname=[wc.roiname '-cF_NonRej_vBestUnchosen']; wc.behname= 'per.cF_NotAccept';   
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.roiname='roi_tc001_Amyg_L_svc'; wc.roiconname='Parameter estimates from Left Amygdala\n(V(Best Unchosen) tracking on Non-Reject trials, Exp condition)';
        wc.betaname=[wc.roiname '-cF_NonRej_vBestUnchosen']; wc.behname= 'per.cF_NotAccept';   
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.roiname='roi_tc001_Amyg_L_svc'; wc.roiconname='Parameter estimates from Left Amygdala\n(V(Best Unchosen) tracking on Non-Reject trials, Exp condition)';
        wc.betaname=[wc.roiname '-cF_NonRej_vBestUnchosen']; wc.behname= 'per.cF_Reject';   
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.nonpar_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
    end
end

  
% r_corr=r_corr(4:6,:);
% close all hidden
% Automated plot all requested correlations
disp('Targeted/requested correlations for plotting:'); 
f.scatter_dotsize=100;   f.scatter_linewidth=4;   f.FontSize=25; f.fontname='PT Sans Caption';   r_corr{:,1}
f.plotcols=1;  f.figwidth= f.plotcols*400; f.figheight=ceil(size(r_corr,1)/ (f.plotcols))*400;
f.subplot_VerHorz=[0.2 0.1]; f.fig_BotTop=[0.25 0.1]; f.fig_LeftRight=[0.1 0.05]; k=1;
figure('Name', 'Beh - beta ', 'NumberTitle', 'off', 'Position', [100 250 f.figwidth f.figheight], 'Color', 'w');  
for c=1:size(r_corr,1)
    disp(r_corr{k,1}) 
    
    % Plot rank 
    input('Requested: Plot rank rather than score. OK?  ');
    wc.os= r_corr{c,5};  
    r_corr{c,5}= sortrows(r_corr{c,5},1 ); 
    r_corr{c,5}(:,1)=1:log.n_subjs;
    r_corr{c,5}= sortrows(r_corr{c,5},2); 
    r_corr{c,5}(:,2)=1:log.n_subjs;
%     r_corr{c,1}=['Ranked: '  lower(r_corr{c,1}(1) ) r_corr{c,1}(2:end)];
%     r_corr{c,2}=['Ranked: '  lower(r_corr{c,2}(1) ) r_corr{c,2}(2:end)] ;
    r_corr{c,1}=['Ranked: '   r_corr{c,1}];
    r_corr{c,2}=['Ranked: '  r_corr{c,2}];
         
    
    subtightplot(ceil(size(r_corr,1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
%     scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2)); lsline
    scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2), f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth)
     
    
    ylabel(r_corr{c,2}, 'FontSize',f.FontSize, 'FontName', f.fontname);  
%     ylabel(sprintf(r_corr{c,2}), 'FontSize',f.FontSize, 'FontName', f.fontname)
%     xlabel(sprintf(['Parameter estimates from\n' r_corr{c,1}]), 'FontSize',15, 'FontName', f.fontname)
    xlabel(sprintf(r_corr{c,1}), 'FontSize',f.FontSize, 'FontName', f.fontname)
    set(gca,'FontSize',f.FontSize, 'FontName', f.fontname, 'LineWidth', 0.8);
    if r_corr{c,6}==1; 
        if request.nonpar_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', f.fontname)    
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', f.fontname)
        end
    else
        if request.nonpar_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', f.fontname) 
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', f.fontname)
        end
    end
    
    
% %     Extra adhoc stuff
%     xlim([-1 2])
xlim('auto')
    
    % 
    k=k+1;
end


%% 

clc, close all hidden
request.comp_rois={
%     'HPC_aL_stc_c13bat'     'HPC_aL_sc_c13bat'  % c13 compare
%     'HPC_aR_stc_c13bat'     'HPC_aR_sc_c13bat'
%     'HPC_aL_stc_c13bat'     'HPC_aL_v36g'       % c13 stc vs v36g
%     'HPC_aR_stc_c13bat'     'HPC_aR_v36g'
    'HPC_aL_stc_c13bat'     'HPC_aL_v3gtv001'       % c13 stc vs v3g
    'HPC_aR_stc_c13bat'     'HPC_aR_v3gtv001'
%     'rCA1_aL' 'rCA3_aL'     % Subfield binary
%     'rCA1_aR' 'rCA3_aR'

% 'HPC_aL_atc_v28g'       'HPC_aL_apc_v28g';    % v28g 
% 'HPC_aR_atc_v28g'       'HPC_aR_apc_v28g';
% 'HPC_aL_s01tc_v28g'       'HPC_aL_apc_v28g';
% 'HPC_aR_s01tc_v28g'       'HPC_aR_apc_v28g';
% 'HPC_aL_s01tc_v28g'       'HPC_aL_v36g'   % v28g vs v36g
% 'HPC_aR_s01tc_v28g'       'HPC_aR_v36g' 
% 'HPC_aL_atc_v28g'       'HPC_aL_v36g' 
% 'HPC_aR_atc_v28g'       'HPC_aR_v36g' 
    };
for o=1:1 % Other attempted ROIs (not useful)
    
    

% 'HPC_aL_atc_v28g'       'HPC_aL_apc_v28g';    % v28g 
% 'HPC_aR_atc_v28g'       'HPC_aR_apc_v28g';
% 'HPC_aL_s01tc_v28g'       'HPC_aL_apc_v28g';
% 'HPC_aR_s01tc_v28g'       'HPC_aR_apc_v28g';
% 'HPC_aL_s01tc_v28g'       'HPC_aL_v36g'   % v28g vs v36g
% 'HPC_aR_s01tc_v28g'       'HPC_aR_v36g' 
% 'HPC_aL_atc_v28g'       'HPC_aL_v36g' 
% 'HPC_aR_atc_v28g'       'HPC_aR_v36g' 

end
request.comp_FLcons={'cF_Rej_vMargCho'; 'cF_NonRej_vMargCho'; 'ct_Rej_vMargCho'; 'ct_NonRej_vMargCho' };

% Compile !
% d_cb: format=request.comp_rois. Order of cons in each cell is request.comp_FLcons
d_cb= cell(size(request.comp_rois));
for rr=1:size(request.comp_rois,1)
    for r=1:2
        d_cb{rr,r}=[];
        for c=1:length(request.comp_FLcons)
            wc.name=[request.comp_rois{rr,r} log.roi_con_connector request.comp_FLcons{c}];
            if sum(strcmp(d_betas(1,:),wc.name))~=1; error('no matching roicons ~=1'); end
            d_cb{rr,r}=[d_cb{rr,r} cell2mat(d_betas(2:end, strcmp(d_betas(1,:),wc.name)))];
        end
    end
end

% Mean across ROIs
k=3; comb=[1 2];
d_cb{k,1}= (d_cb{comb(1),1}  + d_cb{comb(2),1})./2;    request.comp_rois{k,1} =[request.comp_rois{comb(1),1} ' + '  request.comp_rois{comb(2),1}];
d_cb{k,2}= (d_cb{comb(1),2}  + d_cb{comb(2),2})./2;      request.comp_rois{k,2} = [request.comp_rois{comb(1),2} ' + '  request.comp_rois{comb(2),2}];

% Everything from here is manual !!
disp('Manual between-ROI analysis!!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
r_crstats= request.comp_rois;   k=1; 
for rr=1:size(request.comp_rois,1)  % Roi x Task x Choice 
    
%     rr=2
    
    wr.d1=[d_cb{rr,1} d_cb{rr,2}];  openvar wr.d1
    
    % MANUAL !! 
    if rr==1;  disp(' '); disp('Comparison: Roi x Task x RejOr -------------- '); disp(' '); end
    [wr.anova]=teg_repeated_measures_ANOVA(wr.d1 ,[2 2 2], {'ROI' 'Task' 'RejOr'});  % Col: F, df1, df2, p
    r_crstats{rr,k}= [wr.anova.labels'  num2cell(wr.anova.R)];
    disp([request.comp_rois{rr,1} '     vs     ' request.comp_rois{rr,2}   '     -------'])
    for s=1:size(r_crstats{rr,k},1)
        wr.p= [r_crstats{rr,k}{s,1} '  FX:  '];
        if r_crstats{rr,k}{s,5}>0.1;  wr.p1=' nsf';
        else  wr.p1= ['F('  num2str(r_crstats{rr,k}{s,3}) ','  num2str(r_crstats{rr,k}{s,4}) ') = ' num2str(r_crstats{rr,k}{s,2},3) '    p= '  num2str(r_crstats{rr,k}{s,5},3)];
        end
        wr.p(30:29+length(wr.p1))=wr.p1; 
        disp(wr.p)
    end
    disp(' '); k=k+1; 
end


for rr=1:size(request.comp_rois,1)  % Roi x Task
    wr.d2=[d_cb{rr,1}(:, [3 4]) d_cb{rr,2}(:, [3 4])];
     
    % MANUAL !! 
    if rr==1;  disp(' '); disp('Comparison: Roi x Task (Non-Rej only)  -------------- '); disp(' '); end
    [wr.anova]=teg_repeated_measures_ANOVA(wr.d2 ,[2 2], {'ROI' 'Task'});  % Col: F, df1, df2, p
    r_crstats{rr,k}= [wr.anova.labels'  num2cell(wr.anova.R)];
    disp([request.comp_rois{rr,1} '     vs     ' request.comp_rois{rr,2}   '     -------'])
    for s=1:size(r_crstats{rr,k},1)
        wr.p= [r_crstats{rr,k}{s,1} '  FX:  '];
        if r_crstats{rr,k}{s,5}>0.1;  wr.p1=' nsf';
        else  wr.p1= ['F('  num2str(r_crstats{rr,k}{s,3}) ','  num2str(r_crstats{rr,k}{s,4}) ') = ' num2str(r_crstats{rr,k}{s,2},3) '    p= '  num2str(r_crstats{rr,k}{s,5},3)];
        end
        wr.p(30:29+length(wr.p1))=wr.p1; 
        disp(wr.p)
    end
    disp(' '); k=k+1; 
    
    
end


for rr=1:size(request.comp_rois,1)  % Roi x Choice (Task effect)
    wr.d3=[ d_cb{rr,1}(:, 1) - d_cb{rr,1}(:, 3)     d_cb{rr,1}(:, 2) - d_cb{rr,1}(:, 4)  d_cb{rr,2}(:, 1) - d_cb{rr,2}(:, 3)  d_cb{rr,2}(:, 2) - d_cb{rr,2}(:, 4)      ]; 
    % 1-3= Rej, cF-ct.    2-4=Non Rej, cF-ct.
     
    % MANUAL !! 
    if rr==1;  disp(' '); disp('Comparison: Roi x Choice (Task fx)  -------------- '); disp(' '); end
    [wr.anova]=teg_repeated_measures_ANOVA(wr.d3 ,[2 2], {'ROI' 'Choice'});  % Col: F, df1, df2, p
    r_crstats{rr,k}= [wr.anova.labels'  num2cell(wr.anova.R)];
    disp([request.comp_rois{rr,1} '     vs     ' request.comp_rois{rr,2}   '     -------'])
    for s=1:size(r_crstats{rr,k},1)
        wr.p= [r_crstats{rr,k}{s,1} '  FX:  '];
        if r_crstats{rr,k}{s,5}>0.1;  wr.p1=' nsf';
        else  wr.p1= ['F('  num2str(r_crstats{rr,k}{s,3}) ','  num2str(r_crstats{rr,k}{s,4}) ') = ' num2str(r_crstats{rr,k}{s,2},3) '    p= '  num2str(r_crstats{rr,k}{s,5},3)];
        end
        wr.p(30:29+length(wr.p1))=wr.p1; 
        disp(wr.p)
    end
    disp(' '); k=k+1; 
    
    
end
for rr=1:size(request.comp_rois,1)  % Roi x Task (Choice effect)
    wr.d4=[ d_cb{rr,1}(:, 1) - d_cb{rr,1}(:, 2)     d_cb{rr,1}(:, 3) - d_cb{rr,1}(:, 4)  d_cb{rr,2}(:, 1) - d_cb{rr,2}(:, 2)  d_cb{rr,2}(:, 3) - d_cb{rr,2}(:, 4)      ]; 
    % 1-2= cF, Rej-NonRej    2-4=ct,  Rej-NonRej
     
    % MANUAL !! 
    if rr==1;  disp(' '); disp('Comparison: Roi x Task (Choice fx)  -------------- '); disp(' '); end
    [wr.anova]=teg_repeated_measures_ANOVA(wr.d4,[2 2], {'ROI' 'Task'});  % Col: F, df1, df2, p
    r_crstats{rr,k}= [wr.anova.labels'  num2cell(wr.anova.R)];
    disp([request.comp_rois{rr,1} '     vs     ' request.comp_rois{rr,2}   '     -------'])
    for s=1:size(r_crstats{rr,k},1)
        wr.p= [r_crstats{rr,k}{s,1} '  FX:  '];
        if r_crstats{rr,k}{s,5}>0.1;  wr.p1=' nsf';
        else  wr.p1= ['F('  num2str(r_crstats{rr,k}{s,3}) ','  num2str(r_crstats{rr,k}{s,4}) ') = ' num2str(r_crstats{rr,k}{s,2},3) '    p= '  num2str(r_crstats{rr,k}{s,5},3)];
        end
        wr.p(30:29+length(wr.p1))=wr.p1; 
        disp(wr.p)
    end
    disp(' '); k=k+1; 
    
    
end

% Hemisphere as a factor
disp(' '); disp('Comparison: Hemi x Roi x Task x Choice -------------- ');disp(' '); 
rr=1; wr.dl=[d_cb{rr,1} d_cb{rr,2}];
rr=2; wr.dr=[d_cb{rr,1} d_cb{rr,2}];
wr.d= [wr.dl wr.dr];
[wr.anova]=teg_repeated_measures_ANOVA(wr.d ,[2 2 2 2], {'Hemi' 'ROI' 'Task' 'Choice'});  % Col: F, df1, df2, p
rr=4; r_crstats{rr,k}= [wr.anova.labels'  num2cell(wr.anova.R)];
%     disp([request.comp_rois{rr,1} '     vs     ' request.comp_rois{rr,2}   '     -------'])
for s=1:size(r_crstats{rr,k},1)
    wr.p= [r_crstats{rr,k}{s,1} '  FX:  '];
    if r_crstats{rr,k}{s,5}>0.1;  wr.p1=' nsf';
    else  wr.p1= ['F('  num2str(r_crstats{rr,k}{s,3}) ','  num2str(r_crstats{rr,k}{s,4}) ') = ' num2str(r_crstats{rr,k}{s,2},3) '    p= '  num2str(r_crstats{rr,k}{s,5},3)];
    end
    wr.p(30:29+length(wr.p1))=wr.p1;
    disp(wr.p)
end
disp(' '); k=k+1;



rr=1; 
[d_cb{rr,1} d_cb{rr,2}]
openvar ans




for o=1:1 % References plots: Pos/Neg tracking of Pos/Neg values = WHAT ?
    x=-10:1:10; close all hidden
    f.linewidth=5;    f.FontSize=15; f.fontname='PT Sans Caption'; 
    figure('color','w'), k=1; 
    
     
    % Pos value, pos tracking 
    y=2*x+0;  subplot(2,2,k); k=k+1; 
    plot(x,y, 'LineWidth',f.linewidth), xlabel('Value', 'FontSize', f.FontSize, 'FontName', f.fontname), ylabel('Neural activation', 'FontSize', f.FontSize, 'FontName', f.fontname)
    set(gca, 'FontSize', f.FontSize, 'FontName', f.fontname)
    title('Positive values, Positive parameter estimates','FontSize', f.FontSize, 'FontName', f.fontname)
    xlim([0 10]), ylim([0 20])
    
    
    
    % Pos value, Neg tracking 
    y=-2.*x+0;  subplot(2,2,k); k=k+1; 
    plot(x,y, 'LineWidth',f.linewidth), xlabel('Value', 'FontSize', f.FontSize, 'FontName', f.fontname), ylabel('Neural activation', 'FontSize', f.FontSize, 'FontName', f.fontname)
    set(gca, 'FontSize', f.FontSize, 'FontName', f.fontname)
    title('Positive values, Negative parameter estimates','FontSize', f.FontSize, 'FontName', f.fontname)
    xlim([0 10]), ylim([-20 0])
    
    
    % Neg value, pos tracking 
    y=2*x+0;  subplot(2,2,k); k=k+1; 
    plot(x,y, 'LineWidth',f.linewidth), xlabel('Value', 'FontSize', f.FontSize, 'FontName', f.fontname), ylabel('Neural activation', 'FontSize', f.FontSize, 'FontName', f.fontname)
    set(gca, 'FontSize', f.FontSize, 'FontName', f.fontname)
    title('Positive values, Positive parameter estimates','FontSize', f.FontSize, 'FontName', f.fontname)
    xlim([-10  0]), ylim([-20 0])
    
    
    % Neg value, neg tracking 
    y=-2*x+0;  subplot(2,2,k); k=k+1; 
    plot(x,y, 'LineWidth',f.linewidth), xlabel('Value', 'FontSize', f.FontSize, 'FontName', f.fontname), ylabel('Neural activation', 'FontSize', f.FontSize, 'FontName', f.fontname)
    set(gca, 'FontSize', f.FontSize, 'FontName', f.fontname)
    title('Negative values, Negative parameter estimates','FontSize', f.FontSize, 'FontName', f.fontname)
    xlim([-10  0]), ylim([0 20])
    
end 


