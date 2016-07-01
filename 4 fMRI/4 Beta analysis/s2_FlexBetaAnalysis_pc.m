% Flexible beta analysis script 
clear all; close all hidden; clc, request.rois={}; request.whichbat={}; request.roi_rename=[]; 

cd('G:\2 [Explore]\2 Second level results s4Ants')
% request.where_FLrois='m_v28g_ChoicePredChoice_bpji08bpji11_Basic\choice_2x3\ROI';
% request.where_FLrois='m_v28g_ChoicePredChoice_bpji08bpji11_Basic\choice_2x2\ROI';
% request.where_FLrois='m_c14_Choice_Basic\choice_2x3\ROI';
% request.where_FLrois='m_c13_ChoiceFull_ULPEN_Basic\choice_2x3\ROI';
% request.where_FLrois='m_c17_ChCluster6_Basic\choice_2x2\ROI';
% request.where_FLrois='m_v4g_ChoicevChosenAnd_bpji08bpji11_Basic\choice_2x3\ROI';
request.where_FLrois='m_v3g_vChosenAnd_bpji08bpji11_Basic\TaskVal\ROI';
% request.where_FLrois='m_v25g_RejectOrvGamble_bpji08bpji11_Basic\TaskbyRejOr\ROI';
% request.where_FLrois='m_v23g_RejectOrvChosenAnd_bpji08bpji11_Basic\RejOr_vBU\ROI';
% request.where_FLrois='m_v9g_vChosenAndposneg_bpji08bpji11_Basic\vBestUnchosen\ROI';
% request.where_FLrois='m_v30g_cFRejectOrvChosenAnd_bpji08bpji11_Basic\vBestUnchosen\ROI'; 
% request.where_FLrois='m_c20g_ChoicePredChoice_ULPEN_bpji08bpji11_Basic\choice_2x3\ROI';



% WHICH rois? #################
% request.where_rois=[request.where_FLrois filesep 'Anat HPC Amyg'];   request.whichbat='anat_hpcamy';
% request.where_rois=[request.where_FLrois filesep 'Subfields binary']; request.whichbat=[];
% request.where_rois=[request.where_FLrois filesep 'Subfields 80']; request.whichbat=[];% request.where_rois=[request.where_FLrois filesep 'c13 battery'];  request.whichbat='c13battery';
% request.where_rois=[request.where_FLrois filesep 'c13 battery'];  request.whichbat='c13battery';
% request.where_rois=[request.where_FLrois filesep 'v28gChoice battery'];  request.whichbat='v28gbattery';
% request.where_rois=[request.where_FLrois filesep  'v3g HPC' ];    request.rois={'HPC_aL_tv001';'HPC_aR_tv001';'HPC_aL_mev05fwe'; 'HPC_aR_mev05fwe'}; request.whichbat={};
request.where_rois=[request.where_FLrois filesep  'v3g Others' ];  % request.roi_rename=  {'HPC_mR' 'Insula_L' 'Insula_R' 'Parietal_R' 'PosCing' 'Striatum_R' 'Thalamus_R' 'dlPFC' 'dmPFC' 'rlPFC' 'vmPFC'}'; 
% request.where_rois=[request.where_FLrois filesep 'v4gChoice HPC'];  request.whichbat=[]; 
% request.where_rois=[request.where_FLrois ];  request.whichbat=[]; 




% request.where_rois=[request.where_FLrois filesep '1 All\ME Choice']; request.whichbat=[];
% request.where_rois=[request.where_FLrois filesep '1 All\TxC']; request.whichbat=[];



% v28g Choice MEC: request.where_rois=[request.where_FLrois filesep '2 HPC\Mec or similar']; request.whichbat=[];
% c13 % request.where_rois=[request.where_FLrois filesep 'HPC']; request.whichbat=[]; request.rois={'roi_tcc001_HPC_aL';'roi_tcc001_HPC_aR'; 'roi_mec001_HPC_aL'; 'roi_mec001_HPC_aR'};
% c13 ME Choice: request.where_rois=[request.where_FLrois filesep 'MEC'];  request.whichbat=[];  request.rois={ 'HPC_aL_sc' 'HPC_aR_sc' 'roi_mec001_DLPFC_BA10_L' 'roi_mec001_Parietal' 'roi_mec001_Striatum_R'}';

% request.where_rois=[request.where_FLrois filesep 'MEC']; request.whichbat=[];

% request.rois={ 'rCA1_aL_80' 'rCA1_aR_80' 'rCA3_aL_80' 'rCA3_aR_80'};
    
% 
% request.rois={
%     'rCA1_aL'
%     'rCA1_aR'
%     'rCA3_aL'
%     'rCA3_aR'};

% 
% request.rois={
%     'DLPFC_L_c'
%     'DLPFC_R_c'
%     'Parietal_c'
%     'rlPFC_L_c'};

for o=1:1 % MAC ##
    
% cd('/Users/EleanorL/Dropbox/WorkPC/Conflict SL')
% % request.where_FLrois='m_v4g_ChoicevChosenAnd_bpji08bpji11_Basic';
% request.where_FLrois='m_c14_Choice_Basic';
% % request.where_FLrois='m_c13_ChoiceFull_ULPEN_Basic';
% request.where_FLrois='m_v28g_ChoicePredChoice_bpji08bpji11_Basic';

% % cd('/Users/EleanorL/Dropbox/WorkPC/Conflict beta')
% % request.where_rois='m_c14_Choice_Basic/c14 battery';   request.whichbat='c14battery';
% % request.where_rois='m_v28g_ChoicePredChoice_bpji08bpji11_Basic';     request.whichbat='v28gbattery';
% % request.where_rois='c13g/c13g battery';    request.whichbat='c13gbattery';
end


 
request.kendalls_correlation=1; 
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
    
     % [c13g] ##################
     request.c13gbattery_rois={ 'HPC_L_stc';'HPC_R_stc'; 'HPC_L_sc';'HPC_R_sc';'DLPFC_L';'DLPFC_R';'Parietal';'OFC_L';'Striatum_R';};
     request.c13gbattery_roi_rename=[];
     request.c13battery_rois={'HPC_aL_stc';'HPC_aR_stc'; 'HPC_aL_sc'; 'HPC_aR_sc'; 'Striatum_R'; 'BA10';'DLPFC_BA10_L';'DLPFC_BA10_R'; };
     request.c13battery_roi_rename=request.c13battery_rois;
     
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
    request.anat_hpcamy_rois={ 'HPC_aL'; 'HPC_aR'; 'Amygdala_L'; 'Amygdala_R'};
    request.anat_hpcamy_roi_rename = request.v28gbattery_rois;
    
    % [v25g ] ##################
    request.v25ghpc_rois = {'HPC_aL_val_atc'; 'HPC_aR_val_atc'; 'Amyg_L_val_atc'; 'Amyg_R_val_atc'; 'HPC_aL_val_tc'; 'HPC_aR_val_tc'};
    request.v25ghpc_roi_rename = request.v25ghpc_rois;
    
end
% request.whichbat='anat_hpcamy';
if isempty(request.whichbat )==0, eval(['request.rois=request.' request.whichbat '_rois;  request.roi_rename=request.' request.whichbat '_roi_rename;']) , end
% request.rois={'HPC_L_sc';'HPC_R_sc' ;'HPC_L_stc';'HPC_R_stc' };
% request.rois={'HPC_aL_c';'HPC_aL_stc';'HPC_aL_tc';'HPC_aR_c';'HPC_aR_satc';'HPC_aR_stc';'HPC_aR_tc';};  % c14 HPC
request.roi_rename=[];  % Empty to omit
% request.rois={};  % Empty to use all ROIs (without re-ordering)
% % % 



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
    w=pwd;  if strcmp(w(1), '/'); where.where='/Users/EleanorL/Dropbox/SCRIPPS/5 Explore fMRI'; else where.where='D:\Dropbox\SCRIPPS\5 Explore fMRI'; end; addpath(where.where)
    
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
% log.FLtype='predChoice';   % v28g
% log.FLtype='vChovBU_RejOr';

for o=1:1 % Identify FL con types (determines figures, stats, correlations to perform)
    
    log.con_diffs={};
    if isempty(log.FLtype)
        if sum(strcmp(log.cons, 'cF_Accept')) >0 && strcmp(log.cons{find(cellfun(@(x)1-isempty(x),  strfind(log.cons, 'cF_Accept')))+5}, 'ct_Explore')
            log.FLtype='TaskxChoice';
        elseif sum(strcmp(log.cons, 'cF_vChosen'))+sum(strcmp(log.cons, 'ct_vChosen'))+sum(strcmp(log.cons, 'cF_vBestUnchosen'))+ sum(strcmp(log.cons, 'ct_vBestUnchosen')) ==4
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
        else error('FL type not specified yet!');
        end
    end

    % FL Constrast differences scores (requested)
    if sum(strcmp(log.FLtype, {'TaskxChoice'}))==1
        log.con_diffs={
            'cF_Rej-Exp'    {'cF_Reject' '-' 'cF_Explore'}
            'ct_Rej-Exp'    {'ct_Bomb' '-' 'ct_Explore'}
            'Rej_cF-ct'    {'cF_Reject' '-' 'ct_Bomb'}
            };
    elseif sum(strcmp(log.FLtype, {'inTaskxChoice'}))==1
        log.con_diffs={
            'cF_Rej-Exp'    {'in_cF_Reject' '-' 'in_cF_Explore'}
            'ct_Rej-Exp'    {'in_ct_Bomb' '-' 'in_ct_Explore'}
            'Rej_cF-ct'    {'in_cF_Reject' '-' 'in_ct_Bomb'}
            };
    
        
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
            'cF_vCho-RejvBU'    {'cF_vChosen' '-' 'cF_Rej_vBestUnchosen'}
            'ct_vCho-RejvBU'    {'ct_vChosen' '-' 'ct_Rej_vBestUnchosen'}
            'cF_vCho-NonRejvBU'    {'cF_vChosen' '-' 'cF_NonRej_vBestUnchosen'}
            'vChosen_cF-ct'     {'cF_vChosen' '-' 'ct_vChosen'}
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
            'cF_vCho-vBUn'    {'cF_vChosen' '-' 'cF_vBestUnchosen_neg'}
            'cF_vCho-vBUp'    {'cF_vChosen' '-' 'cF_vBestUnchosen_pos'}
            'ct_vCho-vBU'    {'ct_vChosen' '-' 'ct_vBestUnchosen_pos'}
            'vCho_cF-ct'    {'cF_vChosen' '-' 'ct_vChosen'}
            };
    end
    
end

for o1=1:1  % Preprocessing of betas (mean-centre, new beta scores)
    d_newconbetas=cell(log.n_subjs+1,0);
    log.newcons={};
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
                
            case 'inTaskxChoice';
                disp('Calculating new betas for in-zone Choice models')
                in.newconname='Explore'; in.cons={'in_cF_Explore';'in_ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{2}],  2)     );  % Manually compute score
                %         %
                in.newconname='RejmExp'; in.cons={'in_cF_Reject'; 'in_ct_Bomb'; 'in_cF_Explore'; 'in_ct_Explore'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   mean( [in.b{1} in.b{2} ],  2)        -    mean( [in.b{3} in.b{4} ],  2)      );  % Manually compute score
            case  'vChovBU';
                disp('Calculating new betas for vChovBU models')
                in.newconname='vChosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{3}],  2)     );  % Manually compute score
                %
                in.newconname='vBestUnchosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{2} in.b{4}],  2)     );  % Manually compute score
                %
                in.newconname='vBestUnchosen_cFmct'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   in.b{2}  - in.b{4} );  % Manually compute score
                %
                in.newconname='vChosen_cFmct'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   in.b{1}  - in.b{3} );  % Manually compute score
                % 
                in.newconname='vCho-vBU_cFmct'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell((in.b{1}  - in.b{2} ) -  (in.b{3}  - in.b{4}) );  % Manually compute score
                
            case 'vChovBUposneg';
                disp('Calculating new betas for vCho + vBU_posneg models')
                in.newconname='vChosen'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{1} in.b{4}],  2)     );  % Manually compute score
                %
                in.newconname='vBestUnchosen_pos'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(    mean( [in.b{2} in.b{5}],  2)     );  % Manually compute score
                %
                in.newconname='vBestUnchosen_posMneg'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   mean([in.b{2}  in.b{5}], 2) -  in.b{3}    );  % Manually compute score
                %
                in.newconname='vBestUnchosen_negMpos'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{3}  - mean([in.b{2}  in.b{5}], 2)     );  % Manually compute score
                %
                in.newconname='cF_vBestUnchosen_negMpos';  in.cons={'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{2}  - in.b{1}   );  % Manually compute score
                %
                in.newconname='cF_vBestUnchosen_posMneg';  in.cons={'cF_vBestUnchosen_neg'; 'cF_vBestUnchosen_pos';};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{2}  - in.b{1}   );  % Manually compute score
                %
                in.newconname='vBestUnchosenPos_cFmct'; in.cons={'cF_vBestUnchosen_pos'; 'ct_vBestUnchosen_pos'};
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{1} - in.b{2}      );  % Manually compute score
            case 'vChovBU_RejOr';
                disp('Calculating new betas for vCho + vBU (RejOr) models')
                in.newconname='cF_vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{1}  -   in.b{2}  );  % Manually compute score
                %
                in.newconname='ct_vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(     in.b{3}  -   in.b{4}  );  % Manually compute score
                %
                in.newconname='vBestUnchosen_Rej';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(  mean([in.b{1}  in.b{3}],2)   );  % Manually compute score
                %
                in.newconname='vBestUnchosen_NonRej';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(  mean([in.b{2}  in.b{4}],2) );  % Manually compute score
                %
                in.newconname='vBestUnchosen_RejMOr';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(  mean([in.b{1}  in.b{3}],2)  -  mean([in.b{2}  in.b{4}],2) );  % Manually compute score
                %
                in.newconname='RejvBestUnchosen_cFMct';  in.cons= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
                for c=1:length(in.cons)
                    in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                end
                d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                d_newconbetas(2:end, end)=num2cell(   in.b{1}   -  in.b{3}  );  % Manually compute score
                
                
            otherwise,  disp('Simple effects may not be set up yet!')
                

                for o2=1:1 % [v8c]
                    %             if isempty(strfind(request.where_rois, 'm_v8c_vChosenAndposneg2')) ==0
                    
                    if strcmp(request.filename(1:3), 'v8c')==1
                        disp('Calculating new betas for v8')
                        %
                        in.newconname='subEVpos'; in.cons={'cF_subEVpos';  'cF_subEVneg'; 'ct_subEVpos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(    mean([in.b{1}  in.b{3}], 2)     );  % Manually compute score
                        %
                        in.newconname='subEVposMneg'; in.cons={'cF_subEVpos';  'cF_subEVneg'; 'ct_subEVpos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(    mean([in.b{1}  in.b{3}], 2)  - in.b{2}    );  % Manually compute score
                        %
                        in.newconname='cF_subEVposMneg'; in.cons={'cF_subEVpos';  'cF_subEVneg'; 'ct_subEVpos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(    in.b{1}  -  in.b{2}    );  % Manually compute score
                        
                        
                    end
                end
                for o2=1:1 % [v9c] - I think this is a duplicate of vChovBUposneg
                    if isempty(strfind(request.where_rois, 'm_v9c_vChosenAndposneg2')) ==0
                        in.newconname='vBestUnchosen_posMneg'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(   d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(   mean([in.b{2}  in.b{5}], 2) -  in.b{3}    );  % Manually compute score
                        %
                        in.newconname='vBestUnchosen_negMpos'; in.cons={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(     in.b{3}  - mean([in.b{2}  in.b{5}], 2)     );  % Manually compute score
                        %
                        in.newconname='cF_vBestUnchosen_negMpos';  in.cons={'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg';};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(     in.b{2}  - in.b{1}   );  % Manually compute score
                        %
                        in.newconname='cF_vBestUnchosen_posMneg';  in.cons={'cF_vBestUnchosen_neg'; 'cF_vBestUnchosen_pos';};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(     in.b{2}  - in.b{1}   );  % Manually compute score
                        %
                        in.newconname='vBestUnchosenPos_cFmct'; in.cons={'cF_vBestUnchosen_pos'; 'ct_vBestUnchosen_pos'};
                        for c=1:length(in.cons)
                            in.b{c}=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), [log.rois{r} log.roi_con_connector in.cons{c}]))));
                        end
                        d_newconbetas{1, size(d_newconbetas,2)+1}= [log.rois{r} log.roi_con_connector  in.newconname];if r==1; log.newcons=[log.newcons; in.newconname]; end
                        d_newconbetas(2:end, end)=num2cell(     in.b{1} - in.b{2}      );  % Manually compute score
                    end
                end
                
                
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


%% Get all stats (r_stats)
%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dostats=1;

% Get statistics (switch things on and off here!!!)
if dostats
    r_stats=[log.rois cell(log.n_rois, 5)]; k=1;
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
                [wr.anova]=teg_repeated_measures_ANOVA(wr.d(:,1:6) ,[2 3], {'Task' 'Choice'});
                r_stats{r, 3}=wr.anova.R;  % Row=Task, Choice, TxC; Col=F, df1, df2, p
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
                [wr.h wr.p wr.ci wr.stats]=ttest(abs(wr.TaskValtype(:, 2))   - (wr.TaskValtype(:, 1)) );
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
                [wr.h wr.p wr.ci wr.stats]=ttest(    mean([wr.ctval(:, 2)   wr.cFval(:, 2)],2)  -  wr.cFval(:, 3)   );
                r_stats{r, 4}{k,2}=[wr.stats.tstat wr.stats.df wr.p]; k=k+1;
                if wr.p<0.051;  disp([log.rois{r} '  cF_vBU_negmpos:   '  num2str(wr.stats.tstat,2) '  ,  '  num2str(wr.p,3) ' *'])
                else disp([log.rois{r} '  vBU_posMneg:'])
                end
            case 'vChovBU_RejOr';
                disp('Simple effects not set up yet!');
            otherwise disp('Simple effects stats not set up yet!');
        end
        
        %
        wr=[];
    end
end
% 
% for r=1:log.n_rois
% disp([r_stats{r,1} ': interaction p='  num2str(r_stats{r,3}(3, 4))])
% end

%% Figures

%   Col 2: Mean, Std Error, One-sample tstat, df, One-sample p  (row=contrast)
%   Col 3: 2x2 factorial
%   Col 4: Simple-effects ttests

dofig=1;
for o1=1:1  % Figure settings
        fontsize=15;
        fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms
        % fontname='Cambria';
        % fontname='Arial';
end
if dofig
    
    
    close all hidden
    f.plotcols=4;  f.figwidth= 1800; f.figheight=400; f.fontsize=15; f.fontsize_title=15;
    f.subplot_VerHorz=[0.09 0.075]; f.fig_BotTop=[0.06 0.05]; f.fig_LeftRight=[0.1 0.1];
%     f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05]; % Loads of rois
    figure('Name', 'ROI betas', 'NumberTitle', 'off', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for r=1: log.n_rois
        set(gca, 'FontSize',f.fontsize, 'FontName', fontname)
        
        switch log.FLtype
            case 'TaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_Accept'; 'cF_Reject';  'cF_Explore'; 'ct_NoBomb'; 'ct_Bomb'; 'ct_Explore' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 3,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 3,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
                ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
%                 ylabel('Parameter estimates (mean)', 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Accept', 'Reject', 'Explore'), end
                % %                     xticklabel_rotate
            case 'inTaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'in_cF_Reject';  'in_cF_Explore'; 'in_ct_Bomb'; 'in_ct_Explore' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject', 'Explore'), end
                % %                     xticklabel_rotate
                
            case 'outTaskxChoice'; % [Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'out_cF_Accept';  'out_cF_Reject'; 'out_ct_Accept'; 'out_ct_Reject' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Accept', 'Reject'), end
                % %                     xticklabel_rotate
            case 'predChoice'; % [Predicted Task x Choice]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={'cF_predReject';  'cF_predExplore'; 'ct_predReject';'ct_predExplore'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('(Predicted) Reject', '(Predicted) Explore'), end
                % %                     xticklabel_rotate
            case 'vChovBU';  % [vChosen + vBestUnchosen]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_vChosen'; 'cF_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('V(Chosen)', 'V(Best Unchosen)'), end
                % %                     xticklabel_rotate
                
            case 'vChovBUposneg'; % [vChosen + vBestUnchosen_Pos/Neg]  ###################
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
%                 wr.whichcons_name={ 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vBestUnchosen_pos';};
%                 wr.whichcons_name_label={'cF_vbu_pos'; 'cF_vbu_neg'; 'ct_vbu'};
                wr.whichcons_name={'cF_vChosen'; 'cF_vBestUnchosen_pos'; 'cF_vBestUnchosen_neg'; 'ct_vChosen'; 'ct_vBestUnchosen_pos'};
                wr.whichcons_name_label={'cF_vcho'; 'cF_vbu_pos'; 'cF_vbu_neg'; 'ct_vcho'; 'ct_vbu'};
                % 
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0)); 
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr(r_stats{r,2}(wr.whichcons,2)', r_stats{r,2}(wr.whichcons,1)', 'y'); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:length(wr.whichcons_name), 'XTickLabel', wr.whichcons_name_label,  'FontSize',f.fontsize, 'FontName', fontname);  
                % set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Conflict'; 'Control'},  'FontSize',f.fontsize, 'FontName', fontname);  
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
%                 if r==log.n_rois,  legend(wr.whichcons_name), end
            case 'vChovBU_RejOr';   % [RejOr vBU]  ###################
                
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name={ 'cF_vChosen'; 'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen' ;  'ct_vChosen'; 'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen' };
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 3,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 3,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois, legend('V(Chosen)', 'Reject V(Best Unchosen)', 'Non-Reject V(Best Unchosen)'), end
                % %                     xticklabel_rotate
            case 'vGamble_RejOr'
                
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'ct_Rej_vGamble'; 'ct_NonRej_vGamble'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( reshape( r_stats{r,2}(wr.whichcons,2), 2,2)' ,   reshape( r_stats{r,2}(wr.whichcons,1), 2,2)'     ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:2, 'XTickLabel', {'Exp' 'Ctrl' },  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
%                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
%                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject V(Gamble)', 'Non-Reject V(Gamble)'), end
                % %                     xticklabel_rotate
                
                
                
                
            case 'vChovBU_cFRejOr'
                
                subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                wr.whichcons_name= {'cF_vChosen'; 'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_vChosen'; 'ct_vBestUnchosen'};
                wr.whichcons_name_label= {'cF_vCho'; 'cF_Rej_vbu'; 'cF_NonRej_vbu'; 'ct_vCho'; 'ct_vbu'};
                wr.whichcons=cell2mat(cellfun(@(x)find(strcmp(log.cons, x)), wr.whichcons_name, 'UniformOutput',0));
                if length(wr.whichcons)~=length(wr.whichcons_name); error('Could not find requested contrast for plotting!'); end
                barwitherr( r_stats{r,2}(wr.whichcons,2) , r_stats{r,2}(wr.whichcons,1) ); % Mean +/- 1 Std Err
                set(gca, 'FontSize',f.fontsize, 'FontName', fontname, 'TickDir','out',  'XTick', 1:length(wr.whichcons_name_label), 'XTickLabel', wr.whichcons_name_label,  'FontSize',f.fontsize, 'FontName', fontname);  % xlim([0.3 2.7])
                
                
                %                 ylim([min(r_stats{r,2}(wr.whichcons,1))  - max(r_stats{r,2}(wr.whichcons,2)) - 0.5 max(r_stats{r,2}(wr.whichcons,1))  + max(r_stats{r,2}(wr.whichcons,2)) + 0.5])
                title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
                %                 ylabel(sprintf('Parameter estimates\n(mean)'), 'FontSize',f.fontsize, 'FontName', fontname); set(gca,'FontSize',f.fontsize, 'FontName', fontname, 'LineWidth', 0.8);
                if r==log.n_rois,  legend('Reject V(Gamble)', 'Non-Reject V(Gamble)'), end
                % %                     xticklabel_rotate
                
              
                
                
                
                
                
            otherwise error('Which plots? '); 
                
        end
        
        
        % Signicant?
        %         [h p ci stats]=ttest(cell2mat(d_betas(2:end, 1+[129 136])));
        %         stats
        %         p
        
        
        
        
%         ylim([-5 2])
        xlim([0.5 2.5])
        
        
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
            %         title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
            %         ylabel('Parameter estimates', 'FontSize',f.fontsize, 'FontName', fontname)
            %         set(gca,'xticklabel', wr.whichcons_namelabel, 'FontSize',f.fontsize, 'FontName', fontname)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
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
            %         title(log.rois{r},  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', fontname,'Interpreter', 'none')
            %         ylabel('Parameter estimates', 'FontSize',f.fontsize, 'FontName', fontname)
            %         set(gca,'xticklabel', wr.whichcons_namelabel)   % 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
            % %         xticklabel_rotate
            %
        end
        
        
        %
        %         wr.whichcons_name={'vBUpos';'vBUneg'; 'neg-pos'};
        % %         wr.whichcons_name={'+ Counterfactual value';'- Counterfactual value'};
        %         wr.whichcons_name={'Positive values';'Negative values'};
        %         set(gca,'TickDir','out',  'XTick', 1:length(wr.whichcons_name), 'XTickLabel', wr.whichcons_name,  'FontSize',10, 'FontName', fontname);
        %         ylim([0 1])
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
% 	wr.whichcons_namelabel={'Gamble to Info';'Info to Outcome'; 'Gamble to Outcome'; 'Left ventral striatum';'Right ventral striatum';'vmPFC';};
%     
%     %
%     subtightplot(ceil((log.n_rois+1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); axis off;
%     offset=0.15; fontsize=20;  % Legend
%     set(gca,'FontSize',35, 'FontName', fontname, 'LineWidth', 0.8);
%     for l=1:length(wr.whichcons_namelabel)
%         text(0,1-offset*l, wr.whichcons_namelabel{l},'FontSize', fontsize, 'FontName', fontname, 'FontWeight','bold')
%     end
    end
end


% error('Done with stats and fig')

%% Battery correlations (specific)
%       r_batcorr{k}: row= beh, col= roi/con




% Certain ROIs?
request.rois=log.rois;
d_meanbetas=[d_betas(1,2:end)'  num2cell(mean(cell2mat(d_betas(2:end, 2:end)))') num2cell(std(cell2mat(d_betas(2:end, 2:end)))'./sqrt(log.n_subjs))  num2cell( ttest(cell2mat(d_betas(2:end, 2:end)))' )  ]; 
[h p]=ttest(cell2mat(d_betas(2:end, 2:end))); p(p>0.1)=nan; d_meanbetas(:,4)=num2cell(p'); 
openvar d_meanbetas  % mean beta
% error
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

end

    

% request.cor_inhib.rois={'Hippocampus Left (Task x Choice)'; 'Hippocampus Right (Task x Choice)';};
request.cor_inhib.rois=request.rois(1:2);
request.cor_inhib.rois=request.rois;
switch log.FLtype
    case 'TaskxChoice';   request.cor_inhib.con={'cF_Reject'; 'cF_Rej-Exp';'Rej_cF-ct'};
    case 'inTaskxChoice';   request.cor_inhib.con={'in_cF_Reject'; 'cF_Rej-Exp';'Rej_cF-ct'};
    case 'outTaskxChoice';   request.cor_inhib.con={'in_cF_Reject'; 'Rej_cF-ct'};
    case 'predChoice';   
        
        request.cor_inhib.con={'cF_predReject'; 'cF_predExplore';'predReject_cF-ct'; 'cF_predRejMExp'};
    case 'vChovBUposneg';  request.cor_inhib.con={'cF_vBestUnchosen_posMneg'; 'cF_vBestUnchosen_pos';  'cF_vBestUnchosen_neg'; 'vBestUnchosen_posMneg'; 'vBestUnchosenPos_cFmct'; };
    case 'vChovBU'; request.cor_inhib.con={'cF_vBestUnchosen';'cF_vCho-vBU'; 'vChosen_cF-ct'; 'vBestUnchosen_cF-ct';'vCho-vBU_cFmct'; 'cF_vCho-vBUn'; 'cF_vCho-vBUp'; 'ct_vCho-vBU'; 'vCho_cF-ct'};
%         request.cor_inhib.con={'cF_vChosen';'cF_vBestUnchosen'; 'ct_vChosen';'ct_vBestUnchosen'; 'vBestUnchosen_cFmct'; 'vChosen_cFmct'; };            
    case 'vChovBU_RejOr'; request.cor_inhib.con= {'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen';'cF_vBestUnchosen_RejMOr';  'vBestUnchosen_RejMOr'; 'RejvBestUnchosen_cFMct'};
%         request.cor_inhib.con={'cF_vChosen';'cF_Rej_vBestUnchosen'; 'cF_NonRej_vBestUnchosen'; 'ct_vChosen';   'ct_Rej_vBestUnchosen'; 'ct_NonRej_vBestUnchosen'; };
%             request.cor_inhib.con={'cF_vCho-RejvBU'; 'cF_vCho-NonRejvBU' ; 'vChosen_cF-ct'};
    case 'vGamble_RejOr'; request.cor_inhib.con= {'cF_Rej_vGamble'; 'cF_NonRej_vGamble'; 'cF_vGam_RejMNonRej'; 'Rej_vGam_cFMct'; 'cF_Rej_vGam_mOthers'};
    case 'vChovBU_cFRejOr'; request.cor_inhib.con =  {'cF_vCho-RejvBU'; 'cF_vCho-NonRejvBU'; 'cF_vBU-Rej-NonRej'; 'vChosen_cF-ct'};
end
% request.cor_inhib.con= {'cF_Rej_vGamble' 'cF_vGam_RejMNonRej';  'Rej_vGam_cFMct'; 'cF_Rej_vGam_mOthers'};  % v25
% request.cor_inhib.con={    % see log.cons for list  ############
% % 'cF_vBestUnchosen_negMpos'
% 'cF_vBestUnchosen_pos';  'cF_vBestUnchosen_neg';
% 'vBestUnchosen_posMneg'; 'vBestUnchosenPos_cFmct';
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
        
        if request.kendalls_correlation
            [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  disp('Kendalls correlation'); %
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
    end
end
openvar r_batcorr{1,1}  % Put probabilities in 2nd col if necessary

% error('Done with BehInhib correlation battery')


% 
% % AD HOC PLOTS 
% roicon= 'HPC_aL_tc-cF_Rej-Exp';
% beh='Trait anxiety scores';
% % beh='State anxiety scores';
% % % beh='Probability of following null information\n from Exploring in experimental task';
% beh='B.Bis';
% 
% % --- 
% 
% wr.rc= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),roicon)) ));
% wr.b= cell2mat( d_betas(2:end,  find(strcmp(d_betas(1,:),beh)) ));
% 
% 
% figure('color', 'w')
% scatter(wr.rc, wr.b), lsline
% xlabel(roicon,'FontSize',20), ylabel(beh,'FontSize',20)
% if      request.kendalls_correlation==0, [r p]= corr(wr.rc, wr.b); title(['r='  num2str(r) ', p=' num2str(p)],'FontSize',20)
% else      [r p]= corr(wr.rc, wr.b,'type', 'Kendall'); title(['tau='  num2str(r) ', p=' num2str(p)],'FontSize',20)
% end
% set(gca,'FontSize',20)

% (2) Exploration battery  #########################################
request.cor_explore.beh=d_betas(1,  find(cellfun(@(x)~isempty(x), strfind(d_betas(1, :), 'BEH_Explore')))+1:end)';
request.cor_explore.con=log.cons;
request.cor_explore.rois=request.rois;
switch log.FLtype
    case 'TaskxChoice';   request.cor_explore.con={'cF_Explore';'ct_Explore'; 'Explore'; 'RejmExp'};
    case 'vChovBU';  request.cor_explore.con={'cF_vChosen';'cF_vBestUnchosen'; 'ct_vChosen';'ct_vBestUnchosen'; 'vBestUnchosen'; 'vChosen';};
    case 'vChovBUposneg';
    case 'vChovBU_RejOr';
end
request.cor_explore.con={
    %     'vChosen';  'vBestUnchosen_pos';   'cF_vBestUnchosen_neg';
    % 'vBestUnchosen_posMneg';   'vBestUnchosen_negMpos';
    % 'cF_vExplore';'ct_vExplore'
    % 'cF_CueValue';'cF_OutcomeMagnitude';'ct_CueValue';'ct_OutcomeMagnitude'
    };

rc=1; request.cor_explore.roicons={}; % Compile roicons + execute battery
for r=1: length(request.cor_explore.rois)    
    for c=1:length(request.cor_explore.con)
        request.cor_explore.roicons{rc,1}=[request.cor_explore.rois{r}  log.roi_con_connector    request.cor_explore.con{c}];  rc=rc+1;
    end
end
k=2;  r_batcorr{k,1}=[ [{' ' } request.cor_explore.roicons'];   [request.cor_explore.beh  cell(length(request.cor_explore.beh),  length(request.cor_explore.roicons))] ];
for b=2:size(r_batcorr{k,1},1)
    for rc=2:size(r_batcorr{k,1},2)
        wbc.beh=cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{b,1}), 1, 'first')));
        wbc.roiconbeta= cell2mat(d_betas(2:end, find(strcmp(d_betas(1,:), r_batcorr{k,1}{1,rc}))));
        
        
        [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta);  % 
        [wbc.r wbc.p]=corr(wbc.beh, wbc.roiconbeta,'type', 'Kendall');  % 
        wbc.stat=wbc.r;   % Print r or p?5
%         wbc.stat=wbc.p;
        wbc.stat=['r=' num2str(wbc.r,2) ', p=' num2str(wbc.p,2)];

        if wbc.p<0.01;  r_batcorr{k,1}{b,rc}= [wbc.stat  ' **'];;  % [num2str(wbc.stat,3) ' **'];
        elseif wbc.p<=0.05;  r_batcorr{k,1}{b,rc}= [wbc.stat ' *'];  %[num2str(wbc.stat,3) ' *'];
%         elseif wbc.p<0.1; r_batcorr{k,1}{b,rc}=wbc.stat ;
%         else r_batcorr{k,1}{b,rc}=wbc.stat ;
        end
    end
end
openvar r_batcorr{2,1}  % Put probabilities in 2nd col if necessary
% error('Done with Explore correlation battery')

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

request.plotcorr={'o3'}; 
for o1=1:1 % [c13 FL & c13 rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o1', request.plotcorr))==1
        
        % LEFT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_stc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left Inferior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
%         r_corr{k,1}='Parameter estimates from Left Superior \n  Hippocampus (Exp, Reject > Explore)'; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left Inferior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
%         r_corr{k,1}='Parameter estimates from Left Superior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
        %
        % RIGHT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aR_stc'; % wc.roiname='Hippocampus Left (Task x Choice)';
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];
        wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Right Inferior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
%         r_corr{k,1}='Parameter estimates from Right Superior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Right Inferior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;
%         r_corr{k,1}='Parameter estimates from Right Superior \n Hippocampus (Exp, Reject > Explore)'; k=k+1;

    
    
    end
end
for o2=1:1 % [c13 FL & Subfield ROIs: HPC TxC Beh inhibition]
    
    if sum(strcmp('o2', request.plotcorr))==1
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.Trait';   wc.behname='Trait anxiety scores';
        wc.roiname='rCA1_aL';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)'; k=k+1;
        %
        wc.roiname='rCA1_aR'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        %
        wc.roiname='rCA3_aL'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)'; k=k+1;
        %
        wc.roiname='rCA3_aR'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Right \n CA3 (Exp, Reject > Explore)'; k=k+1;
        
    
        
        % Comparing all ROIs, cF Rej-Exp ##############################
        wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.roiname='rCA1_aL';  wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left \n CA1 (Exp, Reject > Explore)'; k=k+1;
        %
        wc.roiname='rCA1_aR'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Right \n CA1 (Exp, Reject > Explore)'; k=k+1;
        %
        wc.roiname='rCA3_aL'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left \n CA3 (Exp, Reject > Explore)'; k=k+1;
        %
        wc.roiname='rCA3_aR'; wc.betaname=[wc.roiname '-cF_Rej-Exp'];  
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Right \n CA3 (Exp, Reject > Explore)'; k=k+1;
      
        
    end
end
for o3=1:1 % [v28g FL & v28g rois: HPC TxC Beh inhibition]
    
    if sum(strcmp('o3', request.plotcorr))==1
        
        % Inferior HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_stc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Left Inferior \n Hippocampus (Exp, Reject > Explore)'; 
        wc.roiconname='Left Inferior Hippocampus betas \n  (Exp, Reject > Explore)';
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        % % % ----------------
        wc.roiname='HPC_aR_stc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Right Inferior \n Hippocampus (Exp, Reject > Explore)'; 
wc.roiconname='Left Inferior Hippocampus betas \n  (Exp, Reject > Explore)';
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
        
        % Superior HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_aL_sc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Left Superior \n Hippocampus (Exp, Reject > Explore)'; 
wc.roiconname='Left Superior Hippocampus betas \n  (Exp, Reject > Explore)';
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        % % % --------------------------------------
        wc.roiname='HPC_aR_sc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_Rej-Exp'];   % observed choice 
%         wc.betaname=[wc.roiname '-cF_predRejMExp'];   % predicted choice 
%         wc.roiconname='Parameter estimates from Right Superior \n Hippocampus (Exp, Reject > Explore)'; 
wc.roiconname='Right Superior Hippocampus betas \n  (Exp, Reject > Explore)';
        %
        wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        %
        wc.behname= 'B.Bis'; wc.behname='BIS scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        if request.kendalls_correlation,[wc.r wc.p]=corr(wc.beta, wc.beh, 'type', 'Kendall');  else  [wc.r wc.p]=corr(wc.beta, wc.beh);  end; r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}=wc.roiconname; k=k+1;
        
        
    
    
    end
end

% Automated plot all requested correlations
f.scatter_dotsize=100; 
f.scatter_dotsize=50; 
f.scatter_linewidth=4;   f.FontSize=10; f.FontName='PT Sans Caption';   r_corr{:,1}
f.plotcols=6;  f.figwidth= f.plotcols*400; f.figheight=ceil(size(r_corr,1)/ (f.plotcols))*400;
f.subplot_VerHorz=[0.2 0.07]; f.fig_BotTop=[0.2 0.1]; f.fig_LeftRight=[0.07 0.05];
figure('Name', 'Beh - beta ', 'NumberTitle', 'off', 'Position', [100 250 f.figwidth f.figheight], 'Color', 'w');  k=1;
for c=1:size(r_corr,1)
    disp(r_corr{k,1})
    subtightplot(ceil(size(r_corr,1)/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
%     scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2)); lsline
    scatter(r_corr{c,5}(:,1),  r_corr{c,5}(:,2), f.scatter_dotsize,'LineWidth', 3); h=lsline; set(h,'LineWidth', f.scatter_linewidth)
    
    
    
%     ylabel(r_corr{c,2}, 'FontSize',15, 'FontName', fontname);  
    ylabel(sprintf(r_corr{c,2}), 'FontSize',f.FontSize, 'FontName', f.FontName)
%     xlabel(sprintf(['Parameter estimates from\n' r_corr{c,1}]), 'FontSize',15, 'FontName', fontname)
    xlabel(sprintf(r_corr{c,1}), 'FontSize',f.FontSize, 'FontName', f.FontName)
    set(gca,'FontSize',f.FontSize, 'FontName', f.FontName, 'LineWidth', 0.8);
    if r_corr{c,6}==1; 
        if request.kendalls_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', fontname)    
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2) ' *'], 'FontSize',15, 'FontName', fontname)
        end
    else
        if request.kendalls_correlation, title(['tau='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', fontname) 
        else title(['r='  num2str(r_corr{c,3},2) ', p=' num2str(r_corr{c,4},2)], 'FontSize',15, 'FontName', fontname)
        end
    end
    
    
% %     Extra adhoc stuff
%     xlim([-1 2])
xlim('auto')
    
    % 
    k=k+1;
end

%%


for oz=1:1  % OLD correlation plots (Pre Manuscript v2-6)
for o1=1:1 % [c3 FL & c3 rois: HPC TxC Beh inhibition]
    this=0;
    if this
        % LEFT HPC, cF Rej-Exp ##############################
        wc.roiname='HPC_L_stc'; %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State';   wc.behname='State anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}='Left Hippocampus (Conflict Reject > Explore)'; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Left Hippocampus (Conflict Reject > Explore)'; k=k+1;
        %
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'B.Drive'; wc.behname='BIS/BAS Drive scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
%         [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}='Left Hippocampus (Conflict Reject > Explore)'; k=k+1;
        
        % HPC TxC  && Modelling loss aversion ##################
        % LEFT HPC, cF Rej-Exp
        wc.roiname='HPC_L_stc'; % wc.roiname='Hippocampus Left (Task x Choice)';
%         wc.betaname=[wc.roiname '-cF_Rej-Exp'];
%         wc.behname= 'st.State';   wc.behname='State anxiety scores'; %
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}='Parameter estimates from Left Hippocampus, Conflict Reject > Explore'; k=k+1;
        %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left Hippocampus, Conflict Reject > Explore'; k=k+1;
        
%         % RIGHT HPC, cF Rej-Exp
        wc.roiname='HPC_R_stc'; % wc.roiname='Hippocampus Right (Task x Choice)';
%         wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.State'; wc.behname='State anxiety scores';
%         wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
%         [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
%         r_corr{k,1}='Parameter estimates from Left Hippocampus, Conflict Reject > Explore'; k=k+1;
%         %
        wc.betaname=[wc.roiname '-cF_Rej-Exp']; wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';%
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.behname)));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Parameter estimates from Left Hippocampus, Conflict Reject > Explore'; k=k+1;
    end
end
for o1=1:1 % [c3 FL & c3 rois: Frontal-striatal Explore]
    this=0;
    if this
        % Winnings on task ##############################
        wc.roiname='Striatum_R';  wc.roiname='Right striatum';  %%%
        wc.betaname=[wc.roiname '-Explore']; wc.behname= 'o.allfMRI';   wc.behname='Overall money won\n on task (points)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Right striatum on Explore trials'; k=k+1;
        %
        wc.roiname='SupMFG';  wc.roiname='Superior/Mid frontal gyrus';  %%%
        wc.betaname=[wc.roiname '-Explore']; wc.behname= 'o.allfMRI';   wc.behname='Overall money won\n on task (points)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Superior/Mid frontal gyrus on Explore trials'; k=k+1;
        %
        wc.roiname='Precuneus';  wc.roiname='Precuneus';  %%%
        wc.betaname=[wc.roiname '-Explore']; wc.behname= 'o.allfMRI';   wc.behname='Overall money won\n on task (points)';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Precuneus on Explore trials'; k=k+1;
        
        
        % Following information ##############################
        wc.roiname='BA10';  wc.roiname='BA10';  %%%
        wc.betaname=[wc.roiname '-Explore']; wc.behname= 'info.allfMRI';   wc.behname='Probability of following information\n from Exploring';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='BA10 on Explore trials'; k=k+1;
        %
        wc.roiname='Precuneus';  wc.roiname='Precuneus';  %%%
        wc.betaname=[wc.roiname '-Explore']; wc.behname= 'info.allfMRI';   wc.behname='Probability of following information\n from Exploring';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Precuneus on Explore trials'; k=k+1;
    end
end
for o1=1:1 % [v6c FL & c3 rois: Value tracking Beh inhibition]
    this=0;
    if this
        % Value tracking ##############################
        wc.roiname='HPC_R_satc';  %  wc.roiname='Hippocampus Right (Task x Choice)';  %%%
        %     wc.betaname=[wc.roiname '-vBestUnchosen_posMneg'];
        wc.betaname=[wc.roiname '-cF_vBestUnchosen_posMneg'];
        wc.behname= 'st.State'; wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        %     r_corr{k,1}='Preferential tracking of + > - counterfactual value \n  in the right inferior hippocampus'; k=k+1;
        r_corr{k,1}='Difference in tracking of positive > negative V(Best Unchosen) \n  in the right inferior hippocampus (experimental task)'; k=k+1;
        % --------------------------------
        wc.roiname='HPC_R_satc';  % wc.roiname='Hippocampus Right (Task x Choice)';  %%%
        %     wc.betaname=[wc.roiname '-vBestUnchosen_posMneg'];
        wc.betaname=[wc.roiname '-cF_vBestUnchosen_posMneg'];
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        %     r_corr{k,1}='Preferential tracking of + > - counterfactual value \n  in the right inferior hippocampus'; k=k+1;
        r_corr{k,1}='Difference in tracking of positive > negative V(Best Unchosen) \n  in the right inferior hippocampus (experimental task)'; k=k+1;
        % --------------------------------
        wc.roiname='HPC_L_satc';  %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        %     wc.betaname=[wc.roiname '-vBestUnchosen_posMneg'];
        wc.betaname=[wc.roiname '-cF_vBestUnchosen_posMneg'];
        wc.behname= 'infoexnull.cF_all'; wc.behname='Probability of following null information\n from Exploring in experimental task';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        %     r_corr{k,1}='Preferential tracking of + > - counterfactual value \n  in the left inferior hippocampus'; k=k+1;
        r_corr{k,1}='Difference in tracking of positive > negative V(Best Unchosen) \n  in the left inferior hippocampus (experimental task)'; k=k+1;
        %         --------------------------------
    end
end
for o1=1:1 % [v6 FL & v4 rois: beh inhibition]
    this=0;
    if this
        wc.roiname='Caudate';  wc.roiname='Caudate';  %%%
        wc.betaname=[wc.roiname '-vBestUnchosen_posMneg'];  wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Preferential tracking of positive over negative \n counterfactual values in the caudate'; k=k+1;
        
        %
        wc.roiname='Caudate';  wc.roiname='Caudate';  %%%
        wc.betaname=[wc.roiname '-vBestUnchosen_posMneg'];  wc.behname= 'st.State'; wc.behname='State anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Preferential tracking of positive over negative \n counterfactual values in the caudate'; k=k+1;
    end
end
for o1=1:1 % [v10c FL & c3 rois: Value tracking Beh inhibition]
    this=0;
    if this
        % Value tracking ##############################
        wc.roiname='HPC_aR_satc';  %  wc.roiname='Hippocampus Right (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_vBestUnchosen_RejMOr'];
        wc.behname= 'st.Trait'; wc.behname='Trait anxiety scores';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Exp - Difference in parameter estimates for V(Best Unchosen) \n tracking (Reject > Non-Reject) '; k=k+1;
          % --------------------------------
        wc.roiname='HPC_aL_satc';  %  wc.roiname='Hippocampus Left (Task x Choice)';  %%%
        wc.betaname=[wc.roiname '-cF_vBestUnchosen_RejMOr'];
        wc.behname= 'infoexnull.cF_all'; wc.behname='Probability of following null information\n from Exploring in experimental task';
        wc.beta=cell2mat(d_betas(2:end,  strcmp(d_betas(1,:), wc.betaname))); wc.beh=cell2mat(d_betas(2:end,  find(strcmp(d_betas(1,:), wc.behname),1,'first')));
        [wc.r wc.p]=corr(wc.beta, wc.beh); r_corr(k,1:2)=[{wc.betaname} {wc.behname}]; r_corr(k,3:4)=num2cell([wc.r wc.p]); r_corr{k,5}=[wc.beta wc.beh]; r_corr{k, 6}=wc.p<=0.05; if r_corr{k, 6}==0;  r_corr{k, 6}=(wc.p<=0.1)*0.5; end;
        r_corr{k,1}='Exp - Difference in parameter estimates for V(Best Unchosen) \n tracking (Reject > Non-Reject) '; k=k+1;
        % --------------------------------
    end
end
end 


% save([request.where_rois filesep 'Stats results(' date ')'],'r_corr', 'r_batcorr', 'log', 'request')
% save([request.where_rois filesep 'Stats results(' date ')'], 'r_stats','r_corr', 'r_batcorr', 'log', 'request')
