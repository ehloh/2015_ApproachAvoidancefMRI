% Extract ROI betas using SPM functions (subject-level) for (flexible) Events model
%       Beta within each ROI is then modelled with GLMs describing each trial
clear all;close all hidden; clc

% Request
log.specificsubjects={};
request.loadextractedbetas= '(10-May-2014)'; % empty to extract
request.betatype=1;   % 1= cons, 2=betas 
%
request.beh_filename='Behavioural event data (09-May-2014)';
request.glmtype='linear';  % linear/interactions

for o1=1:1 % General settings and specifications 
    
    log.FLmodel='f3_1_Event';         % ##### WHICH FL MODEL #########
    switch request.betatype
        case 1; log.betatype='con';
        case 2; log.betatype='beta';
    end
    % Folders
    w.w=pwd; if strcmp(w.w(1), '/')==1; 
        where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI'; where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data';  where.data_brain=[where.experiment_folder filesep '1 Brain data'];  
        where.behmodel='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';    
    else where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]'; 
        where.behmodel='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
        where.data_brain='C:\Users\eloh\Desktop\2 [Explore]\1 Brain data';
    end;   addpath(where.where); addpath(where.behmodel);
    where.data_beh=[where.where filesep '1 Behavioural data'];
    if isdir([where.experiment_folder filesep '2 Second level results' filesep log.FLmodel])==0; mkdir([where.experiment_folder filesep '2 Second level results' filesep log.FLmodel]); end
    log.SLfol=[where.experiment_folder filesep '2 Second level results' filesep log.FLmodel filesep]; if isdir(log.SLfol)==0; mkdir(log.SLfol); end
    log.FLfol=['m_' log.FLmodel ' Contrasted'];
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    request.specificsubjects=log.specificsubjects;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end; disp(' ');
    disp(['             First level model:              ' log.FLmodel])
    disp(['             FL folder:      ' log.FLfol])
    disp(['             SL folder:      ' log.SLfol]); disp(' ');
%     input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% ROI settings

% Where are the ROI images?
% log.roi_fol=[log.SLfol 'c3 ROI battery'];
log.roi_fol=[log.SLfol 'Choice Anat battery'];
% log.roi_fol=[log.SLfol 'Test'];

% ROI details (List of ROIs + requested order)
for o1=1:1
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
    % Battery in order
        'BA10_L';  'BA10_R'; 'BA46_L';'BA46_R'; 'Cingulum_Ant_L';  'Cingulum_Ant_R';               % Cortex
        'Striatum_L'; 'Striatum_R'; 'Caudate_L';  'Caudate_R'; 'Putamen_L'; 'Putamen_R'; % Striatum
        'NAccCore_L';  'NAccCore_R'; 'NAccShell_L'; 'NAccShell_R'; 'NAcc_L'; 'NAcc_R';
        'HPC_aL';  'HPC_aR';   'Hippocampus_L';  'Hippocampus_R';    % HPC/Amygdala
        'Amygdala_L';  'Amygdala_R';
        'Insula_L';  'Insula_R'; 'SNVTA_L'; 'SNVTA_R';  % Exploratory/various
        };

end
% log.batteryname='c3battery';
log.batteryname='choiceanat';
eval(['log.rois=log.roilist.' log.batteryname ';'])
% log.rois={'BA10';'BA46'};
    
if isempty(request.loadextractedbetas)==1
    log.roi_files=cellstr(spm_select('List', log.roi_fol, '.nii$')); log.n_rois=length(log.roi_files); % Display ROI details
    if isfield(log, 'batteryname')==1; disp(['######    Battery name = ' log.batteryname '  ######']); end
    disp('------  ROI files       ---------------      ROI short names ------'); disp([log.roi_files   log.rois]); input('OK?    ');
    d_roimasks=cell(length(log.roi_files),2); % Load ROI images
    for r=1:log.n_rois
        d_roimasks{r,1}=log.rois{r};
        d_roimasks{r,2}=spm_read_vols(spm_vol([log.roi_fol filesep log.roi_files{r}]));
    end
end

for o2=1:1 % % [ Load behavioural and beta data ] ##########

%% (1) Load behavioural data

if isempty(request.beh_filename)==1;
    for o1=1:1 % Columns for behavioural data
        % Columns from functional raw data
        col.f.Trialtype=1;
        col.f.NTokens=2;
        col.f.EnvThreat=3;
        col.f.Task=5;
        col.f.Trialnum=7;
        col.f.Choice=8;
        col.f.Choice2=10;
        col.f.RT1=9;
        col.f.TrialValid=13;
        col.f.OutcomePres=14;
        col.f.OutcomeMagnitude=15;
        col.f.ExploredBomb=18;
        
        % Columns for concise data (Save only details requested)
        col.Choice=1;
        col.EnvThreat=2;
        col.NTokens=3;
        col.pLoss=4;
        col.Entropy=5;
        col.VExplore=6;
        col.EntropyNTok=7;
        col.EV=8;
        col.Task=9;
        %
        col.OutcomePres=10;
        col.OutcomeMagnitude=11;
        col.OutcomeMean=12;
        col.OutcomeVariance=13;
        col.ExploredBomb=14;  % Check: Is this variable actually correct?
        %
        col.TrialValid=15;
        col.Trialnum=17;    % Running trial number - must edit
        col.Trialtype=18;
        %
        col.RT1=19;
        col.Choice2=20;
    end
    subjdata=cell(log.n_subjs,2);
    for s=1:log.n_subjs
        ws=load([where.data_beh filesep log.subjects{s} filesep log.subjects{s} '_file_taskfMRI.mat']);
        ws.rd=ws.alldata;
        ws.rd(:,col.f.Trialnum)=1:size(ws.rd,1);
        
        % Mark variables (include only necessary)
        ws.d=nan(size(ws.rd, 1), 19);
        ws.d(:, col.Trialtype)=ws.rd(:, col.f.Trialtype );
        ws.d(:, col.EnvThreat)=ws.rd(:, col.f.EnvThreat);
        ws.d(:, col.NTokens )=ws.rd(:, col.f.NTokens );
        ws.d(:, col.Task )=ws.rd(:, col.f.Task);
        ws.d(:, col.Trialnum )=ws.rd(:, col.f.Trialnum );
        ws.d(:, col.Choice )=ws.rd(:, col.f.Choice );
        ws.d(:, col.Choice2 )=ws.rd(:, col.f.Choice2 );
        ws.d(:, col.RT1)=ws.rd(:, col.f.RT1);
        ws.d(:, col.TrialValid )=ws.rd(:, col.f.TrialValid  );
        ws.d(:, col.OutcomePres)=ws.rd(:, col.f.OutcomePres);
        ws.d(:, col.OutcomeMagnitude)=ws.rd(:, col.f.OutcomeMagnitude);
        ws.d(:, col.ExploredBomb)=ws.rd(:, col.f.ExploredBomb);
        ws.cf=ws.d(ws.d(:,col.Task)==1, :);
        ws.ct=ws.d(ws.d(:,col.Task)==2, :);
        [ ws.cf] = fpar_conflict( ws.cf, col);
        [ ws.ct] = fpar_control( ws.ct, col);
        ws.d=sortrows([ws.cf; ws.ct], col.Trialnum);
        
        % Check explored bomb?
        if sum(ws.rd(:,col.f.ExploredBomb)==1  & ws.rd(:, col.f.Choice)~=3)~=0; disp(['Check column for ExploredBomb for subject ' num2str(s) '(' log.subjects{s} ')']); end
        
        %
        subjdata{s,1}=log.subjects{s};
        subjdata{s,2}=ws.d;
        ws=[];
    end
    save([log.SLfol filesep 'Behavioural event data (' date ')'], 'subjdata', 'col'); % Save
else
    w=load([log.SLfol filesep  request.beh_filename    '.mat']); subjdata=w.subjdata; col=w.col;
     if isempty(request.specificsubjects)==0    % Implement subject selections requested in master script
         [log.subjects log.n_subjs w.subjdata] = f_selectsubjects([{' ' '  ' }; subjdata], request.specificsubjects, [{'Subject'    'OK'};  [subjdata(:,1) num2cell(ones(size(subjdata,1),1))]], 'OK'); subjdata=w.subjdata(2:end,:);
     end
end
    
%% (2) Extract betas: Load image data for each event and extract 
%       d_event = roi x subject x trial (event allocation within design not specified)

if isempty(request.loadextractedbetas)==1
    d_event=nan(log.n_rois, log.n_subjs, 100);
    for s=1:log.n_subjs
        ws.c=clock;  disp(['Subject ' num2str(s) '  (' log.subjects{s} ')            [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  --------------'])
        ws.FLfol=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.FLfol filesep];
        ws.s=load([ws.FLfol 'SPM.mat']); ws.cons=ws.s.SPM.xCon;
        
        % Load all contrasts + calculate mean beta per roi
        for c=1:size(ws.cons,2)
            
            % Extract from what ? 
            switch log.betatype
                case  'con'; wc=spm_read_vols(spm_vol([ws.FLfol ws.cons(c).Vcon.fname]));  if c==1; disp('Loading cons!'); end
                case  'beta'; wc=spm_read_vols(spm_vol([ws.FLfol  'beta_' char(cellfun(@(x)x(2:end),  {num2str(10000+c)},'UniformOutput',0)) '.img']));  if c==1; disp('Loading betas!'); end
            end
            
            for r=1:log.n_rois % Calculate mean beta for each roi
                d_event(r, s, c)= nanmean( wc(d_roimasks{r,2}~=0));
            end
            
            % Check contrast name
            if strcmp(ws.cons(c).name, ['tr' num2str(c)])~=1; disp('Discrepancy in contrast name!'); end
            if c/100==round(c/100); disp(c); end
            wc=[];
        end

        ws=[];
    end
    
    % Save
    save([log.roi_fol filesep 'Extracted event betas from ' log.betatype ' (' date ')'], 'd_event', 'log');  try % Notify researcher
        f_sendemail('kurzlich', strcat(['Analysis script is complete: extracting betas from all event contrasts (' mfilename ')']), ' ',1);
    end
    disp('======================================================='); w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end; disp(['             First level model:      ' log.FLmodel]); disp('=======================================================')
else disp(['Loading dataset previously extracted ' request.loadextractedbetas ]); 
    load([log.roi_fol filesep 'Extracted event betas from ' log.betatype ' ' request.loadextractedbetas '.mat'])
    if isempty(request.specificsubjects)==0 % Implement subject selections requested in script
        log.orig.subjects=log.subjects;  o.d_event=d_event;    d_event=nan(size(o.d_event,1), length(request.specificsubjects),  size(o.d_event,3));
        for s=1:length(request.specificsubjects)
            if sum(strcmp(log.subjects, request.specificsubjects{s}))~=1; error(['Invalid subject requested: ' request.specificsubjects{s}]); end
            d_event(:, s, 1:size(o.d_event(:,  find(strcmp(log.subjects, request.specificsubjects{s})), :),3))=o.d_event(:,  find(strcmp(log.subjects, request.specificsubjects{s})), :);
        end
        log.subjects=request.specificsubjects; log.n_subjs=length(log.subjects);  clear('o')
    end
end
col.Subject=size(subjdata{1,2},2)+1;
for s=1:log.n_rois
    subjdata{s,2}(:, col.Subject)=s;
end

for o2=1:1 % [DISUSED FOR NOW]  %% (2b) Pre-process betas for analysis

% % Selection settings ############
% request.flex_mintrials4inclusion=5; % empty to omit trialwise filtering. Trialwise filtering applied first, then Subject-wise filtering
% request.minsubs_cell=5;
% log.filter_high=5;
% log.filter_low=-5;
% log.meancentre=0;
% % request.flex_mintrials4inclusion=0; % Trialtype models - dont need filtering
% % request.minsubs_cell=0;
% % log.filter_high=999;
% % log.filter_low=-999;
% % log.meancentre=1;
% 
% 
% % Mean centre + High/Low-pass filter + Other optional transforms (e.g zscore)
% for s=1:log.n_subjs
%     for r=1:log.n_rois
%         
%         % High/Low-pass filter (individual cells excluded only)
%         for e=1:log.nlevels
%             for n=1:log.nlevels
%                 for c=1:log.nchoices
%                     
%                     % cF
%                     if d_roimatrix_cF(e,n,c,r,s)>log.filter_high || d_roimatrix_cF(e,n,c,r,s)<log.filter_low
%                         d_roimatrix_cF(e,n,c,r,s)=nan;
%                     end
%                     
%                     % ct
%                     if d_roimatrix_ct(e,n,c,r,s)>log.filter_high || d_roimatrix_ct(e,n,c,r,s)<log.filter_low
%                         d_roimatrix_ct(e,n,c,r,s)=nan;
%                     end
%                     
%                 end
%             end
%         end
%         
%         if log.meancentre
%             wsr.mean=mean2(nanmean(nanmean(([d_roimatrix_cF(:,:,:,r,s); d_roimatrix_ct(:,:,:,r,s)]))));
%             d_roimatrix_cF(:,:,:,r,s)=d_roimatrix_cF(:,:,:,r,s)-wsr.mean;
%             d_roimatrix_ct(:,:,:,r,s)=d_roimatrix_ct(:,:,:,r,s)-wsr.mean;
%         end
%     end
% end
% 

% % Subject exclusions (data replaced by nans)
% request.excludesubs={}; % 'p04_MP'; 'p41_AC'};
% for s=1:length(request.excludesubs)
%     if sum(strcmp(log.subjects, request.excludesubs{s}))~=1; error(['Invalid subject marked for exclusion (' request.excludesubs{s} ')']); end
%     d_roimatrix_cF(:,:,:,:, find(strcmp(log.subjects, request.excludesubs{s})))=nan;
%     d_roimatrix_ct(:,:,:,:, find(strcmp(log.subjects, request.excludesubs{s})))=nan;
% end
% 
% % Flexible models: no. of trials in each cell?
% if isempty(request.flex_mintrials4inclusion)==0 && strcmp(log.FLmodel(1), 'f')==1
%     % NTrials data to match s_roimatrix and s_roivect below
%     %  s_ntrialmatrix {subject}  = Task x Choice x EnvThreat x NTokens 
%     %  s_ntrialvect= matrix of subject ntrials in horizontal vector  (row=subject, cols = task-choice-envthreat-ntokens)
%     s_ntrialmatrix=cell(log.n_subjs+1, 1);  s_ntrialmatrix{log.n_subjs+1}=zeros(2, log.nchoices, log.nlevels, log.nlevels); s_ntrialvect=nan(log.n_subjs, 2*log.nchoices*(log.nlevels*log.nlevels+1)+ 4);    
%     w.beh=load([where.modellingdata filesep request.beh_filename '.mat']);
%     d_beh=w.beh.subjdata;  % Col 2=cF, col 3=ct, col 4=both, col 5=
%     w.behcol=w.beh.details.col;
%     w.behcol.Design_EnvThreatLevel=size(d_beh{1,2},2)+1;
%     w.behcol.Design_NTokensLevel=size(d_beh{1,2},2)+2;
%     
%     for s=1:log.n_subjs
%         ws.cF=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),2}, [w.behcol.EnvThreat w.behcol.NTokens]);
%         ws.ct=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),3}, [w.behcol.EnvThreat w.behcol.NTokens]);
%         
%         for o2=1:1
%         % Mark design levels (EnvThreat, NTokens)
%         if isempty(strfind(log.FLmodel, 'Chunk'))==1
%             ws.cF(:, w.behcol.Design_EnvThreatLevel) = ws.cF(:, w.behcol.EnvThreat)*6;
%             ws.cF(:, w.behcol.Design_NTokensLevel) = ws.cF(:, w.behcol.NTokens)/2;
%             ws.ct(:, w.behcol.Design_EnvThreatLevel) = ws.ct(:, w.behcol.EnvThreat)*6;
%             ws.ct(:, w.behcol.Design_NTokensLevel) = ws.ct(:, w.behcol.NTokens)/2;
%         else
%             ws.cF(:, w.behcol.Design_EnvThreatLevel) = 10+ws.cF(:, w.behcol.EnvThreat)*6;
%             ws.cF(:, w.behcol.Design_NTokensLevel) = 10+ws.cF(:, w.behcol.NTokens)/2;
%             ws.ct(:, w.behcol.Design_EnvThreatLevel) = 10+ws.ct(:, w.behcol.EnvThreat)*6;
%             ws.ct(:, w.behcol.Design_NTokensLevel) = 10+ws.ct(:, w.behcol.NTokens)/2;
%             ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==11  | ws.cF(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
%             ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==13 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
%             ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==15 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
%             ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==11 |  ws.cF(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
%             ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==13 |  ws.cF(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
%             ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==15 |  ws.cF(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
%             ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==11  | ws.ct(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
%             ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==13 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
%             ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==15 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
%             ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==11 |  ws.ct(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
%             ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==13 |  ws.ct(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
%             ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==15 |  ws.ct(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
%         end
%         end
%         
%         % Scroll through choice x envthreat x ntokens, nan-out trials with insufficient n trials
%         for ch=1:log.nchoices
%             for ee=1:log.nlevels
%                 e=log.nlevels+1-ee;
%                 for n=1:log.nlevels
%                     
%                     % No. of trials
%                     s_ntrialmatrix{s}(1, ch, e, n)= sum(ws.cF(ws.cF(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.cF(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
%                     s_ntrialmatrix{s}(2, ch, e, n)= sum(ws.ct(ws.ct(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.ct(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
%                     
%                     % Nan-out cF and cts cells with insufficient trials
%                     if s_ntrialmatrix{s}(1, ch, e, n)<request.flex_mintrials4inclusion;
%                         d_roimatrix_cF(e,n, ch, :, s)=nan;
%                     end
%                     if s_ntrialmatrix{s}(2, ch, e, n)<request.flex_mintrials4inclusion;
%                         d_roimatrix_cF(e,n, ch, :, s)=nan;
%                     end
%                     
%                 end
%             end
%         end
%         s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}+s_ntrialmatrix{s};
%     end
%     s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}/log.n_subjs;
%     
%     % Assemble vector.  Conversion of EnvThreat x NTokens matrix to a vector
%     %   Vector order: EnvThreat (descending threat) followed by NTokens (ascending N)
%     for s=1:log.n_subjs+1
%         k=1;
%         for t=1:2
%             for ch=1:log.nchoices
%                 wv=(squeeze(s_ntrialmatrix{s}(t,ch,:,:)))';
%                 s_ntrialvect(s, k:k+length(wv(:))-1)=(wv(:))';
%                 k=k+length(wv(:))+1;
%             end
%             k=k+2;
%         end
%     end
% end

end
% Preprocessing? Mean Centre for subject?

% (a) ROI selection/re-ordering
if isfield(log.requestrois, log.batteryname)==1  % If battery, load requested ROI order
    eval(['request.ROIs=log.requestrois.' log.batteryname ';'])
else request.ROIs=[];
end
% request.ROIs={
    %         % [ANAT]  Core only
    %         'BA10_L';  'BA10_R';  'BA46_L';'BA46_R';  'Cingulum_Ant_L';  'Cingulum_Ant_R';               % Cortex
    %         'NAccCore_L';  'NAccCore_R'; 'NAccShell_L'; 'NAccShell_R'; 'NAcc_L'; 'NAcc_R';    % Striatum
    %         'HPC_aL';  'HPC_aR'; 'Amygdala_L';  'Amygdala_R'; 'Insula_L';  'Insula_R'; 'SNVTA_L'; 'SNVTA_R';  % Exploratory/various
    %         % [ANAT] Non core
    %         'Hippocampus_L';  'Hippocampus_R';  'Striatum_L'; 'Striatum_R'; 'Caudate_L';  'Caudate_R'; 'Putamen_L'; 'Putamen_R';
% };  % Comemnt out entirely if not requesting  particulars

for o1=1:1  % Implement ROI requests + subject selection
    if isempty(request.ROIs)==0
        disp('Re-ordering ROI. See log.orig for details of original orders. --- ');
        log.orig.roi_files=log.roi_files; log.orig.rois=log.rois; log.orig.n_rois=log.n_rois;
        o.d_event=d_event;
        
        % New data variables
        log.n_rois=length(request.ROIs); log.rois=request.ROIs; log.roi_files='See log.orig.roi_files (out of order)';
        d_event=nan(size(o.d_event));
        for r=1:log.n_rois  % Re-read variables
            if isempty(strcmp(log.orig.rois,  request.ROIs{r}))==0
                wrr.orign=find(strcmp(log.orig.rois,   request.ROIs{r}));
                d_event(r, :, :)=o.d_event(wrr.orign, :, :);
            else error(['Cannot find requested ROI: ' request.ROIs{r}]);
            end
        end
        disp('ROIs requested:    -----------'); disp(request.ROIs)
        
        % Clear original data (for workspace memory)
        clear('o');
    end
end


end


% error('Done w extraction directly from beta :))')

%% (3) Re-arrange data for GLM
%       d_event = roi x subject x trial (event allocation within design not specified)

% Assign trial-event betas to the trialstats (add as columns to subjdata)
for r=1:log.n_rois
    eval(['col.' log.rois{r} '=size(subjdata{1,2},2)+4+r;'])
end
for r=1:log.n_rois
    eval(['cn=col.' log.rois{r} ';'])
    for s=1:log.n_subjs
        subjdata{s,2}(:, cn)=squeeze(d_event(r,s,:));
    end
end

% d_glm: data for glm (columns in col)
col.Subject=size(subjdata{1,2},2)+1; d_glm=[];
for s=1:log.n_subjs
    subjdata{s,2}(:, col.Subject)=s;
    d_glm=[d_glm;subjdata{s,2}];
end

% d_trialtype: roi x subject x task x envthreat x ntokens x trial-iteration
d_trialtype=nan(log.n_rois, log.n_subjs, 2, 6,6, 18);
for s=1:log.n_subjs
    ws.d=d_glm(d_glm(:, col.Subject)==s,:);
    ws.d=sortrows(ws.d, [col.EnvThreat col.NTokens col.Trialnum]);
    for t=1:2
        for e=1:6
            for n=1:6
                ws.cell=ws.d(find(ws.d(:, col.Task)==t  &  ws.d(:, col.EnvThreat)==e/6  & ws.d(:, col.NTokens)==n*2),:);
                eval(['ws.r=ws.cell(:, col.' log.rois{1}  ' : col.' log.rois{end} ');'])
                for r=1:log.n_rois
                    d_trialtype(r,s, t,e,n,:)=shiftdim(ws.r(:,r), -5);
                end
            end
        end
    end
    ws=[];
end
m_trialtype=mean(d_trialtype, 6); % mean over all trial-iterations
for r=1:log.n_rois; for t=1:2 % Create group-level mean  
        m_trialtype(r, log.n_subjs+1,t,:,:)=squeeze(mean(m_trialtype(r,:,t,:,:)));
end; end
openvar m_trialtype

% Plot trial-type betas (recompile roi x task x envthreat x ntokens from events)
f.figwidth= 600;  f.figheight=1200; f.subplot_VerHorz=[0.005 0.0005];f.fig_BotTop=[0.001 0.02];  f.fig_LeftRight=[0.05 0.05]; f.subplotcols=3; 
for r=22:22% log.n_rois
    f.f{r}=figure('Name', ['[' num2str(r) '  ' log.rois{r} ']  Task x EnvThreat x NTokens'], 'NumberTitle', 'off', 'Position',[50,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    for s=1:log.n_subjs+1
        if s==log.n_subjs+1; ss=1; else ss=s+1;  end
        subtightplot(log.n_subjs+1, f.subplotcols,  (ss-1)* f.subplotcols+1, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        if s==log.n_subjs+1; text(1, 0.6, 'Group'); 
            
            text(1, 1, log.rois{r},'FontSize', 16,'FontWeight','bold'); 
        else text(1, 0.6, log.subjects{s});  end; axis off
        
        % cF
        subtightplot(log.n_subjs+1, f.subplotcols,  (ss-1)* f.subplotcols+2, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);        
        imagesc((fliplr(squeeze(m_trialtype(r,s,1,:,:))'))'); axis off; axis square; colorbar
        if ss==1; title('cF'); end
        
        % cF
        subtightplot(log.n_subjs+1, f.subplotcols,  (ss-1)* f.subplotcols+3, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);        
        imagesc((fliplr(squeeze(m_trialtype(r,s,2,:,:))'))'); axis off; axis square; colorbar
        if ss==1; title('ct'); end
    end
    
    
end

% Mark new variables (data or betas)
for o1=1:1 
    
    % Columns
    col.TaskChoice = size(subjdata{1,2},2)   +1;
    col.Accept =size(subjdata{1,2},2)   +2;
    col.Reject=size(subjdata{1,2},2)   +3;
    col.Explore =size(subjdata{1,2},2)   +4;
    
    % Create new variables
    d_glm(:,col.Accept)=d_glm(:,col.Choice)==1; %#ok<SAGROW>
    d_glm(:,col.Reject)=d_glm(:,col.Choice)==2; %#ok<SAGROW>
    d_glm(:,col.Explore)=d_glm(:,col.Choice)==3; %#ok<SAGROW>
    d_glm(:,col.TaskChoice)=(d_glm(:,col.Task)-1).*3 +d_glm(:,col.Choice); %#ok<SAGROW> % 1=cF Accept, 2= cF Reject, 3= cF Explore, 4= ct Accept, 5= ct Bomb, 6= ct Explore
end

%% (3b) Checks

% Collect betas in b (row = sub x trial
eval(['b=d_glm(:, col.' log.rois{1}  ' : col.' log.rois{end} ');'])

% Extreme values?
f.figwidth= 1800;  f.figheight=1000; f.subplot_VerHorz=[0.05 0.055];f.fig_BotTop=[0.001 0.035];  f.fig_LeftRight=[0.05 0.01]; f.subplotcols=5; 
f.f=figure('Name', 'Subject checks', 'NumberTitle', 'off', 'Position',[50,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
for s=1: log.n_subjs  %log.n_rois; % log.n_subjs
    subtightplot(ceil(log.n_subjs/f.subplotcols), f.subplotcols,  s, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
    imagesc(squeeze(d_event(:,s,:)));  title(log.subjects{s}); 
%     imagesc(squeeze(d_event(s,:,:))); title(log.rois{s}); 
    colorbar
%     caxis([-10 10]);
%     axis off; axis square; caxis(log.betarange{r,1});
end
    
%% (4) Manual GLM (manually specify IVs)

% Restrict to certain datas?
task=1; % 1=cF, 2=ct, 3=cFct;
switch task
    case 1;  dd=d_glm(d_glm(:,col.Task)==1,:); disp('cF only!!!'); taskname='cF';
    case 2; dd=d_glm(d_glm(:,col.Task)==2,:); disp('ct only!!!'); taskname='ct';
    case 3; dd=d_glm; taskname='cFct';
end

% Settings ####
request.all_ivs={'Subject';'Choice';
                    'EnvThreat';'NTokens';'pLoss';'Entropy';'VExplore';'EntropyNTok';
                    'EV';'OutcomeMagnitude';'OutcomeMean';'OutcomeVariance';
                    'Trialtype';'RT1';'Task'};
request.rl={'EnvThreat';'NTokens';'pLoss';'Entropy';'VExplore';'EntropyNTok';};
request.outcome={'EV';'OutcomeMagnitude';'OutcomeMean';'OutcomeVariance';};
%
request.ivs=request.all_ivs;
% request.ivs=request.rl;
d_ivs=[]; for i=1:length(request.ivs) % Construct IV columns
    eval(['d_ivs=[d_ivs   dd(:, col.'  request.ivs{i}  ')];'])
end


% Run GLM!! 
r_glmr= [[{' '} {'Constant'} request.ivs'];  [log.rois cell(log.n_rois, length(request.ivs)+1)   ]]; r_glmsr=r_glmr; r_glmp=r_glmr;
for r=1:log.n_rois
    eval(['d_dv=dd(:, col.' log.rois{r} ');'])
    [wrr.b wrr.d wrr.stats]=glmfit(d_ivs, d_dv);
    
    % Record
    r_glmr(r+1, 2:end)=num2cell(wrr.b)';
    r_glmp(r+1, 2:end)=num2cell(wrr.stats.p)';
    for i=2:size(r_glmr,2)
        if r_glmp{r+1, i}<0.001;
            r_glmsr{r+1, i}= [num2str(r_glmr{r+1, i}, 3) '  ***'];
        elseif r_glmp{r+1, i}<=0.05;
            r_glmsr{r+1, i}= [num2str(r_glmr{r+1, i}, 3) '  *'];
        elseif r_glmp{r+1, i}<0.1;
            r_glmsr{r+1, i}= [num2str(r_glmr{r+1, i}, 2) '  (t)'];
        end
    end
    
    wr=[];
end
% openvar r_glmsr

%% (5) Stepwise GLM

d_stepglm=[[{' '};log.rois] cell(log.n_rois+1, 1)]; r_stepglm=[[{' '};log.rois] [[{'Constant'} request.ivs' {'Interactions'}];  cell(log.n_rois, length(request.ivs)+2)]]; % all interactions in last column
for r=1: log.n_rois
    eval(['d_dv=dd(:, col.' log.rois{r} ');'])
    w.c=clock; disp(['[ROI # ' num2str(r) ':  ' log.rois{r} ']        (' num2str(w.c(4)) ':'  num2str(w.c(5)) ')'])
    
    % Run and record
    wrr=GeneralizedLinearModel.stepwise(d_ivs,d_dv,request.glmtype,  'VarNames', [request.ivs' log.rois(r)]);     % do NOT open this variable in workspace - will crash!
    wr.sig_ivs=wrr.PredictorNames;
    wr.sig_coeffs=wrr.CoefficientNames;
    wr.coeff=wrr.Coefficients;
    wr.ivs_pool=wrr.VariableNames(1:end-1);
    wr.modfit_LL=wrr.LogLikelihood;
    wr.modfit_predict=wrr.predict ;
    wr.modfit_formula=wrr.Formula;
    wr.modfit_residuals=wrr.Residuals;
    wr.modfit_criteria=wrr.ModelCriterion;
    wr.modfit_sse=wrr.SSE;
    wr.modfit_fullvsconstantp=wrr.coefTest;
    d_stepglm{r+1,2}=wr;
    
    % Log in significance table
    for p=1: length(wr.sig_coeffs')
        if strcmp(wr.sig_coeffs{p}, '(Intercept)')==1;  c=find(strcmp(r_stepglm(1,:) , 'Constant'));
        else c=find(strcmp(r_stepglm(1,:) ,wr.sig_coeffs{p}));
        end
        
        if isempty(c)==0 % INTERACTION/ Cannot find specified condition
            if wr.coeff{p,end}<0.001;  r_stepglm{r+1,c}=[num2str(wr.coeff{p,1},2) ' ***'];
            elseif wr.coeff{p,end}<0.01; r_stepglm{r+1,p+1}=[num2str(wr.coeff{p,1},2) '  **'];
            elseif wr.coeff{p,end}<=0.05; r_stepglm{r+1,p+1}=[num2str(wr.coeff{p,1},2) '  *'];
            elseif wr.coeff{p,end}<0.1; r_stepglm{r+1,p+1}=[num2str(wr.coeff{p,1},2) '  (t)'];
            end
        else
            c=find(strcmp(r_stepglm(1,:) ,'Interactions'));
            r_stepglm{r+1, c}(size(r_stepglm{r+1, c},1)+1, 1:3)={wr.sig_coeffs{p}   wr.coeff{p,1}  wr.coeff{p,end}}; % Name, beta, p
        end
    end
        
    %
    wr=[]; wrr=[];
end
sound(randn(4000, 1), 4192) % Alert done
openvar r_stepglm


%%  Save

cd(log.roi_fol)
input('Save completed GLM results?');
save(['Glm_results (' date ') '  taskname ])


