% Extract ROI betas using SPM functions (subject-level) for FLEXIBLE models only
%       Flex models: Task x Choice x EnvThreat x NTokens i.e. empty cells for some subjects
%       Trialtype models: Task x EnvThreat x NTokens 
%       Extract betas for all available contrasts (or specified numbers). Identity of each contrast varies from subject to subject.
clear all;close all hidden; clc
% clear all;clc

% Request
log.specificsubjects={};
request.loadextractedbetas=[]; % 'find';    % empty to extract, 'find' to look for file, specify date otherwise
request.plotacrossdesign=0;
request.plotsubjecchecks=0;
%
request.behfilename='All data (09-May-2014)'; % From modelling folder
log.AntsType='_Basic';

% log.FLmodel='t1_1_Trialtype';  % ##### WHICH FL MODEL #########
log.FLmodel='t2_1_TrialtypeNc';

for o1=1:1 % General settings and specifications 
    
    % Folders
    w.w=pwd; if strcmp(w.w(1), '/')==1; where.where='/Users/EleanorL/Dropbox/SANDISK/5 Explore fMRI'; where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data';  where.data_brain=[where.experiment_folder filesep '1 Brain data'];  where.modellingdata=     '/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs';   where.modfol='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
    else where.where='D:\Dropbox\SANDISK\5 Explore fMRI'; where.experiment_folder='C:\Users\eloh\Desktop\2 [Explore]';  where.modfol='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
%         where.data_brain=[where.experiment_folder filesep '1 Brain data']; 
        where.data_brain='I:\1 Explore fMRI'; addpath(where.where);  where.modellingdata='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs';   end;  addpath(where.where); addpath(where.modfol)
    if isdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType])==0; mkdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType]); end
    log.SLfol=[where.experiment_folder filesep '2 Second level results s4Ants' filesep log.FLmodel log.AntsType filesep]; if isdir(log.SLfol)==0; mkdir(log.SLfol); end
    log.FLfol=['m_' log.FLmodel log.AntsType  ' Contrasted'];
    
    % Subjects
    log.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=log.w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    request.specificsubjects=log.specificsubjects;
    
    %
    if isempty(strfind(log.FLmodel, 'Chunk'))==1; log.nlevels=6;  log.trialchar={'t';[]};
    else log.nlevels=3;  log.trialchar={'e';'n'};
    end
    if strcmp(log.FLmodel(1), 'f')==1; log.nchoices=3; else log.nchoices=1; end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end; disp(' ');
    disp(['             First level model:              ' log.FLmodel])
    disp(['             No. of EnvThreat/NTokens levels:       ' num2str(log.nlevels)]); 
    disp(['             Cells divided by how many choices:     ' num2str(log.nchoices)]); disp(' ')
    disp(['             FL folder:      ' log.FLfol])
    disp(['             SL folder:      ' log.SLfol]); disp(' ');
%     input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% (1) ROI settings + Preprocessing 
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per typical plot)

% Selection settings ############
log.meancentre=1;

% Where are the ROI images/beta file  ############
% log.roi_fol=[log.SLfol 'c13g battery'];
% log.roi_fol=[log.SLfol 'Anat HPC Amyg'];
log.roi_fol=[log.SLfol 'c14 HPC'];

% ROI details (List of ROIs + requested order)
for o=1:1  % Battery details
    for o1=1:1 % From before beh modelling done (April 2015)
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
            'BA10_L';  'BA10_R'; 'BA46_L';'BA46_R'; 'Cingulum_Ant_L';  'Cingulum_Ant_R';               % Cortex
            'Striatum_L'; 'Striatum_R'; 'Caudate_L';  'Caudate_R'; 'Putamen_L'; 'Putamen_R'; % Striatum
            'NAccCore_L';  'NAccCore_R'; 'NAccShell_L'; 'NAccShell_R'; 'NAcc_L'; 'NAcc_R';
            'HPC_aL';  'HPC_aR';   'Hippocampus_L';  'Hippocampus_R';    % HPC/Amygdala
            'Amygdala_L';  'Amygdala_R';
            'Insula_L';  'Insula_R'; 'SNVTA_L'; 'SNVTA_R';  % Exploratory/various
            };
        
    end
    
    %     log.roilist.c13gbattery='c13g battery';
    log.roilist.c13gbattery={'DLPFC_L';'DLPFC_R';'HPC_L_sc';'HPC_L_stc';'HPC_R_satc';'HPC_R_sc';'Parietal';'OFC_L';'Striatum_R';};
    log.roilist.anatHPCamy={'Amyg_L';'Amyg_R';'HPC_aL';'HPC_aR';};
    log.roilist.c14HPC={
        'HPC_aL_c';
        'HPC_aL_stc';
        'HPC_aL_tc'; 
        'HPC_aR_c';
        'HPC_aR_satc';
        'HPC_aR_stc';
        'HPC_aR_tc'; };
end
% log.batteryname='c13gbattery';
% log.batteryname='anatHPCamy';
log.batteryname='c14HPC';
eval(['log.ROInames=log.roilist.' log.batteryname ';'])
% log.ROInames={'BA10_L';'BA10_L'};

for o=1:1 % Load or Extract
    
    if isempty(request.loadextractedbetas)==1
        % Display + Load
        log.roi_files=cellstr(spm_select('List', log.roi_fol, '.nii$')); log.n_rois=length(log.roi_files);
        if isfield(log, 'batteryname')==1; disp(['######    Battery name = ' log.batteryname '  ######']); end, disp('------  ROI files       ---------------      ROI short names ------'); disp([log.roi_files   log.ROInames]); input('OK?    ');
        d_roimasks=cell(length(log.roi_files),2);
        for r=1:log.n_rois
            d_roimasks{r,1}=log.ROInames{r};
            d_roimasks{r,2}=spm_read_vols(spm_vol([log.roi_fol filesep log.roi_files{r}]));
        end
    elseif  strcmp(request.loadextractedbetas, 'find')==1
        request.loadextractedbetas = spm_select('List', log.roi_fol, 'Extracted beta matrix ');
        if size(request.loadextractedbetas ,1)~=1; error(['Requested FIND extracted beta file, but none/>1 found (' log.roi_fol ')']); end
        request.loadextractedbetas = request.loadextractedbetas(strfind(request.loadextractedbetas,  '('): end-4);
    end
    
    % (1) Full flexible extraction: Load all images and extract (conditionally per subject)
    %  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
    %                       EnvThreat in reverse order (i.e. as per typical plot)
    
    if isempty(request.loadextractedbetas)==1
        d_roimatrix_cF= nan(log.nlevels,log.nlevels,log.nchoices, length(log.roi_files), log.n_subjs); d_roimatrix_ct=d_roimatrix_cF; l_roimatrix_cF=d_roimatrix_cF; l_roimatrix_ct=d_roimatrix_cF; % log: each subject ok or not?
        if log.nchoices==3; log.choicecF={'Accept_';'Reject_';'Explore_'};  log.choicect={'NoBomb_';'Bomb_';'Explore_'};
        else log.choicecF={[]};  log.choicect={[]};
        end
        for s=1:log.n_subjs
            ws.c=clock;  disp(['Subject ' num2str(s) '  (' log.subjects{s} ')            [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  --------------'])
            ws.FLfol=[where.data_brain filesep log.subjects{s} filesep '2 First level s4Ants'  filesep log.FLfol  filesep];
            ws.s=load([ws.FLfol 'SPM.mat']);
            ws.mask=spm_read_vols(spm_vol([ws.FLfol 'mask.img']));
            
            % Load all contrasts for this subject
            ws.con=cell(size(ws.s.SPM.xCon,2),3);
            for c=1:size(ws.s.SPM.xCon,2)
                ws.con{c,1}=ws.s.SPM.xCon(c).name;
                ws.con{c,2}=ws.s.SPM.xCon(c).Vcon.fname;
                ws.con{c,3}=spm_read_vols(spm_vol([ws.FLfol  ws.con{c,2}]));
                %             % Fake contrast inputs;
                %             disp('Contrast extraction turned off! Fake cons!');
                %             ws.con{c,3}=nan(79,95,68);
                %             ws.con{c,3}(1:20,1:20, 1:20)=rand(20,20,20);
                %             ws.con{c,3}(:,:,:)=rand(size(ws.con{c,3}));
            end
            
            % Read betas for all available cF contrast
            disp('        Loading cF')
            for ee=1:log.nlevels
                e=log.nlevels+1-ee;
                for n=1:log.nlevels
                    for ch=1:log.nchoices
                        wc.conname=['cF_' log.choicecF{ch} log.trialchar{1} num2str(ee) '-' log.trialchar{2} num2str(n)];
                        if isempty(find(strcmp(ws.con, wc.conname)))==0
                            wc.con=ws.con{find(strcmp(ws.con(:,1),  wc.conname)), 3};
                            for r=1: log.n_rois
                                % Read all rois for this contrast. If crashes here, check that ROI image and contrast images have the same dimensions!
                                d_roimatrix_cF(e,n,ch,r,s)= nanmean( wc.con(d_roimasks{r,2}~=0));
                                %
                                
                                l_roimatrix_cF(e,n,ch,r,s)=1;
                            end
                        else
                            l_roimatrix_cF(e,n,ch,:,s)=0;
                        end
                        wc=[];
                    end
                end
            end
            
            % Read betas for all available ct contrast
            disp('        Loading ct')
            for ee=1:log.nlevels
                e=log.nlevels+1-ee;
                for n=1:log.nlevels
                    for ch=1:log.nchoices
                        wc.conname=['ct_' log.choicect{ch} log.trialchar{1} num2str(ee) '-' log.trialchar{2} num2str(n)];
                        if isempty(find(strcmp(ws.con, wc.conname)))==0
                            wc.con=ws.con{find(strcmp(ws.con(:,1),  wc.conname)), 3};
                            for r=1: log.n_rois % Read all rois for this contrast
                                d_roimatrix_ct(e,n,ch,r,s)=nanmean( wc.con(d_roimasks{r,2}~=0));
                                %
                                l_roimatrix_ct(e,n,ch,r,s)=1;
                            end
                        else
                            l_roimatrix_ct(e,n,ch,:,s)=0;
                        end
                        wc=[];
                    end
                end
            end
            
            %
            ws=[];
        end
        save([log.roi_fol filesep 'Extracted beta matrix (' date ')'], 'd_roimatrix_cF', 'd_roimatrix_ct', 'l_roimatrix_cF', 'l_roimatrix_ct', 'log')
        try % Notify researcher
            f_sendemail('kurzlich', strcat(['Analysis script is complete: extracting betas from all fully flexible contrasts (' mfilename ')']), ' ',1);
        end
        disp('======================================================='); w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
        if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
        disp(['             First level model:      ' log.FLmodel]); disp('=======================================================')
    else disp(['Loading dataset previously extracted ' request.loadextractedbetas ]); load([log.roi_fol filesep 'Extracted beta matrix ' request.loadextractedbetas '.mat'])
        if isempty(request.specificsubjects)==0 % Implement subject selections requested in script
            w.dcf=d_roimatrix_cF; w.dct=d_roimatrix_ct;  w.lcf=d_roimatrix_cF; w.lct=d_roimatrix_ct; d_roimatrix_cF= nan(log.nlevels,log.nlevels,log.nchoices, length(log.roi_files), length(request.specificsubjects)); d_roimatrix_ct=d_roimatrix_cF; l_roimatrix_cF=d_roimatrix_cF; l_roimatrix_ct=d_roimatrix_cF;
            for s=1:length(request.specificsubjects)
                if sum(strcmp(log.subjects, request.specificsubjects{s}))~=1; error(['Invalid subject requested: ' request.specificsubjects{s}]); end
                d_roimatrix_cF(:,:,:,:, s)= w.dcf(:,:,:,:, find(strcmp(log.subjects, request.specificsubjects{s})));
                d_roimatrix_ct(:,:,:,:, s)= w.dct(:,:,:,:, find(strcmp(log.subjects, request.specificsubjects{s})));
                l_roimatrix_cF(:,:,:,:, s)= w.lcf(:,:,:,:, find(strcmp(log.subjects, request.specificsubjects{s})));
                l_roimatrix_ct(:,:,:,:, s)= w.lct(:,:,:,:, find(strcmp(log.subjects, request.specificsubjects{s})));
            end
            log.subjects=request.specificsubjects; log.n_subjs=length(log.subjects); clear('w.dcf','w.dct','w.lcf','w.lct')
        end
    end
end
if isfield(log.requestrois, log.batteryname)==1;  eval(['request.ROIs=log.requestrois.' log.batteryname ';'])   else request.ROIs=[]; end
request.ROIs={
%     'HPC_L_stc';'HPC_R_satc';  % c13 battery ################
%     'HPC_L_sc';'HPC_R_sc';
%     'Striatum_R';'DLPFC_L';'DLPFC_R';
'HPC_aL';'HPC_aR';'Amyg_L';'Amyg_R';       % Anat ########
'HPC_aL_c';   'HPC_aR_c';  %  c14 hippocampus ########
'HPC_aL_tc';   'HPC_aR_tc'; 
'HPC_aL_stc';   'HPC_aR_stc';  'HPC_aR_satc';
};  % Comemnt out entirely if not requesting  particulars
for o=1:1 % Preprocessing
    
    request.flex_mintrials4inclusion=[]; % empty to omit trialwise filtering. Trialwise filtering applied first, then Subject-wise filtering
    request.minsubs_cell=1;
    log.filter_high=500;  % Set to preposterous to omit
    log.filter_low=-500;
    for o1=1:1  % Implement ROI requests
        if isempty(request.ROIs)==0
            disp('Re-ordering ROI. See log.orig for details of original orders. --- ');
            log.orig.roi_files=log.roi_files;
            log.orig.ROInames=log.ROInames;
            log.orig.n_rois=log.n_rois;
            o.d_cF=d_roimatrix_cF;
            o.d_ct=d_roimatrix_ct;
            o.l_cF=l_roimatrix_cF;
            o.l_ct=l_roimatrix_ct;
            
            % New data variables
            log.n_rois=length(request.ROIs);
            log.ROInames=request.ROIs;
            log.roi_files='See log.orig.roi_files (out of order)';
            d_roimatrix_cF= nan(log.nlevels,log.nlevels,log.nchoices, log.n_rois, log.n_subjs); d_roimatrix_ct=d_roimatrix_cF; % e x n x c x r x s
            l_roimatrix_cF=d_roimatrix_cF; l_roimatrix_ct=d_roimatrix_cF;
            
            % Re-read variables
            for r=1:log.n_rois
                if isempty(strcmp(log.orig.ROInames,  request.ROIs{r}))==0
                    wr.orign=find(strcmp(log.orig.ROInames,   request.ROIs{r}));
                    d_roimatrix_cF(:,:,:,r,:)=o.d_cF(:,:,:,wr.orign,:);
                    d_roimatrix_ct(:,:,:,r,:)=o.d_ct(:,:,:,wr.orign,:);
                    l_roimatrix_cF(:,:,:,r,:)=o.l_cF(:,:,:,wr.orign,:);
                    l_roimatrix_ct(:,:,:,r,:)=o.l_ct(:,:,:,wr.orign,:);
                else
                    error(['Cannot find requested ROI: ' request.ROIs{r}]);
                end
            end
            disp('ROIs requested:    -----------'); disp(request.ROIs)
            
            % Clear original data (for workspace memory)
            clear('o');
        end
    end
    
    % Mean centre + High/Low-pass filter + Other optional transforms (e.g zscore)
    for s=1:log.n_subjs
        for r=1:log.n_rois
            
            % High/Low-pass filter (individual cells excluded only)
            if s==1 && r==1,  disp(['High/low-pass fltering between: ' num2str(log.filter_low) ' to ' num2str(log.filter_high) ]), end
            for e=1:log.nlevels
                for n=1:log.nlevels
                    for c=1:log.nchoices
                        
                        % cF
                        if d_roimatrix_cF(e,n,c,r,s)>log.filter_high || d_roimatrix_cF(e,n,c,r,s)<log.filter_low
                            d_roimatrix_cF(e,n,c,r,s)=nan;
                        end
                        
                        % ct
                        if d_roimatrix_ct(e,n,c,r,s)>log.filter_high || d_roimatrix_ct(e,n,c,r,s)<log.filter_low
                            d_roimatrix_ct(e,n,c,r,s)=nan;
                        end
                        
                    end
                end
            end
            
            if log.meancentre
                wsr.mean=mean2(nanmean(nanmean(([d_roimatrix_cF(:,:,:,r,s); d_roimatrix_ct(:,:,:,r,s)]))));
                d_roimatrix_cF(:,:,:,r,s)=d_roimatrix_cF(:,:,:,r,s)-wsr.mean;
                d_roimatrix_ct(:,:,:,r,s)=d_roimatrix_ct(:,:,:,r,s)-wsr.mean;
            end
        end
    end
    
    % Subject exclusions (data replaced by nans)
    request.excludesubs={}; % 'p04_MP'; 'p41_AC'};
    for s=1:length(request.excludesubs)
        if sum(strcmp(log.subjects, request.excludesubs{s}))~=1; error(['Invalid subject marked for exclusion (' request.excludesubs{s} ')']); end
        d_roimatrix_cF(:,:,:,:, find(strcmp(log.subjects, request.excludesubs{s})))=nan;
        d_roimatrix_ct(:,:,:,:, find(strcmp(log.subjects, request.excludesubs{s})))=nan;
    end
    
    % Flexible models: no. of trials in each cell?
    if isempty(request.flex_mintrials4inclusion)==0 && strcmp(log.FLmodel(1), 'f')==1
        for o=1:1 % Flex
            if s==1, input('FILTERING by no. trials etc. Continue?');  end
            
            % NTrials data to match s_roimatrix and s_roivect below
            %  s_ntrialmatrix {subject}  = Task x Choice x EnvThreat x NTokens
            %  s_ntrialvect= matrix of subject ntrials in horizontal vector  (row=subject, cols = task-choice-envthreat-ntokens)
            s_ntrialmatrix=cell(log.n_subjs+1, 1);  s_ntrialmatrix{log.n_subjs+1}=zeros(2, log.nchoices, log.nlevels, log.nlevels); s_ntrialvect=nan(log.n_subjs, 2*log.nchoices*(log.nlevels*log.nlevels+1)+ 4);
            w.beh=load([where.modellingdata filesep request.behfilename '.mat']);
            d_beh=w.beh.subjdata;  % Col 2=cF, col 3=ct, col 4=both, col 5=
            w.behcol=w.beh.details.col;
            w.behcol.Design_EnvThreatLevel=size(d_beh{1,2},2)+1;   w.behcol.Design_NTokensLevel=size(d_beh{1,2},2)+2;
            
            for s=1:log.n_subjs
                ws.cF=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),2}, [w.behcol.EnvThreat w.behcol.NTokens]);
                ws.ct=sortrows(d_beh{find(strcmp( d_beh(:,1),  log.subjects{s})),3}, [w.behcol.EnvThreat w.behcol.NTokens]);
                
                for o2=1:1
                    % Mark design levels (EnvThreat, NTokens)
                    if isempty(strfind(log.FLmodel, 'Chunk'))==1
                        ws.cF(:, w.behcol.Design_EnvThreatLevel) = ws.cF(:, w.behcol.EnvThreat)*6;
                        ws.cF(:, w.behcol.Design_NTokensLevel) = ws.cF(:, w.behcol.NTokens)/2;
                        ws.ct(:, w.behcol.Design_EnvThreatLevel) = ws.ct(:, w.behcol.EnvThreat)*6;
                        ws.ct(:, w.behcol.Design_NTokensLevel) = ws.ct(:, w.behcol.NTokens)/2;
                    else
                        ws.cF(:, w.behcol.Design_EnvThreatLevel) = 10+ws.cF(:, w.behcol.EnvThreat)*6;
                        ws.cF(:, w.behcol.Design_NTokensLevel) = 10+ws.cF(:, w.behcol.NTokens)/2;
                        ws.ct(:, w.behcol.Design_EnvThreatLevel) = 10+ws.ct(:, w.behcol.EnvThreat)*6;
                        ws.ct(:, w.behcol.Design_NTokensLevel) = 10+ws.ct(:, w.behcol.NTokens)/2;
                        ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==11  | ws.cF(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
                        ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==13 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
                        ws.cF(ws.cF(:, w.behcol.Design_EnvThreatLevel)==15 | ws.cF(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
                        ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==11 |  ws.cF(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
                        ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==13 |  ws.cF(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
                        ws.cF(ws.cF(:, w.behcol.Design_NTokensLevel)==15 |  ws.cF(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
                        ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==11  | ws.ct(:, w.behcol.Design_EnvThreatLevel)==12, w.behcol.Design_EnvThreatLevel)=1;
                        ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==13 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==14, w.behcol.Design_EnvThreatLevel)=2;
                        ws.ct(ws.ct(:, w.behcol.Design_EnvThreatLevel)==15 | ws.ct(:, w.behcol.Design_EnvThreatLevel)==16, w.behcol.Design_EnvThreatLevel)=3;
                        ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==11 |  ws.ct(:, w.behcol.Design_NTokensLevel)==12, w.behcol.Design_NTokensLevel)=1;
                        ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==13 |  ws.ct(:, w.behcol.Design_NTokensLevel)==14, w.behcol.Design_NTokensLevel)=2;
                        ws.ct(ws.ct(:, w.behcol.Design_NTokensLevel)==15 |  ws.ct(:, w.behcol.Design_NTokensLevel)==16, w.behcol.Design_NTokensLevel)=3;
                    end
                end
                
                % Scroll through choice x envthreat x ntokens, nan-out trials with insufficient n trials
                for ch=1:log.nchoices
                    for ee=1:log.nlevels
                        e=log.nlevels+1-ee;
                        for n=1:log.nlevels
                            
                            % No. of trials
                            s_ntrialmatrix{s}(1, ch, e, n)= sum(ws.cF(ws.cF(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.cF(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
                            s_ntrialmatrix{s}(2, ch, e, n)= sum(ws.ct(ws.ct(:,w.behcol.Design_EnvThreatLevel)==ee &  ws.ct(:,w.behcol.Design_NTokensLevel)==n , w.behcol.Choice)==ch);
                            
                            % Nan-out cF and cts cells with insufficient trials
                            if s_ntrialmatrix{s}(1, ch, e, n)<request.flex_mintrials4inclusion;
                                d_roimatrix_cF(e,n, ch, :, s)=nan;
                            end
                            if s_ntrialmatrix{s}(2, ch, e, n)<request.flex_mintrials4inclusion;
                                d_roimatrix_cF(e,n, ch, :, s)=nan;
                            end
                            
                        end
                    end
                end
                s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}+s_ntrialmatrix{s};
            end
            s_ntrialmatrix{log.n_subjs+1}=s_ntrialmatrix{log.n_subjs+1}/log.n_subjs;
            
            % Assemble vector.  Conversion of EnvThreat x NTokens matrix to a vector
            %   Vector order: EnvThreat (descending threat) followed by NTokens (ascending N)
            for s=1:log.n_subjs+1
                k=1;
                for t=1:2
                    for ch=1:log.nchoices
                        wv=(squeeze(s_ntrialmatrix{s}(t,ch,:,:)))';
                        s_ntrialvect(s, k:k+length(wv(:))-1)=(wv(:))';
                        k=k+length(wv(:))+1;
                    end
                    k=k+2;
                end
            end
        end
    end
    
    % Check betas for nans
    docheck=1;
    if docheck
        for s=1:log.n_subjs, k=1;
            for c=1:log.nchoices
                for r=1:log.n_rois
                    ws.r(k)=sum(sum(isnan(d_roimatrix_cF(:,:,1,r,s)))); k=k+1;
                end
            end
            disp([log.subjects{s} '  : '  num2str(sum(ws.r)) ' bad cells'])
            ws=[];
        end
    end
end
for o2=1:1 % Subject betas (sanity check)
    
    %  s_roimatrix {subject, roi}  = Task x Choice x EnvThreat x NTokens
    %  s_roivect{roi}= matrix of subject betas in horizontal vector  (row=subject, cols = task-choice-envthreat-ntokens)
    
    % Assemble subject beta vectors
    s_roimatrix=cell(log.n_subjs, log.n_rois); s_roivect=repmat({nan(log.n_subjs, 2*log.nchoices*(log.nlevels*log.nlevels+1)+ 4)}, log.n_rois,1);
    for s=1:log.n_subjs
        for r=1:log.n_rois
            s_roimatrix{s, r}=nan(2, log.nchoices, log.nlevels, log.nlevels); s_roimatrix{s, r}=nan(2, log.nchoices, log.nlevels, log.nlevels);
            for ch=1:log.nchoices
                s_roimatrix{s, r}(1, ch,:,:)=d_roimatrix_cF(:,:, ch, r, s);
                s_roimatrix{s, r}(2, ch,:,:)=d_roimatrix_ct(:,:, ch, r, s);
            end
            
            % Assemble vector.  Conversion of EnvThreat x NTokens matrix to a vector
            %   Vector order: EnvThreat (descending threat) followed by NTokens (ascending N)
            k=1;
            for t=1:2
                for ch=1:log.nchoices
                    wv=(squeeze(s_roimatrix{s,r}(t,ch,:,:)))';
                    s_roivect{r}(s, k:k+length(wv(:))-1)=(wv(:))';
                    k=k+length(wv(:))+1;
                end
                k=k+2;
            end
        end
    end
    
    % Plot
    if request.plotsubjecchecks
        log.plotsubjbetarange=[]; % -10 10]; % Plotting range only, not actually data filtering
        %
        f.subplotcols=3; k=1; % How many columns?
        f.figwidth= 1400; f.figheight=1000;  f.subplot_VerHorz=[0.05 0.1];
        f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.01];  f.nancol=[0 0 0];
        fs=figure('Name', [log.FLmodel ': Subject betas in each cell'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
        if isempty(log.plotsubjbetarange)==0; disp(['Beta range for subject plot: ' num2str(log.plotsubjbetarange(1)) ' to '  num2str(log.plotsubjbetarange(2)) ]); end
        for r=1:log.n_rois
            subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            if isempty(log.plotsubjbetarange)==1;  imagescnan(s_roivect{r}, 'NaNColor', f.nancol);   else imagescnan(s_roivect{r}, log.plotsubjbetarange, 'NaNColor', f.nancol);  end
            colorbar; set(gca, 'XTick', []); title(log.ROInames{r})
            k=k+1;
        end
        
        if log.nchoices ==3
            % Total no. of trials in each cell for each subject
            subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(s_ntrialvect, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
            title('Total no. of trials in each cell'); k=k+1;
            
            % Mean no. of trials in each cell for each subject
            subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(s_ntrialvect(log.n_subjs+1,:), 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
            title('Mean no. of trials in each cell'); k=k+1;
            
            % Cell inclusion
            subtightplot(ceil((log.n_rois+3)/f.subplotcols), f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(s_ntrialvect>=request.flex_mintrials4inclusion, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
            title(['Cell included or not for each subject (min ntrials=' num2str(request.flex_mintrials4inclusion) ')']); k=k+1;
        end
    end
end

%% (2) Plot betas
% Calculate means (and log no. of subjects, according to minimum)
%       m_roimatrix: EnvThreat x NTokens x Choice x ROI
m_roimatrix_cF=nan(log.nlevels,log.nlevels,log.nchoices, log.n_rois); m_roimatrix_ct=m_roimatrix_cF;  m_roimatrix_comb=m_roimatrix_cF; % e x n x choice x roi
n_roimatrix=nan(log.nlevels,log.nlevels,log.nchoices,2); % e x n x choice x task
for e=1:log.nlevels
    for n=1:log.nlevels
        for c=1:log.nchoices
            for r=1: log.n_rois
                wr.cf=squeeze(d_roimatrix_cF(e,n,c,r,:));
                wr.ct=squeeze(d_roimatrix_ct(e,n,c,r,:));
                wr.cf=wr.cf(isnan(wr.cf)==0);  % Remove nans
                wr.ct=wr.ct(isnan(wr.ct)==0);
                
                % Calculate mean ##### 
                %       Other options for mean-ning: trimmean, z score 
                %       Exclusions applied here are across-subject, not within-subject!
                m_roimatrix_cF(e,n,c,r)= mean(wr.cf);    % Typical mean
                m_roimatrix_ct(e,n,c,r)= mean(wr.ct);
                m_roimatrix_comb(e,n,c,r)=mean([wr.cf; wr.ct]);
                if e==1 && n==1 && c==1 && r==1;  disp('Mean = mean'); end
                
                
                % Apply minimums
                n_roimatrix(e,n,c,1)=length(wr.cf);
                n_roimatrix(e,n,c,2)=length(wr.ct);
                if n_roimatrix(e,n,c,1) <request.minsubs_cell;  m_roimatrix_cF(e,n,c,r)=nan; end % Empty cell (nan) if too few subjects
                if n_roimatrix(e,n,c,2) <request.minsubs_cell;  m_roimatrix_ct(e,n,c,r)=nan; end
                
                %
                wr=[];
                
            end
        end
    end
end

% Plot mean betas
if request.plotacrossdesign==1
    log.betarange=[]; %[-1 1]; %[]; %[0 10];
    if strcmp(log.FLmodel(1), 'f')==1 &&  isempty(log.betarange)==1   % Beta range? (standardize within ROIs; overridden by requests)
        log.meanbrange=cell(log.n_rois, 2);   for r=1:log.n_rois % tweak here to calculate ranges for cF and ct separately
            %             log.meanbrange{r,1}=  [min( min(min(m_roimatrix_cF(:,:, :,r))))    max( max(max(m_roimatrix_cF(:,:, :,r)))) ];
            %             log.meanbrange{r,2}=  [min( min(min(m_roimatrix_ct(:,:, :,r))))     max( max(max(m_roimatrix_ct(:,:, :,r)))) ];
            log.meanbrange{r,1}=   [min( min(min([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))    max( max(max([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))  ];
            log.meanbrange{r,2}=log.meanbrange{r,1};
            disp('Natural range for betas (for means plot): '); disp(['   '  log.ROInames{r} ':    '   num2str(log.meanbrange{r,1}(1)) ' to ' num2str(log.meanbrange{r,1}(2))])
        end
        log.betarange =log.meanbrange;
    elseif strcmp(log.FLmodel(1), 't')==1 &&  isempty(log.betarange)==1
        log.meanbrange=cell(log.n_rois, 2); for r=1:log.n_rois
            log.meanbrange{r,1}=   [min( min(min([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))    max( max(max([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))  ];
            log.meanbrange{r,2}=log.meanbrange{r,1};
            disp('Natural range for betas (for means plot): '); disp(['   '  log.ROInames{r} ':    '   num2str(log.meanbrange{r,1}(1)) ' to ' num2str(log.meanbrange{r,1}(2))])
        end ; log.betarange =log.meanbrange;
    elseif strcmp(log.FLmodel(1), 'f')==1;  log.betarange = repmat({log.betarange}, log.n_rois, 2);
    end
    
    % Plot
    f.figwidth= 1400; f.figheight=log.n_rois*250+150;   f.figheight=1000; f.nancol=[0 0 0]; if log.nchoices==3; f.subplotcols=8; f.plotspace=3;  f.plotrowsadd=2; f.subplot_VerHorz=[0.01 0.005];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.15 0.01];  else f.subplotcols=4; f.plotspace=2; f.plotrowsadd=0;   f.subplot_VerHorz=[0.01 0.005]; f.fig_BotTop=[0.001 0.025]; f.fig_LeftRight=[0.3 0.3];   end 
    fm=figure('Name', [log.FLmodel ': Mean betas for each ROI [Accept, Reject, Explore, NoBomb, Bomb, Explore]   (min. no. subjects per cell = ' num2str(request.minsubs_cell) ')'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    for r=1:log.n_rois % Betas for Roi x task
        % cF
        subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        text(-0.5,0.5, log.ROInames{r},'FontSize', 15); axis off
        for c=1:log.nchoices
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_cF(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.nchoices==1; title('cF');  elseif r==1; title(['cF ' log.choicecF{c}(1:end-1)]); end
            if isempty(log.betarange)==0; caxis(log.betarange{r,1}); end
%             colormap('summer')
        end
        
        % ct
%         if log.nchoices==3
%             subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+f.plotspace,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
%             text(0.3,0.5, log.ROInames{r},'FontSize', 8); axis off
%         end
        for c=1:log.nchoices
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+f.plotspace+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_ct(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.nchoices==1; title('ct');  elseif r==1; title(['ct ' log.choicect{c}(1:end-1)]); end
            if isempty(log.betarange)==0; caxis(log.betarange{r,2}); end
%             colormap('summer')
        end
    end
    for o1=1:1 % No. of subjects (add to figure)
        log.range_nsubs=[8 log.n_subjs];
        log.range_ntrials=[8 14];
        if log.nchoices==3;
            t=1; c=1;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'No. of subjects included' ,'FontSize', 8); axis off
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagesc(n_roimatrix(:,:,c,t), log.range_nsubs);  axis off; axis square; colorbar
            if log.nchoices==3
                t=1; c=2;
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(n_roimatrix(:,:,c,t),log.range_nsubs);  axis off; axis square; colorbar
                t=1; c=3;
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(n_roimatrix(:,:,c,t), log.range_nsubs);  axis off; axis square; colorbar
            end
            
            t=2; c=1;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+5,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'No. of subjects included' ,'FontSize', 8); axis off
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+f.plotspace+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagesc(n_roimatrix(:,:,c,t),log.range_nsubs);  axis off; axis square; colorbar
            if log.nchoices==3
                t=2; c=2;
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+f.plotspace+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(n_roimatrix(:,:,c,t),log.range_nsubs);  axis off; axis square; colorbar
                t=2; c=3;
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+f.plotspace+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(n_roimatrix(:,:,c,t), log.range_nsubs);  axis off; axis square; colorbar
            end
            
            % Mean number of trials in each cell
            t=1;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'Mean no. of trials per cell' ,'FontSize', 8); axis off
            for c=1:log.nchoices
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(squeeze(s_ntrialmatrix{log.n_subjs+1}(t,c,:,:)), log.range_ntrials);  axis off; axis square; colorbar
            end
            t=2;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+5,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'Mean no. of trials per cell' ,'FontSize', 8); axis off
            for c=1:log.nchoices
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+f.plotspace+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(squeeze(s_ntrialmatrix{log.n_subjs+1}(t,c,:,:)), log.range_ntrials);  axis off; axis square; colorbar
            end
        end
    end
    f.figheight=log.n_rois*250+150;   f.figheight=1000; f.nancol=[0 0 0];  % Betas x roi (combine tasks)
    if log.nchoices==3; f.figwidth= 1400; f.subplotcols=4; f.plotspace=3; f.subplot_VerHorz=[0.01 0.005];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.15 0.01]; f.nancol=[0 0 0];
    else f.figwidth= 700; f.subplotcols=2;  f.subplot_VerHorz=[0.01 0.005]; f.fig_BotTop=[0.001 0.025]; f.fig_LeftRight=[0.3 0.3];  
    end 
    fm=figure('Name', [log.FLmodel ': Mean betas for each ROI, combined across tasks    (min. no. subjects per cell = ' num2str(request.minsubs_cell) ')'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);    
    for r=1:log.n_rois
        subtightplot(log.n_rois, f.subplotcols,  (r-1)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        text(0.3,0.5, log.ROInames{r},'FontSize', 15); axis off
        for c=1:log.nchoices
            subtightplot(log.n_rois, f.subplotcols,  (r-1)*f.subplotcols+1+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_comb(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.nchoices==1; title('Combined cF ct');  elseif r==1; title(['cFct ' log.choicecF{c}(1:end-1)]); end
            %             if isempty(log.betarange)==0; caxis(log.betarange{r,1}); end
%             colormap('summer')
        end
    end    
end

for o=1:1 % ARCHIVED analysis  
    
    
    % (4) ANOVA analysis on betas (Task x EnvThreat x NTokens ANOVA)
    %  d_roibeta - reformatted betas from d_roimatrix for analysis with ANOVA
    %                   cell, col 1=roi name, col 2=cF betas, col 3= ct betas, col 4=anova stats table
    %   r_anova - significance table for each roi x comparison
    
    request.anova_tasktrialtype=0;
    if request.anova_tasktrialtype
        d_roibeta=cell(log.n_rois, 3); r_anova=cell(log.n_rois+1, 8); r_anova(1,2:end)={'Task';'Env';'NTok';'Task x Env';'Task x NTok';'Env x NTok';'Task x Env x NTok';};
        
        for r=1:log.n_rois
            wr.cf=[]; wr.ct=[]; r_anova{r+1,1}= log.ROInames{r};
            
            % Compile sample
            for ee=1:6
                e=log.nlevels+1-ee;
                for n=1:6
                    wr.cf=[wr.cf squeeze(d_roimatrix_cF(e, n, 1, r, :))];
                    wr.ct=[wr.ct squeeze(d_roimatrix_ct(e, n, 1, r, :))];
                end
            end
            d_roibeta{r,2}=wr.cf; d_roibeta{r,3}=wr.ct;
            
            %
            wr.anova=teg_repeated_measures_ANOVA([d_roibeta{r,2} d_roibeta{r,3}],  [2 6 6], {'Task'; 'EnvThreat'; 'NTokens'});
            d_roibeta{r,4}=[{'Fac' 'F' 'df1' 'df2' 'p' ' ' ' '}; [cellstr(wr.anova.labels') num2cell(wr.anova.R)]];
            for i=1:7
                if d_roibeta{r,4}{1+i, 5}<0.001
                    r_anova{r+1,i+1}='  ***';
                elseif d_roibeta{r,4}{1+i, 5}<0.05
                    r_anova{r+1,i+1}='  *';
                elseif d_roibeta{r,4}{1+i, 5}<0.1
                    r_anova{r+1,i+1}='t';
                end
            end
            wr=[];
        end
        openvar r_anova
    end
end

%% (3) On-sample stats on cell betas

dostats=1;
if dostats
    r_tcf=nan(6,6,log.n_rois); r_tct=r_tcf; r_tcfmct=r_tcf; ms_cf=r_tcf; ms_ct=r_tcf;   ms_cfmct=r_tcf;  % t stat, m_s: mean for sig diff from 0 only
    psig=0.1;       psig=psig/36;
    for r=1: log.n_rois % Compile
        for e=1:6
            for n=1:6
                % cF
                [h p ci stats]=ttest(squeeze(d_roimatrix_cF(e,n,1, r,:)));
                r_tcf(e,n,r)  = stats.tstat;
                if p<psig
                    ms_cf(e,n,r)  = nanmean(squeeze(d_roimatrix_cF(e,n,1, r,:)));
                end
                
                % ct
                [h p ci stats]=ttest(squeeze(d_roimatrix_ct(e,n,1, r,:)));
                r_tct(e,n,r)  = stats.tstat;
                if p<psig
                    ms_ct(e,n,r)  = nanmean(squeeze(d_roimatrix_ct(e,n,1, r,:)));
                end
                
                
                % cF - ct
                [h p ci stats]=ttest(squeeze(d_roimatrix_cF(e,n,1, r,:)) - squeeze(d_roimatrix_ct(e,n,1, r,:)));
                r_tcfmct(e,n,r)  = stats.tstat;
                tp(e,n,r)  =p;
                if p<psig
                    ms_cfmct(e,n,r)  = nanmean( squeeze(d_roimatrix_cF(e,n,1, r,:)) - squeeze(d_roimatrix_ct(e,n,1, r,:))  );
                end
            end
        end
    end
    
    % Plot
%     d_plot={ms_cf ms_ct  ms_cfmct}; % Mean significant from 0
    d_plot={r_tcf r_tct r_tcfmct }; % T statistic
    d_roinames=request.ROIs;
    %
    k=1; f.subplotcols=4;  f.subplot_VerHorz=[0.005 0.03]; f.fig_BotTop=[0.001 0.025]; f.fig_LeftRight=[0.1 0.1];  f.nancol=[0 0 0];
    f.FontSize=25;
    figure('color','w', 'Position', [0 50 600 1200])
    for r=1:log.n_rois
        wr=  [d_plot{2}(:,:,r) d_plot{1}(:,:,r)];
        wr.range=[min(wr(:)) max(wr(:))];
        
        subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        %     text(0.2,0.5, log.ROInames{r}, 'FontSize', f.FontSize);
        text(-0.2,0.5, d_roinames{r}, 'FontSize', f.FontSize);
        axis off, k=k+1;
        
        % Plot cF
        subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(d_plot{1}(:, :,r), 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        if r==1; title('cF', 'FontSize', f.FontSize);end
        axis square, axis off, caxis(wr.range),k=k+1;
        
        
        % Plot ct
        subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(d_plot{2}(:, :,r), 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        if r==1; title('ct', 'FontSize', f.FontSize); axis square, end
        axis square, axis off, caxis(wr.range), k=k+1;
        
        
        % cF - ct
        subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagescnan(d_plot{3}(:, :,r), 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
        if r==1;  title('cF - ct', 'FontSize', f.FontSize); axis square, end
        axis square, axis off, k=k+1;
        
    end
    
end


%% (4) GLM on betas
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
% d_betas:       d_betas{r,1}{task}: all data for single glm
%                     d_betas{r,2}{subject, task}:subject data for 2-level glms
%                     d_betas{r,3}{task}: betas from subject-specific first level glms (subject = row, col=iv)
% d_glmres{task}{roi+1, iv+1}: t stat if sig
% d_glmresp: p values

% Load design (canonical)
desvars={'EnvThreat';'NTokens';'pLoss';'Entropy';'EV';'EVLoss';};
col.EnvThreat=1; 
col.NTokens=2; 
col.pLoss=3;
col.Entropy=4;
col.EV=5;
col.Subject=6;
col.b=7;
col.EVLoss=8;
col.pLossAversive=9;
d_par=nan(36,3); d_par(:, [col.EnvThreat col.NTokens])=[(repmat(1:6, 1,6)'/6) sortrows(repmat(1:6, 1,6)')*2]; 
dcf=fpar_conflict(d_par, col);
dct=fpar_conflict(d_par, col);




% WHICH IVS?
% ivs={'Subject';'Entropy';'pLoss';};
ivs={'Subject';'Entropy'};
% ivs={'Subject'; 'EnvThreat';'NTokens';'pLoss';'Entropy';'EV';};

% Get betas + combine w design variables    
d_glmres{1} = [[{' '}; log.ROInames] [ivs(2:end)';  cell(length(log.ROInames),  length(ivs(2:end)))]];   d_glmres{2} = d_glmres{1} ; d_glmresp =d_glmres ;
f.subplotcols=3;  f.subplot_VerHorz=[0.05 0.2];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.15 0.01];  f.figwidth= 800;   f.figheight=1000; f.nancol=[0 0 0];
figure('color', 'w', 'Position',[200,70,f.figwidth,f.figheight]);  f.FontSize=15; k=1;
for r=1:length(log.ROInames)
    disp([log.ROInames{r} ' -------------']);
    
    % Compile betas and 
    d_bcf=[]; d_bct=[];
    for s=1:log.n_subjs
        ws=d_roimatrix_cF(:,:,1, r,s);
        ws=flipud(ws);  % Revert EnvThreat
        d_bcf=[d_bcf; [repmat(s, 36,1) ws(:)]];
        
        ws=d_roimatrix_ct(:,:,1, r,s);
        ws=flipud(ws);  % Revert EnvThreat
        d_bct=[d_bct; [repmat(s, 36,1) ws(:)]]; ws=[];
    end
    d_bcf(:, [col.Subject col.b])= d_bcf;
    for v=1:length(desvars)
        eval(['d_bcf(:, col.' desvars{v} ') = repmat(dcf(:, col.' desvars{v}  '),log.n_subjs,1);']);
    end
    d_bct(:, [col.Subject col.b])= d_bct;
    for v=1:length(desvars)
        eval(['d_bct(:, col.' desvars{v} ') = repmat(dct(:, col.' desvars{v}  '),log.n_subjs,1);']);
    end
    d_betas{r,1}{1}=d_bcf;
    d_betas{r,1}{2}=d_bct;
    
    % One big glm per ROI
    disp('One Big GLM per ROI')
    ivcols=[];  for i=1:length(ivs),  eval(['ivcols=[ivcols col.'  ivs{i} '];']),  end
    [beta, d, stats]  =glmfit(d_bcf(:, ivcols),  d_bcf(:, col.b));
    disp('cF:'), disp([beta  stats.p ]);
    [beta, d, stats]  =glmfit(d_bct(:, ivcols),  d_bct(:, col.b));
    disp('ct:'), disp([beta  stats.p ])
    
    
%     FOR residuals
%     [b1 bint rr rint stats1]=regress(d_bcf(:, col.b),  [ones(size(d_bcf,1),1)   d_bcf(:, col.Entropy)]); stats1    
%     [b2 bint rr2 rint stats2]=regress(rr,  [ones(size(d_bcf,1),1)   d_bcf(:, col.pLoss)]); stats2
    
%     FOR residuals
    disp('Residuals? ');
    d_both=[d_bcf; d_bct];
    d_both(:, col.pLossAversive)=d_both(:, col.pLoss);
    d_both(size(d_bcf,1)+1:end , col.pLossAversive) =0;
    [b1 bint rr rint stats1]=regress(d_both(:, col.b),  [ones(size(d_both,1),1)   d_both(:, col.Entropy)]); stats1    
    [b2 bint rr2 rint stats2]=regress(rr,  [ones(size(d_both,1),1)   d_both(:, col.pLoss)]); stats2
 
    
    
    
    
    % First level glms
    for s=1:log.n_subjs
        % cF 
        d_betas{r,2}{s, 1} = d_bcf( d_bcf(:, col.Subject)==s, :);
        [ws.b ws.d ws.stats]= glmfit(d_betas{r,2}{s, 1}(:, ivcols(2:end)), d_betas{r,2}{s, 1}(:,col.b));
        d_betas{r,3}{1}(s,:)=  ws.b' ; 
        
        % ct
        d_betas{r,2}{s, 2} = d_bct( d_bct(:, col.Subject)==s, :);
        [ws.b ws.d ws.stats]= glmfit(d_betas{r,2}{s, 2}(:, ivcols(2:end)), d_betas{r,2}{s, 2}(:,col.b));
        d_betas{r,3}{2}(s,:)=  ws.b' ; 
    end
    
    % Second level: one-sample ttests per iv
    [h p ci st]=ttest(d_betas{r,3}{1}); 
    st.tstat = num2cell(st.tstat); st.tstat( p>0.05) = repmat({[]}, 1, sum(p>0.05)); 
    d_glmres{1}(r+1, 2:end) =  st.tstat(2:end);
    d_glmresp{1}(r+1, 2:end) =  num2cell(p(2:end));
    [h p ci st]=ttest(d_betas{r,3}{2}); 
    st.tstat = num2cell(st.tstat); st.tstat( p>0.05) = repmat({[]}, 1, sum(p>0.05)); 
    d_glmres{2}(r+1, 2:end) =  st.tstat(2:end);
    d_glmresp{2}(r+1, 2:end) =  num2cell(p(2:end));
    
    % PLOT 
    subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
    text(0.5, 0.8,  log.ROInames{r},'FontSize', f.FontSize), axis off
    subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
    barwitherr(std(d_betas{r,3}{1})./sqrt(log.n_subjs), mean(d_betas{r,3}{1}), 'y')
    title(['[cF] ' ] ,'FontSize', f.FontSize),  set(gca, 'xticklabel', [{'Cst'}; ivs(2:end)])
    subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
    barwitherr(std(d_betas{r,3}{2})./sqrt(log.n_subjs), mean(d_betas{r,3}{2}), 'y')
    title(['[ct] '],'FontSize', f.FontSize), set(gca, 'xticklabel', [{'Cst'}; ivs(2:end)])
end
openvar d_glmres{1} , openvar d_glmres{2} 
% openvar d_glmresp{1} , openvar d_glmresp{2}         
        
% Systematically check residuals

    
    

error('Done! :)')

% Check distribution
do_checklineup=0;
if do_checklineup
    d_checkcf={}; d_checkct={};
    for e=1:6
        for n=1:6
            wc =d_bcf(  d_bcf(:, col.EnvThreat).*6==e & d_bcf(:, col.NTokens)/2==n, :);
            d_checkcf{col.EnvThreat}(7-e,n)= unique(wc(:, col.EnvThreat));
            d_checkcf{col.NTokens}(7-e,n)= unique(wc(:, col.NTokens));
            d_checkcf{col.b}(7-e,n)=mean( wc(:, col.b));
            
            
            wc =d_bct(  d_bcf(:, col.EnvThreat).*6==e & d_bct(:, col.NTokens)/2==n, :);
            
            d_checkct{col.EnvThreat}(7-e,n)= unique(wc(:, col.EnvThreat));
            d_checkct{col.NTokens}(7-e,n)= unique(wc(:, col.NTokens));
            d_checkct{col.b}(7-e,n)=mean( wc(:, col.b));
        end
    end
    figure('color','w'), k=1; 
    subplot(2,4,k),  imagesc(d_checkcf{col.EnvThreat}), axis square, colorbar, axis off, k=k+1;
    subplot(2,4,k),  imagesc(d_checkcf{col.NTokens}), axis square, colorbar, axis off, k=k+1;
    subplot(2,4,k),  imagesc(d_checkcf{col.b}), axis square, colorbar, axis off, k=k+1;
    
    k=k+1; % ct
    subplot(2,4,k),  imagesc(d_checkct{col.EnvThreat}), axis square, colorbar, axis off, k=k+1;
    subplot(2,4,k),  imagesc(d_checkct{col.NTokens}), axis square, colorbar, axis off, k=k+1;
    subplot(2,4,k),  imagesc(d_checkct{col.b}), axis square, colorbar, axis off, k=k+1;
end



%% (7) Apply contrast weights
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject

% where.modres='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs\3 Hierarchical';
where.modres='/Users/EleanorL/Dropbox/sandisk/4 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs/3 Hierarchical';
cf=load([where.modres filesep 'res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968.mat']);
ct= load([where.modres filesep 'res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293.mat']);
cfee=load([where.modres filesep 'ETwarp conflict.mat']); ctee=load([where.modres filesep 'ETwarp control.mat']); % Assumed to be the correct model
cfev=cfee.d_etwarpv; ctev=ctee.d_etwarpv;
cfe= cfee.d_etwarp; cte= ctee.d_etwarp;

for r=1:log.n_rois
    for s=1:log.n_subjs
        ws.r=squeeze(d_roimatrix_cF(:,:,1,r,s));
        
        error
        
%         cfe{s,2}.EnvThreat
        
        
    end
end


%% Tracking of all quantities separated by rejecting


