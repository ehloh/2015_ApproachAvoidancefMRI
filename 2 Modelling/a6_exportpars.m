% Read results from model fitting and set them up for fMRI
%   Every new parameter is kept in a different workspace-variable. Prefixed d_xxx variables 
%   hold parameters in 6x6 matrix format, where dt_xxx variables hold parameters 
%   in trialstats format (columns defined in scol)
%   Use this to set up/analyse/create plots for ad hoc analyses as well
% clear all; close all hidden; clc
clear all; clc; request.where_res=[];

for o=1:1 % Which results?
% request.where_res=['2 High iter' filesep];  % Typicalfits (bpmi) -----
% request.cF_res='res_fitmodels_cF (22-Jul-2014) bpmi16 2000iter recomb';  % Best 
% request.ct_res='res_fitmodels_ct (22-Jul-2014) bpmi11 bpm11 1000iter';
% 
% request.where_res=['2 High iter' filesep];  % Typicalfits (bpm) -----
% request.cF_res='res_fitmodels_cF (22-Jul-2014) bpm16 1000iter';  % Best 
% request.ct_res='res_fitmodels_ct (22-Jul-2014) bpm11 1000iter';

% request.where_res=['Det basic fits' filesep];  % Typicalfits (bpm) -----
% request.cF_res='res_fitmodels_cF (18-Jul-2014) all20iter';  % Best 
% request.ct_res='res_fitmodels_ct (18-Jul-2014) all20iter';
end
% request.where_res=['3 Hierarchical' filesep '1 Manuscript fits' filesep];  % Hierarchical ------- 
% request.cF_res='res_hierarfitmodels_cF (21-Jul-2014) top10';  % [ Currently in manuscript v2-3 ]
% request.ct_res='res_hierarfitmodels_ct (21-Jul-2014) top10';

% % Hierar --------------------
request.where_res=['3 Hierarchical' filesep ];  
request.cF_res='res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968';  % [ Current, Manuscript v2-6 ]
request.ct_res='res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293';

% % Fit --------------------
% request.cF_res='res_fitmodels_cF (19-Mar-2015) all til now'; 
% request.ct_res='res_fitmodels_ct (19-Mar-2015) all til now';

for o1=1:1 %% Set up
    Valfxn_type=[]; 
    w=pwd;  if strcmp(w(2), ':')==1; 
        where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI';  
        where.expt_fol='G:\2 [Explore]'; 
        where.data_brain=[where.expt_fol filesep '1 Brain data']; 
        where.beh_mods='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis\4 Fit computational models';
    else where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; 
        where.expt_fol='/Users/EleanorL/Desktop/2 EXPLORE fMRI data'; 
        where.data_brain=[where.expt_fol filesep '1 Brain data'];
        where.beh_mods='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis/4 Fit computational models';  
    end
    where.data_beh =[where.where filesep '1 Behavioural data'];
    path(pathdef);  addpath([where.beh_mods filesep '1 Value functions' Valfxn_type]); addpath(where.beh_mods)
    addpath(genpath([where.beh_mods filesep '1 Value functions' Valfxn_type]));
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'b']);   
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bp']);
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bm']);
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bjm']);
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bi']);   
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bpi']); 
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bmi']); 
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bpm']);
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bpmi']); 
    addpath([where.beh_mods filesep '1 Value functions' Valfxn_type filesep 'bpjm']); 
    
    % Load subjects    
    logg.specificsubjects=[];
    logg.w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); logg.datalog=logg.w.datalog;
    [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.datalog, logg.specificsubjects, [logg.datalog vertcat('include_all', num2cell(ones(size(logg.datalog,1)-1,1)))], 'include_all');
    
    % Requests & Details that don't change much
    logg.scan=load([where.where filesep '3 Scripts - Preprocessing' filesep 'i_scanningdetails.mat']); logg.scan.TRseconds=logg.scan.TRms/1000; 
    
    % Mapping ct models: if best fit includes f param, what is the true name of the model?
    [details.model_defaults  details.par_specs] =f_modelsettings([]); 
    details.fitfiles={request.cF_res request.ct_res};
    
    rc.modname=1;
    rc.subpars=2;
    rc.sp.bic=1;
    rc.sp.nll=2;
    rc.sp.p1=4;
    rc.modelbic=3;
    rc.hessians=4;
    rc.mean_p1=6;
%     disp('=============================================================='); w.c=clock; w.c1=w.c;
%     disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ');
%     disp(['Requested: ' num2str(details.n_subjs) ' subjects, ' num2str(details.n_iterations) ' iterations per model'])
%     disp(['Task type:  ' task ]); disp(' ');
%     disp('Models:'); disp(details.whichmodels); disp(' ');
%     disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
%     input('Hit enter to start                                   ');
%     disp('==============================================================')
end

%% What is the name of this fMRI model (composite of winning cF and ct models)

% Load results and raw behavioural data
%       subjdata: subject, cF only, ct only, all trials, misc (labelled)
%       Note: By the end of this script, j parameter has been applied to
%       all trialstats in subjdata
cd([where.beh_mods filesep '2 Analysis inputs' filesep request.where_res])
d_fits.cf=load([request.cF_res '.mat']);   d_fits.ct=load([request.ct_res '.mat']);
details.n_subjs=d_fits.cf.details.n_subjs; details.subjects=d_fits.cf.details.subjects; d_fits.cf.rc=rc; d_fits.ct.rc=rc;

for o=1:1   % Cols (data from raw trialstats)
    % scol is the same as the scol settings further down in the script. 
    
    scol.Choice=1;
    scol.EnvThreat=2;
    scol.NTokens=3;
    scol.pLoss=4;
    scol.Entropy=5;
    scol.VExplore=6;
    scol.EntropyNTok=7;
    scol.EV=8;
    scol.Task=9;
    %
    scol.OutcomePres=10;
    scol.OutcomeMagnitude=11;
    scol.OutcomeMean=12;
    scol.OutcomeVariance=13;
    scol.ExploredBomb=14;  
    %
    scol.TrialValid=15;
    scol.EnvThreatOriginal=16;
    scol.Trialnum=17;    
    scol.Trialtype=18;
    scol.RT1=19;
    scol.Choice2=20;
    scol.Block=21; 
    %
    scol.Ons_Fixation=22;
    scol.Ons_Gamble=23;
    scol.Ons_Resp1=24;
    scol.Ons_ExploreInfo=25;
    scol.Ons_Resp2=26;
    scol.Ons_Outcome=27;
    scol.Ons_EndTrial=28;
end
    
for o=1:1 % Load behavioural data
    subjdata=cell(d_fits.cf.details.n_subjs,2);
    
    % Columns from functional raw data  (i.e. not used further in data)
    scol.f.Trialtype=1;
    scol.f.NTokens=2;
    scol.f.EnvThreat=3;
    scol.f.Task=5;
    scol.f.Trialnum=7;
    scol.f.Choice=8;
    scol.f.Choice2=10;
    scol.f.RT1=9;
    scol.f.TrialValid=13;
    scol.f.OutcomePres=14;
    scol.f.OutcomeMagnitude=15;
    scol.f.ExploredBomb=18;
    scol.f.OnsFixation=21;
    scol.f.OnsGamble=22;
    scol.f.OnsResp1=23;
    scol.f.OnsExploreInfo=25;
    scol.f.OnsResp2=26;
    scol.f.OnsOutcome=28;
    scol.f.OnsEndTrial=29;
    
    for s=1:details.n_subjs
        ws=load([where.data_beh filesep details.subjects{s} filesep details.subjects{s} '_file_taskfMRI.mat']);
        ws.rd=ws.alldata; ws.startblock=find(ws.rd(:,scol.f.Trialnum)==1);   
        if numel(ws.startblock)~=6, ws.startblock=[1;  find(ws.rd(:,scol.f.OnsFixation) - [0; ws.rd(1:end-1,scol.f.OnsFixation)]<0)]; end
        
        for r=1:6  % Timings 
            try 
                wr.where=['G:\2b Explore fMRI archive\1 Brain data Preproc PPIold\' logg.subjects{s} '\1 Preprocessed\Func_r' num2str(r) '\'];
                ws.nscans(r)= length(spm_select('List', wr.where, ['^f' logg.datalog{s+1, find(strcmp(logg.datalog(1,:), 'Func_scanID'))} '.*.hdr'])); 
            catch;  disp(['Artificial n runs!!    ' logg.subjects{s} '  run ' num2str(r) ]);  if s==1; input('Continue?  '); end
            end
        end
        if sum(ws.nscans)==0, disp(['Artificial n runs!!    ' logg.subjects{s} '  run ' num2str(r) ]);  if s==1; input('Continue?  '); end; end 
        
        % Mark variables (include only necessary)
        ws.d=nan(size(ws.rd, 1), 19);
        ws.d(:, scol.Trialtype)=ws.rd(:, scol.f.Trialtype );
        ws.d(:, scol.EnvThreat)=ws.rd(:, scol.f.EnvThreat);
        ws.d(:, scol.NTokens )=ws.rd(:, scol.f.NTokens );
        ws.d(:, scol.Task )=ws.rd(:, scol.f.Task);
        ws.d(:, scol.Trialnum )=ws.rd(:, scol.f.Trialnum );
        ws.d(ws.startblock(1):ws.startblock(2)-1, scol.Block) =1;  ws.d(ws.startblock(2):ws.startblock(3)-1, scol.Block) =2; 
        ws.d(ws.startblock(3):ws.startblock(4)-1, scol.Block) =3;  ws.d(ws.startblock(4):ws.startblock(5)-1, scol.Block) =4; 
        ws.d(ws.startblock(5):ws.startblock(6)-1, scol.Block) =5;  ws.d(ws.startblock(6):end, scol.Block) =6; 
        ws.d(:,scol.Trialnum)=1:size(ws.d,1);
        ws.d(:, scol.Choice )=ws.rd(:, scol.f.Choice );
        ws.d(:, scol.Choice2 )=ws.rd(:, scol.f.Choice2 );
        ws.d(:, scol.RT1)=ws.rd(:, scol.f.RT1);
        ws.d(:, scol.TrialValid )=ws.rd(:, scol.f.TrialValid  );
        ws.d(:, scol.OutcomePres)=ws.rd(:, scol.f.OutcomePres);
        ws.d(:, scol.OutcomeMagnitude)=ws.rd(:, scol.f.OutcomeMagnitude);
        ws.d(:, scol.ExploredBomb)=ws.rd(:, scol.f.ExploredBomb);
        ws.d(:, scol.Ons_Fixation)=ws.rd(:, scol.f.OnsFixation);
        ws.d(:, scol.Ons_Gamble)=ws.rd(:, scol.f.OnsGamble);
        ws.d(:, scol.Ons_Resp1)=ws.rd(:, scol.f.OnsResp1);
        ws.d(:, scol.Ons_ExploreInfo)=ws.rd(:, scol.f.OnsExploreInfo);
        ws.d(:, scol.Ons_Resp2)=ws.rd(:, scol.f.OnsResp2); 
        ws.d(:, scol.Ons_Outcome)=ws.rd(:, scol.f.OnsOutcome);
        ws.d(:, scol.Ons_EndTrial)=ws.rd(:, scol.f.OnsEndTrial);
        ws.d(:,  [scol.Ons_Resp1 scol.Ons_EndTrial]) = ws.d(:,  [scol.Ons_Resp1 scol.Ons_EndTrial])./1000;
        ws.timings=[scol.Ons_Fixation scol.Ons_Gamble scol.Ons_Resp1 scol.Ons_ExploreInfo scol.Ons_Resp2 scol.Ons_Outcome scol.Ons_EndTrial]; 
        
        % Correct timings 
        ws.nscans=[0 ws.nscans];  ws.nscans= cumsum(ws.nscans); 
        do_4fMRI=0; 
        if do_4fMRI
        for r=1:6
            ws.d(find(ws.d(:, scol.Block)==r),  ws.timings) = ws.d(find(ws.d(:, scol.Block)==r),  ws.timings)-ws.settings.times{r}.start/1000; % Correct within block
            ws.d(find(ws.d(:, scol.Block)==r),  ws.timings) =  ws.d(find(ws.d(:, scol.Block)==r),  ws.timings)  + ws.nscans(r)* logg.scan.TRms .* logg.scan.nSlicesPerVol./1000;
        end
        else if s==1, input('Not sorting timings for fMRI. Continue?  ');  end 
        end 
        ws.cf=ws.d(ws.d(:,scol.Task)==1, :);
        ws.ct=ws.d(ws.d(:,scol.Task)==2, :);
        [ ws.cf] = fpar_conflict( ws.cf, scol);
        [ ws.ct] = fpar_control( ws.ct, scol);
        ws.d=sortrows([ws.cf; ws.ct], scol.Trialnum);
        
        
        % Check explored bomb?
        if sum(ws.rd(:,scol.f.ExploredBomb)==1  & ws.rd(:, scol.f.Choice)~=3)~=0; disp(['Check column for ExploredBomb for subject ' num2str(s) '(' details.subjects{s} ')']); end
        
        %
        subjdata{s,1}=details.subjects{s};
        subjdata{s,5}.errortrials=ws.d(ws.d(:, scol.TrialValid)~=1,:); 
        %
        ws.d=ws.d(ws.d(:, scol.TrialValid)==1,:); % Valid trials only
        subjdata{s,2}=ws.d(ws.d(:,scol.Task)==1,:);
        subjdata{s,3}=ws.d(ws.d(:,scol.Task)==2,:);
        subjdata{s,4}=ws.d;
        
        % for re-insertion back into trialstats (to match onsets)
        ws=[];
    end
end

% % % Artificially request certain models/alter params!!!
% disp('Altering requested model/params!')
% d_fits.cf.r_res=d_fits.cf.r_res(strcmp(d_fits.cf.r_res(:,1), 'b01'),:);   % Request certain models!!
% d_fits.ct.r_res=d_fits.ct.r_res(strcmp(d_fits.ct.r_res(:,1), 'bpmi11_euw'),:);
% p=1;
% d_fits.cf.r_res{1,2}(:,3+p)=20;
% d_fits.ct.r_res{1,2}(:,3+p)=20;

% Checks and setup
d_fits.cf.r_res=sortrows(d_fits.cf.r_res, d_fits.cf.rc.modelbic); d_fits.ct.r_res=sortrows(d_fits.ct.r_res, d_fits.ct.rc.modelbic);
if sum(strcmp(d_fits.cf.details.subjects, d_fits.ct.details.subjects))~=d_fits.cf.details.n_subjs || sum(strcmp(d_fits.cf.details.subjects, subjdata(:,1)))~=20;  error('Check subjects order for fits for cF, ct and subjdata'); end
details.subjects=d_fits.cf.details.subjects; details.n_subjs=length(details.subjects);
d_par.cf=d_fits.cf.r_res{1,2}(:, d_fits.cf.rc.sp.p1:end);   d_par.ct=d_fits.ct.r_res{1,2}(:, d_fits.ct.rc.sp.p1:end); % Read params
details.cf_mod=d_fits.cf.r_res{1,1}; details.ct_mod=d_fits.ct.r_res{1,1}; 
details.modname=[details.cf_mod(1:strfind(details.cf_mod, '_')-1) details.ct_mod(1:strfind(details.ct_mod, '_')-1)];
if strcmp(details.cf_mod, 'b01')==1; details.modname= [details.cf_mod details.modname]; end
if strcmp(details.ct_mod, 'b01')==1; details.modname= [details.modname details.ct_mod]; end
details.cf_moddetails=d_fits.cf.details.models(strcmp(d_fits.cf.details.models, details.cf_mod),:); details.cf_moddetails=details.cf_moddetails(1,:);
details.ct_moddetails=d_fits.ct.details.models(strcmp(d_fits.ct.details.models, details.ct_mod),:); details.ct_moddetails=details.ct_moddetails(1,:);
details.whichmodel=[details.cf_moddetails(1) details.ct_moddetails(1)];
disp(['Name of model: ' details.modname]);
d_par.partable=[ [{'Subjects'}; d_fits.cf.details.subjects]   [cellfun(@(x)['cF_' x], (details.model_defaults{find(strcmp(details.model_defaults(:,1),details.cf_mod)),3})', 'UniformOutput', 0);  num2cell(d_par.cf)]     [cellfun(@(x)['ct_' x], (details.model_defaults{find(strcmp(details.model_defaults(:,1),details.ct_mod)),3})', 'UniformOutput', 0);  num2cell(d_par.ct)]];    % subjects - cF -  ct

%% Set up variables for model


% Columns (for new variables in subjdata space)
scol.vChosen=30;
scol.vBestUnchosen= scol.vChosen+1;
scol.vWorstUnchosen= scol.vChosen+2;
scol.vBestUnchosen_pos= scol.vChosen+3;
scol.vBestUnchosen_neg= scol.vChosen+4;
scol.BestUnchosen_Is= scol.vChosen+5;       % Whats the best unchosen? (0=draw)
scol.vBestUnchosen_UnchoChoice(1)= scol.vChosen+6;  % Best unchosen is an Accept choice
scol.vBestUnchosen_UnchoChoice(2)= scol.vChosen+7;
scol.vBestUnchosen_UnchoChoice(3)= scol.vChosen+8;
% scol.EVpos= scol.vChosen+9;
% scol.EVneg= scol.vChosen+10;
scol.vGamble= scol.vChosen+11;
scol.vGamblepos= scol.vChosen+12;
scol.vGambleneg= scol.vChosen+13;
scol.Modalchoice=scol.vChosen+14;  % If draw
scol.vModalchoice=scol.vChosen+15;
scol.vSee=scol.vChosen+16;
scol.vNoSee=scol.vChosen+17;
scol.vBest=scol.vChosen+18;
scol.vWorst=scol.vChosen+19;
scol.vBesttoWorst=scol.vChosen+20;
scol.MargvChosen=scol.vChosen+21;
scol.EVgain=scol.vChosen+22;
scol.EVgain=scol.vChosen+23;
scol.EVConflict=scol.vChosen+24;
scol.PavConflict=scol.vChosen+25;
scol.ppActBomb=scol.vChosen+26;  % posterior p(ActBomb|~See)
% scol.pSee=scol.vChosen+27;  % p(See Bomb)  * THIS S 
        
        
% Compile values for each trialtype/cell (in matrix format)
%       d_xxx{s,task+1}{choice} =V(Choice),  EnvThreat x NTokens (6x6, EnvThreat is NOT flipped inverse, i.e. data is NOT like imagesc plots)
disp('Compiling values in matrix format --------------')
d_trialvals=[details.subjects cell(details.n_subjs, 2)];  % V(Accept/Reject/Explore) in trialstats form (36 x 3)
d_vchoice=[details.subjects cell(details.n_subjs, 2)];               % [1a] V(Accept/Reject/Explore), in 6x6 
        d_vchoice2=[details.subjects cell(details.n_subjs, 2)];               % [1b] V(Accept/Reject/Explore), in 6x6  (task x envthreat x ntokens x choice) i.e. multi-indexing. just for latter processing. 
d_predchoice=[details.subjects cell(details.n_subjs, 2)];        %  [3] p(Chosen==Best/2nd Best), in 6x6  (task x choice x envthreat x ntokens)
d_vbestchoice=[details.subjects cell(details.n_subjs, 2)];        %  [2] V(1stChoice/2ndChoice/Worst), in 6x6  (task x choice x envthreat x ntokens)
d_pbestchoice=[details.subjects cell(details.n_subjs, 2)];        %  [3] p(Chosen==Best/2nd Best), in 6x6  (task x choice x envthreat x ntokens)
d_modalchoice_andval= [details.subjects repmat({repmat({nan(6,6)}, 1,2)}, details.n_subjs, 2)];  % [4] Modal choice and values thereof
% d_vmodalchoice=[details.subjects cell(details.n_subjs, 2)];
d_vseenosee=[details.subjects repmat({repmat({nan(6,6)}, 1,2)}, details.n_subjs, 2)];
d_evgainloss=[details.subjects repmat({repmat({nan(6,6)}, 1,3)}, details.n_subjs, 2)];  % EV Gain, EV Loss, EV Gain - EV Loss
% d_pavconflict=[details.subjects repmat({repmat({nan(6,6)}, 1,1)}, details.n_subjs, 2)];  % Pav Conflict, pLoss*Ntok
d_prob=[details.subjects cell(details.n_subjs, 2)];        % Prior p, posterior p (see & no see, weighted by their objective probability), Info Gain (KL diverg), Info gain (Entropy change)
for t=1:2
    wt.design=nan(6*6, 10);  
    wt.design(:,  [d_fits.cf.details.col.EnvThreat d_fits.cf.details.col.NTokens])=[sortrows(repmat((1:6)',6,1))/6 2*repmat((1:6)',6,1)]; 
    wt.design(:,d_fits.cf.details.col.Task)=t;
    switch t
        case 1; [ wt.design] = fpar_conflict( wt.design, d_fits.cf.details.col);  wt.mod=details.cf_mod ; wt.par=d_par.cf; wt.moddetails=details.cf_moddetails; wt.moddetails=wt.moddetails(1,:);
        case 2; [ wt.design] = fpar_control( wt.design, d_fits.cf.details.col);  wt.mod=details.ct_mod ; wt.par=d_par.ct; wt.moddetails=details.ct_moddetails; wt.moddetails=wt.moddetails(1,:);
    end
    if isempty(strfind(details.cf_mod,'j'))==0 | isempty(strfind(details.ct_mod,'j'))==0;    disp('README: Model calls for EnvThreat distortion. All values and trialstats are warped in the frame script'); end
    
    for s=1:details.n_subjs
        ws.behchoice=subjdata{s,t+1}; % ws.behchoice 
        ws.behchoice(:, scol.EnvThreatOriginal)= ws.behchoice(:, scol.EnvThreat);
        ws.par=wt.par(s,:); % Models params
        [ ws.transpar] = f_transpar(wt.moddetails{3}, ws.par, 'from');  % Apply inverse transformation to params
                
        % Get values
        eval(['[ ws.vchoice] = ' wt.mod '(ws.transpar, {[] wt.design d_fits.cf.details.fixedpar, d_fits.cf.details.col});'])
        ws.vchoice= [squeeze(ws.vchoice)  nan(6*6,1) wt.design(:, [d_fits.cf.details.col.EnvThreat d_fits.cf.details.col.NTokens])  ];
        ws.vchoice(:,5)=ws.vchoice(:,5).*6;  % col 5= EnvThreat level
        ws.vchoice(:,6)=ws.vchoice(:,6)./2;  % col 6= NTokens level
        d_trialvals{s, t+1}=ws.vchoice;
        
        % [1] V(Accept/Reject/Explore)
        d_vchoice{s,t+1}=repmat({nan(6,6)}, 1,3);
        for i=1:size(ws.vchoice,1)
            for c=1:3
                d_vchoice{s,t+1}{c}(ws.vchoice(i, 5), ws.vchoice(i,6))=ws.vchoice(i, c);
                d_vchoice2{s,t+1}(ws.vchoice(i, 5), ws.vchoice(i,6), c)=ws.vchoice(i, c);
            end
        end
        
        % [2] V(Best/2nd Best/Worst) (theoretically best)
        d_vbestchoice{s,t+1}=repmat({nan(6,6)}, 1,3);
        for i=1:size(ws.vchoice,1)
            d_vbestchoice{s,t+1}{1}(ws.vchoice(i, 5), ws.vchoice(i, 6))= ws.vchoice(i,  find(ws.vchoice(i,1:3)==max(ws.vchoice(i, 1:3)), 1)); 
            d_vbestchoice{s,t+1}{2}(ws.vchoice(i, 5), ws.vchoice(i, 6))= ws.vchoice(i, find(ws.vchoice(i,1:3)==median(ws.vchoice(i, 1:3)), 1)); 
            d_vbestchoice{s,t+1}{3}(ws.vchoice(i, 5), ws.vchoice(i, 6))= ws.vchoice(i, find(ws.vchoice(i,1:3)==min(ws.vchoice(i, 1:3)), 1)); 
        end
        
        % [3] Predicted choice
        d_predchoice{s,t+1}=repmat({nan(6,6)}, 1,3);
        cc=d_fits.cf.details.col;  % wt.design goes with cc/d_fits.cf.details.col, NOT scol
        for c=1:3
            wt.s=wt.design;
            wt.s(:, cc.Choice)=c; wt.s(:, cc.Task)=t;            
            if strcmp(details.whichmodel{t}(2), 'p')==1         
                [nll wt.predchoice{c}]=f_nllsoftmax_lapse(ws.transpar, {details.whichmodel{t}    wt.s      d_fits.cf.details.fixedpar    cc});
            else [nll wt.predchoice{c}]=f_nllsoftmax(ws.transpar,        {details.whichmodel{t}    wt.s      d_fits.cf.details.fixedpar    cc});
            end
            d_predchoice{s,t+1}{c}=reshape(wt.predchoice{c}(:),6,6)'; 
        end
        
        % [4] p(Best choice) - % of real choices tt conform to the theoretical optimal
        %       Note: draws in value are dealt with somewhat arbitrarily. See code for details.
        d_pbestchoice{s,t+1}=repmat({nan(6,6)}, 1,3);
        for e=1:6
            for n=1:6
                wc.choice=ws.behchoice(ws.behchoice(:, scol.EnvThreatOriginal)*6==e & ws.behchoice(:, scol.NTokens)/2==n, scol.Choice);
                wc.vchoices=ws.vchoice(find(ws.vchoice(:,5)==e & ws.vchoice(:,6)==n), 1:3);
                if length(unique(wc.vchoices))==3 % no draws in value
                    d_pbestchoice{s,t+1}{1}(e,n)=mean(wc.choice==find(wc.vchoices==max(wc.vchoices)));
                    d_pbestchoice{s,t+1}{2}(e,n)=mean(wc.choice==find(wc.vchoices==median(wc.vchoices)));
                    d_pbestchoice{s,t+1}{3}(e,n)=mean(wc.choice==find(wc.vchoices==min(wc.vchoices)));
                elseif sum(wc.vchoices==median(wc.vchoices))==1 % only best is unique
                    % in the event of a draw (in value), both choices count as an EV hit  and runner up (i.e. divide the hits evenly between 1st and 2nd)
                    wc.evhit=find(ws.vchoice(ws.vchoice(:,5)==e & ws.vchoice(:,6)==n, 1:3)==max(ws.vchoice(ws.vchoice(:,5)==e & ws.vchoice(:,6)==n, 1:3)));
                    d_pbestchoice{s,t+1}{1}(e,n)=mean(wc.choice==wc.evhit(1) | wc.choice==wc.evhit(2)) /2;
                    d_pbestchoice{s,t+1}{2}(e,n)=mean(wc.choice==wc.evhit(1) | wc.choice==wc.evhit(2)) /2;
                    d_pbestchoice{s,t+1}{3}(e,n)=mean(wc.choice==find(wc.vchoices==min(wc.vchoices)));
                elseif sum(wc.vchoices==max(wc.vchoices))==1
                    % in the event of a draw (in value) for runner up, both choices count as an EV runner up and loser (i.e. divide the hits evenly between 2nd and 3rd)
                    d_pbestchoice{s,t+1}{1}(e,n)=mean(wc.choice==find(wc.vchoices==max(wc.vchoices)));
                    wc.ev2nd=find(wc.vchoices==median(wc.vchoices));
                    d_pbestchoice{s,t+1}{2}(e,n)=mean(wc.choice==wc.ev2nd(1) | wc.choice==wc.ev2nd(2)) /2;
                    d_pbestchoice{s,t+1}{3}(e,n)=mean(wc.choice==wc.ev2nd(1) | wc.choice==wc.ev2nd(2)) /2;
                end 
                wc=[];
            end
        end
        
        % [5] Modal choice + value of modal choice
        for e=1:6
            for n=1:6
                wc.choice=ws.behchoice(ws.behchoice(:, scol.EnvThreatOriginal)*6==e & ws.behchoice(:, scol.NTokens)/2==n, scol.Choice);
                if length(mode(wc.choice))==1
                    d_modalchoice_andval{s,t+1}{1}(e,n)=mode(wc.choice);
                    d_modalchoice_andval{s,t+1}{2}(e,n)= d_vchoice{s,t+1}{mode(wc.choice)}(e,n);
%                     d_vmodalchoice{s,t+1}(e,n)= d_vchoice{s,t+1}{mode(wc.choice)}(e,n);
                else error('goddamnit >1 modal now what'); 
                end
            end
        end
        
        % [6] vGamble before and after exploring
        switch t
            case 1
                for k=1:length( details.cf_moddetails{3});  eval(['ws.subvar.' details.cf_moddetails{3}{k} '='  num2str(ws.par(k)) ';']);  end
                if isempty(strfind(details.cf_mod,'w'))==0;  eval(['ws.subvar.' details.cf_mod(strfind(details.cf_mod, 'w')-1:end) '='  num2str(ws.par(strcmp(details.cf_moddetails{3}, 'w')) ) ';']);     end
                [ ws.behchoice] = fcf_ExploreVal(ws.behchoice, scol, ws.subvar);
            case 2
                for k=1:length( details.ct_moddetails{3});   eval(['ws.subvar.' details.ct_moddetails{3}{k} '='  num2str(ws.par(k)) ';']);   end
                if isempty(strfind(details.ct_mod,'w'))==0;   eval(['ws.subvar.' details.ct_mod(strfind(details.ct_mod, 'w')-1:end) '='  num2str(ws.par(strcmp(details.ct_moddetails{3}, 'w')) ) ';']);     end
                [ ws.behchoice] = fct_ExploreVal(ws.behchoice, scol, ws.subvar);
        end
        d_vseenosee{s,t+1}{1}=nan(6,6); d_vseenosee{s,t+1}{2}=nan(6,6);
        for e=1:6
            for n=1:6
                wc.choice=ws.behchoice(ws.behchoice(:, scol.EnvThreatOriginal)*6==e & ws.behchoice(:, scol.NTokens)/2==n, :);
                d_vseenosee{s,t+1}{1}(e,n)=wc.choice(1, scol.vSee);
                d_vseenosee{s,t+1}{2}(e,n)=wc.choice(1, scol.vNoSee);
            end
        end
        
        % EV Gain and Loss
        switch t
            case 1
                if isempty(strfind(details.cf_mod,'j'))==0; 
                    ws.ET_effective=  power(wt.design(:, d_fits.cf.details.col.EnvThreat),  ws.subvar.j);
                else ws.ET_effective=  wt.design(:, d_fits.cf.details.col.EnvThreat);
                end
                ws.ntok=wt.design(:, d_fits.cf.details.col.NTokens);
                ws.mv.FixedLoss=ws.subvar.f;
                [ ws.ov] = fcf_changeEnvThreat(ws.ET_effective, ws.ntok, ws.mv);
                ws.out= [wt.design(:, d_fits.cf.details.col.EnvThreat) wt.design(:, d_fits.cf.details.col.NTokens) ws.ov.EVGain ws.ov.EVLoss ws.ov.EVConflict]; % et, ntok, evgain, evloss, evconflict
            case 2
                 if isempty(strfind(details.ct_mod,'j'))==0; 
                    ws.ET_effective=  power(wt.design(:, d_fits.ct.details.col.EnvThreat),  ws.subvar.j);
                else ws.ET_effective=  wt.design(:, d_fits.ct.details.col.EnvThreat);
                end
                ws.ntok=wt.design(:, d_fits.ct.details.col.NTokens);
                [ ws.ov] = fct_changeEnvThreat(ws.ET_effective, ws.ntok, []);
                ws.out= [wt.design(:, d_fits.ct.details.col.EnvThreat) wt.design(:, d_fits.ct.details.col.NTokens) ws.ov.EVGain ws.ov.EVLoss ws.ov.EVConflict]; 
        end
        d_evgainloss{s,t+1}{1} = (reshape(ws.out(:,3),6,6))';   % evgain
        d_evgainloss{s,t+1}{2} = (reshape(ws.out(:,4),6,6))';   % evloss
        d_evgainloss{s,t+1}{3} = (reshape(ws.out(:,5),6,6))';   % evconflict 

        % Sample Info 
        % Dkl(P||Q) = P.*ln(P/Q): KL divergence between prob distributions P & Q = ..
        % Q should be the start! See Wiki
        ws.des(:, [scol.EnvThreat scol.NTokens scol.pLoss scol.EnvThreatOriginal scol.Entropy])= unique(ws.behchoice(:, [scol.EnvThreat scol.NTokens scol.pLoss scol.EnvThreatOriginal scol.Entropy]), 'rows'); % 1: ET, 2:N, 3: pLoss, 4: ETorg
        ws.ppActBomb= (ws.des(:,scol.EnvThreat).*ws.des(:,scol.NTokens ))./ (24 - ws.des(:,scol.EnvThreat).*ws.des(:,scol.NTokens )); % 5: posterior pActBomb|~See
        ws.pSee= ws.des(:,scol.pLoss).*0.5; 
        % Separately compute InfoGain from each Explore-Outcome, sum weighted by p(Explore-Outcome occurs)
%         ws.InfoG_See=  ws.des(:, [scol.pLoss ]).*log( ws.des(:, [scol.pLoss ])./1);  % p(AB|See)=1;
%         ws.InfoG_NoSee=  ws.des(:, [scol.pLoss ]).*log( ws.des(:, [scol.pLoss ])./ws.ppActBomb);  % p(AB|See)=1;
        ws.InfoG_See=  1.*log( 1./ws.des(:, [scol.pLoss ]));  % p(AB|See)=1;
        ws.InfoG_NoSee=  ws.ppActBomb.*log( ws.ppActBomb./ws.des(:,scol.pLoss));  % p(AB|See)=1;
        ws.InfoG_weightsum=  ws.InfoG_See.* ws.pSee + ws.InfoG_NoSee.*(1-ws.pSee); 
        d_prob{s,t+1}{1} = reshape(   ws.des(:, [scol.pLoss ])    ,6,6)'        ; 
        d_prob{s,t+1}{2} =  reshape(   ws.pSee.*1 +   (1- ws.pSee) .*ws.ppActBomb ,6,6)'    ; 
        d_prob{s,t+1}{3} = reshape(   ws.InfoG_weightsum     ,6,6)'        ; 
        ws.posEntropy_NoSee= - ws.ppActBomb.*log(ws.ppActBomb)  - (1-ws.ppActBomb) .*log(1-ws.ppActBomb); 
        ws.posEntropy_NoSee(ws.ppActBomb==1)= - (1-0.00001).*log(1-0.00001) - (1-(1-0.00001)).*log(1-(1-0.00001)); 
        % d_prob{s,t+1}{4} = reshape(      ws.posEntropy_NoSee - ws.des(:, scol.Entropy)    ,6,6)'        ;   % posterior entropy not weighted!! 
        ws.Expected_ChangeEntropy= (ws.posEntropy_NoSee - ws.des(:, [scol.Entropy])).*(1-ws.pSee) + (0 - ws.des(:, [scol.Entropy])).*ws.pSee; 
        d_prob{s,t+1}{4} = reshape(     ws.Expected_ChangeEntropy    ,6,6)'        ;
       
        
        % NOTE: EnvThreat at this point is already warped by j term (if
        % requested). Use EnvThreatOriginal to index. BE CAREFUL if you are
        % directly using EnvThreat at any point!
        %    For this reason, best not to calculate PavConflict here. Complicated. 
        
        %
        subjdata{s,t+1}=ws.behchoice;
        ws=[];
    end
    wt=[];
end

% Compile quantities that are in trialstats format (scol) into plotting format and vice versa
%       dt_xxx{s,task+1}{choice} =V(Choice), Trialstats format (original trialstats order as experienced by subject is preserved)
d_trialstats=[details.subjects cell(details.n_subjs, 3)];        % Hold all old and new variables in trialstats form
    d_vchosenand=[details.subjects repmat({repmat({nan(6,6)}, 1,3)}, details.n_subjs, 2)];      %     vChosen, vBestUnchosen, vCho-vBU, but in matrix format (for plotting)
    d_vbestunchosen_posneg=[details.subjects repmat({repmat({nan(6,6)}, 1,2)}, details.n_subjs, 2)];      %  vBestUnchosen (positive and negative), in matrix format (for plotting)
    d_vBestUnchosenXUnchosChoice=[details.subjects repmat({repmat({nan(6,6)}, 1,3)}, details.n_subjs, 2)];            %  vBestUnchosen as a function of what the unchosen choice is (Acc/Rej/Exp), in matrix format (for plotting)
    d_vgamble=[details.subjects repmat({repmat({nan(6,6)}, 1,1)}, details.n_subjs, 2)];
    d_vgamble_posneg=[details.subjects repmat({repmat({nan(6,6)}, 1,2)}, details.n_subjs, 2)];
    d_bestchoice=[details.subjects  repmat({ {nan(6,6)} }, details.n_subjs, 2) ];  
    d_choiceinoptim=[details.subjects  repmat({ repmat({nan(6,6)},1,4) }, details.n_subjs, 2) ];  % vCho<vBU: Accept, Reject, Explore, All trials 
    d_outcome=[details.subjects  repmat({ repmat({nan(6,6)},1,4) }, details.n_subjs, 2) ]; % Outcome overall, % Pos, % 0, % Neg
    % Do NOT add things here that were already compiled in the previous block !!
    
disp('Compiling values in trialstats format --------------')
for s=1:details.n_subjs
    disp(['Subject ' num2str(s) '  - ' details.subjects{s}])
    ws.d=subjdata{s,4};    

    % Setup for j term
    ws.d(:,scol.EnvThreatOriginal)= ws.d(:,scol.EnvThreat); ws.dcf= ws.d(ws.d(:,scol.Task)==1, :); ws.dct= ws.d(ws.d(:,scol.Task)==2, :);
    if isempty(strfind(details.cf_mod, 'j'))==0; % Alter j for cF 
        if isempty(strfind(details.cf_mod, 'f'))==0; ws.mv.FixedLoss= d_par.cf(s, find(strcmp(details.cf_moddetails{3}, 'f'))); else ws.mv=[]; end
        ws.dcf(:, scol.EnvThreat)= power(ws.dcf(:, scol.EnvThreat),  d_par.cf(s,  find(strcmp(details.cf_moddetails{3}, 'j'))));
        [ ws.nv] = fcf_changeEnvThreat(ws.dcf(:, scol.EnvThreat), ws.dcf(:, scol.NTokens), ws.mv);
        ws.changepars={'pLoss'; 'Entropy'; 'EntropyNTok'; 'VExplore'; 'EV'};
        for c=1:length(ws.changepars); eval(['ws.dcf(:, scol.'  ws.changepars{c} ')=ws.nv.' ws.changepars{c} ';']); end
    end
    if isempty(strfind(details.ct_mod, 'j'))==0; % Alter j for cF
        ws.dct(:, scol.EnvThreat)= power(ws.dct(:, scol.EnvThreat),  d_par.ct(s,  find(strcmp(details.ct_moddetails{3}, 'j'))));
        [ ws.nv] = fct_changeEnvThreat(ws.dct(:, scol.EnvThreat), ws.dct(:, scol.NTokens), []);
        ws.changepars={'pLoss'; 'Entropy'; 'EntropyNTok'; 'VExplore'; 'EV'};
        for c=1:length(ws.changepars); eval(['ws.dct(:, scol.'  ws.changepars{c} ')=ws.nv.' ws.changepars{c} ';']); end
    end
    ws.d=sortrows([ws.dcf; ws.dct], scol.Trialnum);
    %
    ws.voptions=d_vchoice(s,2:3);  % x{task}{choice}(envthreat, ntokens)
    ws.voptions2=d_vchoice2(s,2:3);  %  x{task}(env threat, ntokens, choice)
    ws.voptionsother{1}{1}=ws.voptions2{1}(:,:,2:3);  % x{task}{choice}(envthreat, ntokens) - this just fetches value of alternative options, doens't specify which options those remaining values map to
    ws.voptionsother{1}{2}=ws.voptions2{1}(:,:,[1 3]);
    ws.voptionsother{1}{3}=ws.voptions2{1}(:,:,1:2);
    ws.voptionsother{2}{1}=ws.voptions2{2}(:,:,2:3);  % x{task}{choice}(envthreat, ntokens) - this just fetches value of alternative options, doens't specify which options those remaining values map to
    ws.voptionsother{2}{2}=ws.voptions2{2}(:,:,[1 3]);
    ws.voptionsother{2}{3}=ws.voptions2{2}(:,:,1:2);
    
    % [1a] d_trialstats: Value of Chosen & Next best
    for tn=1: size(ws.d,1);
        ws.d(tn, scol.vChosen)=ws.voptions{ws.d(tn, scol.Task)}{ws.d(tn, scol.Choice)}(ws.d(tn, scol.EnvThreatOriginal)*6, ws.d(tn, scol.NTokens)/2);
        ws.d(tn, scol.vBestUnchosen)=max(ws.voptionsother{ws.d(tn, scol.Task)}{ws.d(tn, scol.Choice)}(ws.d(tn, scol.EnvThreatOriginal)*6, ws.d(tn, scol.NTokens)/2, :));  % Fuck yeah 1 line of code. Also, draws are unproblematic since we don't care which choice it corresponds to - here we just care about the value
        ws.d(tn, scol.vBestUnchosen_pos)=(ws.d(tn, scol.vBestUnchosen)>=0) .* ws.d(tn, scol.vBestUnchosen);
        ws.d(tn, scol.vBestUnchosen_neg)=(ws.d(tn, scol.vBestUnchosen)<=0) .* ws.d(tn, scol.vBestUnchosen);
        ws.d(tn, scol.vWorstUnchosen)=min(ws.voptionsother{ws.d(tn, scol.Task)}{ws.d(tn, scol.Choice)}(ws.d(tn, scol.EnvThreatOriginal)*6, ws.d(tn, scol.NTokens)/2, :));         
        ws.d(tn, scol.vBest) = max(ws.voptions2{ws.d(tn, scol.Task)}(ws.d(tn, scol.EnvThreatOriginal)*6, ws.d(tn, scol.NTokens)/2,:));
        ws.d(tn, scol.vWorst) = min(ws.voptions2{ws.d(tn, scol.Task)}(ws.d(tn, scol.EnvThreatOriginal)*6, ws.d(tn, scol.NTokens)/2,:));
        ws.d(tn, scol.vBesttoWorst)=  ws.d(tn, scol.vBest) -ws.d(tn, scol.vWorst) ; 
        ws.d(tn, scol.MargvChosen)=ws.d(tn, scol.vChosen)-ws.d(tn, scol.vBestUnchosen); 
        
        % Best unchosen as a fxn of choice
        if sum(ws.d(tn, scol.vBestUnchosen)==squeeze(d_vchoice2{s, ws.d(tn, scol.Task)+1}( ws.d(tn, scol.EnvThreatOriginal)*6,  ws.d(tn, scol.NTokens)/2, :)))==1
            ws.d(tn, scol.BestUnchosen_Is)=find(ws.d(tn, scol.vBestUnchosen)==squeeze(d_vchoice2{s, ws.d(tn, scol.Task)+1}( ws.d(tn, scol.EnvThreatOriginal)*6,  ws.d(tn, scol.NTokens)/2, :)));
            ws.d(tn,  scol.vBestUnchosen_UnchoChoice( ws.d(tn, scol.BestUnchosen_Is)) )=ws.d(tn, scol.vBestUnchosen);
        else ws.d(tn, scol.BestUnchosen_Is)=0; % Draw = count the value in both regressors
            ws.tn_bu_choices=find(ws.d(tn, scol.vBestUnchosen)==squeeze(d_vchoice2{s, ws.d(tn, scol.Task)+1}( ws.d(tn, scol.EnvThreatOriginal)*6,  ws.d(tn, scol.NTokens)/2, :)));
            ws.d(tn,  scol.vBestUnchosen_UnchoChoice( ws.tn_bu_choices(1) ))=ws.d(tn, scol.vBestUnchosen);
            ws.d(tn,  scol.vBestUnchosen_UnchoChoice( ws.tn_bu_choices(2) ))=ws.d(tn, scol.vBestUnchosen);
        end
        
    end
    
    % [1b] Value of gambles
    if sum(strcmp(d_par.partable(1,:), 'cF_f'))==1
        ws.SubjectiveLoss=d_par.partable{find(cellfun(@(x)~isempty(x), strfind(d_par.partable(:,1), details.subjects{s}))), find(strcmp(d_par.partable(1,:), 'cF_f'))};
    else ws.SubjectiveLoss=-12;    
    end
    ws.d(ws.d(:,scol.Task)==1,scol.vGamble)= ws.d(ws.d(:,scol.Task)==1,scol.pLoss).*ws.SubjectiveLoss + (1-ws.d(ws.d(:,scol.Task)==1,scol.pLoss)).*ws.d(ws.d(:,scol.Task)==1,scol.NTokens);
    ws.d(ws.d(:,scol.Task)==2,scol.vGamble)= (1-ws.d(ws.d(:,scol.Task)==2, scol.Entropy)).*ws.d(ws.d(:,scol.Task)==2,scol.NTokens);
    ws.d(:,scol.vGamblepos)=ws.d(:,scol.vGamble) .* double(ws.d(:,scol.vGamble)>0);
    ws.d(:,scol.vGambleneg)=ws.d(:,scol.vGamble) .* double(ws.d(:,scol.vGamble)<0);
    
    % [2] Assemble misc plots (for plotting or otherwise)
    for t=1:2
        for e=1:6
            for n=1:6
                wc=ws.d(ws.d(:, scol.Task)==t & ws.d(:, scol.EnvThreatOriginal)*6==e & ws.d(:, scol.NTokens)/2 ==n, :);
                if length(mode(wc(:,scol.Choice)))==1
                    d_vchosenand{s,t+1}{1}(e,n)=ws.voptions{t}{mode(wc(:,scol.Choice))}(e,n);
                    d_vchosenand{s,t+1}{2}(e,n)=max(ws.voptionsother{t}{mode(wc(:,scol.Choice))}(e,n, :));
                    % d_vchosenand{s,t+1}{3}(e,n)=d_vchosenand{s,t+1}{1}(e,n) - d_vchosenand{s,t+1}{2}(e,n); % This is wrong because this should be on a trial by trial basis! 
                    d_vchosenand{s,t+1}{3}(e,n)= mean(wc(:,scol.MargvChosen)); 
                    
                    % Best unchosen split by valence
                    d_vbestunchosen_posneg{s,t+1}{1}(e,n)=mean(wc(:, scol.vBestUnchosen_pos));
                    d_vbestunchosen_posneg{s,t+1}{2}(e,n)=mean(wc(:, scol.vBestUnchosen_neg));
                    
                    % Best unchosen split by what the (best-unchosen) choice is (note: double-counting for draws is already implemented in trialstats)
                    for ch=1:3   
                        d_vBestUnchosenXUnchosChoice{s,t+1}{ch}(e,n)=mean(wc(:, scol.vBestUnchosen_UnchoChoice(ch)));
                    end
                else disp ('ugh. draws in modal choice not worked out yet. nan it out, and change plots to imagescnan')
                end
                
                % Best choice
                if sum(squeeze(ws.voptions2{t}(e,n,:)) == max(ws.voptions2{t}(e,n,:)))==1
                    d_bestchoice{s,t+1}{1}(e,n)= find(squeeze(ws.voptions2{t}(e,n,:)) == max(ws.voptions2{t}(e,n,:)));
                else
                    if sum(find(squeeze(ws.voptions2{t}(e,n,:)) == max(ws.voptions2{t}(e,n,:))))==3; % A + R
                        d_bestchoice{s,t+1}{1}(e,n)= 1.5;
                    elseif sum(find(squeeze(ws.voptions2{t}(e,n,:)) == max(ws.voptions2{t}(e,n,:))))==4; % A + E
                        d_bestchoice{s,t+1}{1}(e,n)= 1.9;
                    elseif sum(find(squeeze(ws.voptions2{t}(e,n,:)) == max(ws.voptions2{t}(e,n,:))))==5; % R + E
                        d_bestchoice{s,t+1}{1}(e,n)= 2.5;
                    end
                end

                % EVs
                d_vgamble{s,t+1}{1}(e,n)=wc(1,scol.vGamble);
                d_vgamble_posneg{s,t+1}{1}(e,n)=wc(1,scol.vGamblepos);
                d_vgamble_posneg{s,t+1}{2}(e,n)=wc(1,scol.vGambleneg);
                
                % Reject trials(assembled plot only!)
                d_choiceinoptim{s,t+1}{1}(e,n)=sum(wc(wc(:, scol.Choice)==1, scol.vBestUnchosen)>wc(wc(:, scol.Choice)==1, scol.vChosen));
                d_choiceinoptim{s,t+1}{2}(e,n)=sum(wc(wc(:, scol.Choice)==2, scol.vBestUnchosen)>wc(wc(:, scol.Choice)==2, scol.vChosen));
                d_choiceinoptim{s,t+1}{3}(e,n)=sum(wc(wc(:, scol.Choice)==3, scol.vBestUnchosen)>wc(wc(:, scol.Choice)==3, scol.vChosen));
                d_choiceinoptim{s,t+1}{4}(e,n)=sum(wc(:, scol.vBestUnchosen)>wc(:, scol.vChosen));

                % Outcome
                d_outcome{s,t+1}{1}(e,n)= mean(wc(:, scol.OutcomeMagnitude));
                d_outcome{s,t+1}{2}(e,n)=  mean(wc(:, scol.OutcomeMagnitude)>0);
                d_outcome{s,t+1}{3}(e,n)=  mean(wc(:, scol.OutcomeMagnitude)==0);
                d_outcome{s,t+1}{4}(e,n)=  mean(wc(:, scol.OutcomeMagnitude)<0);
                
                
                
                wc=[];
            end
            
        end
    end
    
    
    % NOTE: EnvThreat from this column (subjdata{s,4} is already warped by j term (if
    % requested). Use EnvThreatOriginal to index. BE CAREFUL if you are
    % directly using EnvThreat at any point! 
    
    
    % Re-insert error trials and save
    subjdata{s,5}.errortrials(:, size(subjdata{s,5}.errortrials,2):size(ws.d,2))=nan;
    ws.d=sortrows([ws.d; subjdata{s,5}.errortrials], scol.Trialnum);
    d_trialstats{s,2}=ws.d(ws.d(:,scol.Task)==1,:);
    d_trialstats{s,3}=ws.d(ws.d(:,scol.Task)==2,:);
    d_trialstats{s,4}=ws.d;
    ws=[];
end

% Save + Set up 2nd level model folder 
request.logfol=[where.expt_fol filesep '2 Second level results' filesep '2 Behavioural model details' filesep];
request.modfol=[request.logfol details.modname];
if isdir(request.logfol)==0; mkdir(request.logfol); end; if isdir(request.modfol)==0; mkdir(request.modfol); end
% save([request.modfol filesep 'Model values (' date ')'],  'details', 'd_fits', 'd_vseenosee', 'd_par','d_trialvals',  'd_predchoice',    'd_vchoice', 'd_vbestchoice', 'd_pbestchoice', 'd_vchosenand', 'd_vchosenand', 'd_modalchoice_andval', 'd_trialstats', 'd_evgainloss', 'scol');

% error('done :)')

%% Ad-hoc Plots

% input('Plot?');

% close all hidden

% Plot what
d_plotwhat=d_vchoice;  details.titles={'[cF] V(Accept)';'[cF] V(Reject)';'[cF] V(Explore)';'[ct] V(Accept)';'[ct] V(Reject)';'[ct] V(Explore)'};
details.titles={'V(Accept)';'V(Reject)';'V(Explore)';'V(Accept)';'V(Reject)';'V(Explore)'};
d_plotwhat=d_predchoice;  details.titles={'[cF] predicted %Accept';'[cF] predicted %Reject';'[cF] predicted %Explore';'[ct] predicted %Accept';'[ct] predicted %Reject';'[ct] predicted %Explore'};
% d_plotwhat=d_bestchoice; details.titles={'[cF] Best choice';'[ct] Best choice'};
% d_plotwhat=d_vbestchoice;  details.titles={'[cF] vBest';'[cF] vSecondBest'; '[cF] vWorst'; '[ct] vBest';'[ct] vSecondBest'; '[ct] vWorst'};
% d_plotwhat=d_pbestchoice; details.titles={'[cF] pChooseBest';'[cF] pChoose2ndBest';'[cF] pChooseWorst'; '[ct] pChooseBest';'[ct] pChoose2ndBest'; '[ct] pChooseWorst';};
% d_plotwhat=d_vchosenand;  details.titles={'[cF] vChosen';'[cF] vBestUnchosen';'[cF] vCho-vBU (trialwise)'; '[ct] vChosen';'[ct] vBestUnchosen';'[ct] vCho-vBU (trialwise)'};
% d_plotwhat=d_vgamble;  details.titles={'[cF] V(Gamble)';'[ct] V(Gamble)'};
% d_plotwhat=d_vgamble_posneg;  details.titles={'[cF] Subjective EV+';  '[cF] Subjective EV-'; '[ct] Subjective EV+'; '[ct] Subjective EV-';};
% d_plotwhat=d_vbestunchosen_posneg;  details.titles={'[cF] vBestUncho+';'[cF] vBestUncho-';  '[ct] vBestUncho+';'[ct] vBestUncho-';};
% d_plotwhat=d_vBestUnchosenXUnchosChoice; details.titles={'[cF] vBestUncho_A';'[cF] vBestUncho_R';  '[cF] vBestUncho_E'; '[ct] vBestUncho_A';'[ct] vBestUncho_R';  '[ct] vBestUncho_E';  };
% d_plotwhat=d_modalchoice_andval;  details.titles={'[cF] Modal choice';'[cF] V(Modal choice)';  '[ct] Modal choice';'[ct] V(Modal choice)'};
% d_plotwhat=d_vseenosee;   details.titles={'[cF] V(See)';  '[cF] V(No See)'; '[ct] V(See)'; '[ct] V(No See)';};
% d_plotwhat=d_evgainloss;  details.titles={ '[cF] EVGain';'[cF] EVLoss';'[cF] EVConflict'; '[ct] EVGain';'[ct] EVaLoss';'[ct] EVConflict'; };
% d_plotwhat=d_choiceinoptim; details.titles={'[cF] No. Acc, vCho < vBU';'[cF] No. Rej, vCho < vBU';   '[cF] No. Expl, vCho < vBU';   '[cF] No. trials vCho < vBU'; '[ct] No. Acc, vCho < vBU';'[ct] No. Rej, vCho < vBU'; '[ct] No. Expl, vCho < vBU'; '[ct] No. trials vCho < vBU'; };
% d_plotwhat=d_outcome; details.titles={'[cF] Outcome';'[cF] pPos';'[cF] pNeg'; '[ct] Outcome';'[ct] pPos'; '[ct] pNeg';};
% d_plotwhat=d_prob; details.titles={'[cF] prior p';'[cF] pos p';'[cF] Expected Info Gain KLd'; '[cF] Expected Entropy Change';  '[ct] prior p';'[ct] pos p';'[ct] Expected Info Gain';  '[ct] Expected Entropy Change'; };

% details.titles={'V(Chosen)';'V(Best Unchosen)';'V(Worst)';'V(Chosen)';'V(Best Unchosen)';'V(Worst)';};



% ##############
% details.titles={'V(See)';  'V(No See)'; 'V(See)'; 'V(No See)';};
% details.titles={'Exp';'Exp';'Exp';'Control';'Control';'Control'};
% details.titles={'V(Chosen)';'V(Best Unchosen)';'[cF] vWorstUnchosen'; 'V(Chosen)';'V(Best Unchosen)';'[ct] vWorstUnchosen'};
% details.titles={'V(Best Unchosen) +'; 'V(Best Unchosen) -'; 'V(Best Unchosen) +'; 'V(Best Unchosen) -'}'

% Compile
d_plotmatrix=[[details.subjects; {'Group'}] cell(details.n_subjs+1, length(details.titles))]; % subject, [cF] vA, vR, vE, [ct] vA, vR, vE
nc=length(d_plotwhat{1,2});  d_plotmatrix(details.n_subjs+1,2:end)=repmat({ones(6,6)},1,2*nc);
for t=1:2
    for c=1:nc
        for s=1:details.n_subjs
            d_plotmatrix{s, (t-1)*nc+c+1}=fliplr(d_plotwhat{s,t+1}{c}')';
            d_plotmatrix{details.n_subjs+1, (t-1)*nc+c+1}=d_plotmatrix{details.n_subjs+1, (t-1)*nc+c+1}+d_plotmatrix{s, (t-1)*nc+c+1};
        end
        d_plotmatrix{details.n_subjs+1, (t-1)*nc+c+1}=d_plotmatrix{details.n_subjs+1, (t-1)*nc+c+1}/details.n_subjs;
    end
end

% Plot group
f.subplotcols=(2*nc)+1;  f.figwidth= 1200;  f.figheight=600; f.subplot_VerHorz=[0.16 0.055];f.fig_BotTop=[0.1 0.1];  f.fig_LeftRight=[0.06 0.03]; 
f.f=figure('Color',[1 1 1], 'Name', 'Group value plots', 'NumberTitle', 'off', 'Position',[50,70,f.figwidth,f.figheight]);  s=details.n_subjs  +1;
subtightplot(1,  f.subplotcols,  1, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); text(0.8, 0.5, d_plotmatrix{s,1}); axis off
f.FontSize=40; f.FontName='PT Sans';k=2; 
for c=1:2*nc
%     subtightplot(1,  f.subplotcols,  k, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1; 
    if c==1; k=1; end; subtightplot(2,  nc,  k, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1; 
    imagesc(d_plotmatrix{s, c+1}); axis square;
    colorbar;  

% %     %         % Standardize ranges within task
%     if c<=nc; caxis([min(cell2mat(cellfun(@(x)min(min(x)), d_plotmatrix(s,2:2+nc-1), 'UniformOutput', 0)))     max(cell2mat(cellfun(@(x)max(max(x)), d_plotmatrix(s,2:2+nc-1), 'UniformOutput', 0)))])
%     else caxis([min(cell2mat(cellfun(@(x)min(min(x)), d_plotmatrix(s,nc+2:end), 'UniformOutput', 0)))     max(cell2mat(cellfun(@(x)max(max(x)), d_plotmatrix(s, nc+2:end), 'UniformOutput', 0)))])
%     end

    % Manual caxis
%     if c==1 || c==2
%     caxis([0 12]); disp('Manual axes')
%     end
    
    
%     title(details.titles{c}, 'FontSize', f.FontSize,'FontName', f.FontName);
    set(gca, 'FontSize', f.FontSize,'FontName', f.FontName);
    
    
    ylabel('EnvThreat');  xlabel('No. Tokens');  set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/6' '2/6' '3/6' '4/6' '5/6' '6/6'}),'XTick',1:6, 'XTickLabel',2:2:2*6)
    
end


% Optimal choice (group level)
if sum(sum(d_plotwhat{1,2}{1}-d_vchoice{1,2}{1}))==0; % Only for vChoices
    f.figwidth= 1000;  f.figheight=400; f.subplot_VerHorz=[0.005 0.055];f.fig_BotTop=[0.08 0.035];  f.fig_LeftRight=[0.05 0.05]; f.subplotcols=7;
    f.f=figure('Name', 'Optimal choice (group level) ', 'NumberTitle', 'off', 'Position',[50,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]); s=details.n_subjs  +1;
    subtightplot(1,  f.subplotcols,  1, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); axis off; k=2; pr={'cF';'ct'};
    text(0.5, 0.5, d_plotmatrix{s,1}); 
    for t=1:2
        tc=(t-1)*3+1;
        subtightplot(1,  f.subplotcols,  k, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagesc(d_plotmatrix{s, 1+tc}> max(d_plotmatrix{s, 2+tc}, d_plotmatrix{s, 3+tc})); colorbar; k=k+1; axis square;  title(['[' pr{t} '] Optimal Accept'])
        subtightplot(1,  f.subplotcols,  k, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagesc(d_plotmatrix{s, 2+tc}> max(d_plotmatrix{s, 1+tc}, d_plotmatrix{s, 3+tc})); colorbar; k=k+1; axis square; title(['[' pr{t} '] Optimal Reject'])
        subtightplot(1,  f.subplotcols,  k, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        imagesc(d_plotmatrix{s, 3+tc}> max(d_plotmatrix{s, 1+tc}, d_plotmatrix{s, 2+tc})); colorbar; k=k+1; axis square; title(['[' pr{t} '] Optimal Explore'])
    end
end

% Plot individuals
plotindiv=0;
if plotindiv
    f.figwidth= 1000;  f.figheight=1200; f.subplot_VerHorz=[0.005 0.055];f.fig_BotTop=[0.001 0.035];  f.fig_LeftRight=[0.05 0.01]; f.subplotcols=7;
    f.f=figure('Name', 'Individual value plots', 'NumberTitle', 'off', 'Position',[50,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    for s=1:  1 % details.n_subjs
        subtightplot(details.n_subjs+1,  f.subplotcols,  (s-1)*f.subplotcols+1, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        text(1, 0.6, d_plotmatrix{s,1}); axis off
        
        for c=1:2*nc
            subtightplot(details.n_subjs+1,  f.subplotcols,  (s-1)*f.subplotcols+1+c, f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagesc(d_plotmatrix{s, c+1});
            colorbar; axis off; axis square;
            
            % Standardize ranges within task
            if c<=nc; caxis([min(cell2mat(cellfun(@(x)min(min(x)), d_plotmatrix(s,2:2+nc-1), 'UniformOutput', 0)))     max(cell2mat(cellfun(@(x)max(max(x)), d_plotmatrix(s,2:2+nc-1), 'UniformOutput', 0)))])
            else caxis([min(cell2mat(cellfun(@(x)min(min(x)), d_plotmatrix(s,nc+2:end), 'UniformOutput', 0)))     max(cell2mat(cellfun(@(x)max(max(x)), d_plotmatrix(s, nc+2:end), 'UniformOutput', 0)))])
            end
            
            % Manual caxis
            %         caxis([0.4 1])
            
            
            if s==1; title(details.titles{c}); end
        end
    end
end

error

%% Plot static quantities (i.e. fpar_conflict etc)



col=rmfield(scol, {'OutcomeMagnitude','OutcomeMean','OutcomeVariance','OutcomePres','VExplore'}); col.StanDev=10;  col.vMeanVar=11; col.BinomVar=12;
col.EVGain=38;
col.EVLoss=39;
col.EVConflict=40;
col.Misc=41;
d_par=nan(36,3); d_par(:, [col.EnvThreat col.NTokens])=[repmat(1:6, 1,6)'/6 sortrows(repmat(1:6, 1,6)')*2]; 
%

dcf=fpar_conflict(d_par, col);  % Go turn on plotting here if you want the scripts 
dct=fpar_control(d_par, col);
[d_cfvarcor dp_cfvarcor]=corr(dcf(:, [col.EnvThreat col.NTokens col.pLoss col.Entropy col.EntropyNTok col.EV]));
d_cfvarcor= d_cfvarcor(:) ; d_cfvarcor(dp_cfvarcor(:)>0.5)=nan;  d_cfvarcor= reshape(d_cfvarcor, 6,6);
[d_ctvarcor dp_ctvarcor]=corr(dct(:, [col.EnvThreat col.NTokens col.pLoss col.Entropy col.EntropyNTok col.EV]));
d_ctvarcor= d_ctvarcor(:) ; d_ctvarcor(dp_ctvarcor(:)>0.5)=nan;  d_ctvarcor= reshape(d_ctvarcor, 6,6);
%
[d_varcor dp_varcor]=corr([dcf(:, [col.EnvThreat col.NTokens col.pLoss col.Entropy col.EntropyNTok col.EV]) dct(:, col.EV)]);
d_varcor= d_varcor(:) ; d_varcor(dp_varcor(:)>0.5)=nan;  d_varcor= reshape(d_varcor, 7,7);
figure('color','w'), f.FontName='PT Sans'; f.FontSize=20;  f.FontSize_title=25;
imagescnan(d_varcor, 'nancolor', [0.3 0.3 0.3]), axis square, colorbar
% set(gca, 'xticklabel', {'EnvThreat';'NTokens';'pLoss';'Uncertainty';'Value scaled uncertainty';'EV (Exp)';'EV (Ctrl)'},'yticklabel', {'EnvThreat';'NTokens';'pLoss';'Entropy';'Value scaled uncertainty';'EV (Exp)';'EV (Ctrl)'})
set(gca, 'yticklabel', {'1 EnvThreat';'2 No. Tokens';'3 p(ActBomb)';'4 Entropy';'5 Value-scaled uncertainty';'6 EV (Exp)';'7 EV (Ctrl)'},'FontName', f.FontName,  'FontSize', f.FontSize); 
set(gca,'FontName', f.FontName,  'FontSize', f.FontSize); 
title('Correlation of variables in task space','FontName', f.FontName,  'FontSize', f.FontSize_title); 
% legend('EnvThreat','No. Tokens', 'p(ActBomb)','Uncertainty','Value-scaled uncertainty','EV (Exp)','EV (Ctrl)')
% rotateXLabels( gca(), 45)
% set(gca,'XAxisLocation', 'bottom')


%% Check onsets files (saved in second level modelling folder)

docheckonsets=0; 
if docheckonsets
    c=load([request.modfol filesep 'CheckETwarp_s1.mat']);
    request.check = c.request.check;    s=1;
    
    ws.cfo  = c.d_check{3};  ws.cto  = c.d_check{4};   % From onsets
    ws.cft= d_trialstats{s,2};  ws.ctt= d_trialstats{s,3};    % From this script
    ws.cft_et = sortrows(unique( ws.cft(:, scol.EnvThreat))); ws.ctt_et = sortrows(unique( ws.ctt(:, scol.EnvThreat)));
    
    
    % Assemble grid for current variables (i.e. from this script)
    for e=1:6
        for n=1:6
            ws.cfc= ws.cft(ws.cft(:, scol.EnvThreat) ==  ws.cft_et (e) & ws.cft(:, scol.NTokens) ==  n*2, :);
            ws.ctc= ws.ctt(ws.ctt(:, scol.EnvThreat) ==  ws.ctt_et (e) & ws.ctt(:, scol.NTokens) ==  n*2, :);
            
            for v=1:length(request.check )
                eval(['ws.cf{v}(7-e, n)= unique(ws.cfc(:, col.' request.check{v} '));'])
                eval(['ws.ct{v}(7-e, n)= unique(ws.ctc(:, col.' request.check{v} '));'])
            end
        end
    end
    
    % Check current variables against onsets 
    for v=1:length(request.check )
        if sum(abs(ws.cf{v}(:) - ws.cfo{v}(:)))==0
            disp([request.check{v} ' warped ok']);
        else disp([request.check{v} ' warped WRONG!!']);
        end
    end
    
    ws=[];
end




%% Parameters correlated?

% Trade off between parameters 
whichmods=[details.cf_moddetails(1)  details.ct_moddetails(1)];
modpars={details.cf_moddetails{3} details.ct_moddetails{3}};
d_parcor= [{d_fits.cf.r_res{1,2}(:, 4:end)} {d_fits.ct.r_res{1,2}(:, 4:end)}]; 
for t=1:2
    d_r{t}=nan(length(modpars{t}),length(modpars{t}));
    for p1=1:length(modpars{t})
        for p2=1:length(modpars{t})
            [r p]=corr(d_parcor{t}(:, p1), d_parcor{t}(:, p2));
            
            if p<0.05  % SIGNIFICANT
                stat=1;
                stat=r;
            elseif p<0.1
                stat=0.5;
                stat=r;
            else stat=nan;
            end
            
            
            
            % Print
            d_r{t}(p1,p2)=stat;
        end
        
    end
end
figure('color','w')
for t=1:2
subplot(1,2,t)
% imagesc(d_r{t}); 
imagescnan(d_r{t},'nancolor',[0.9 0.9 0.9]);  axis square; colorbar


set(gca, 'xtick', 1:length( modpars{t}), 'xticklabel', modpars{t}, 'ytick', 1:length( modpars{t}), 'yticklabel', modpars{t}, 'FontSize',15)
title(whichmods{t},'FontSize',20)
end

% Parameters correlated with other information? 
for o=1:1
%     STOPPED HERE!!!
%     where.behprofile='/Users/EleanorL/Dropbox/SCRIPPS/5 Explore fMRI/1 Behavioural data/Group behaviour files';
%     [n t r]=xlsread([where.behprofile filesep 'Behavioural profile for correlation.xlsx']);
    
    % State anxiety, trait, Bis, Bas
    subpersonality=[7,14,23,44;0,6,20,42;3,8,16,43;9,15,24,41;2,8,18,39;1,4,22,47;6,8,14,47;4,8,19,43;1,8,16,39;30,23,22,39;13,26,26,32;12,16,19,37;16,10,15,31;0,7,18,48;4,16,18,43;11,9,12,44;15,20,13,33;15,21,27,38;3,10,24,36;13,20,19,41;]; 
    spers={'state','trait','bis','bas'};
    
    subpar=d_parcor{1}(:, 3:5);
    subparnames=modpars{1}(3:5);
    for p=1:size(subpar,2)
        for pers=1:size(subpersonality,2)
            [r pp]=corr(subpar(:,p),  subpersonality(:, pers));
            
            disp([subparnames{p} ' and ' spers{pers} ': r=' num2str(r) ',  p=' num2str(pp)]);
        end
    end
    
    
end



%% No. Trials with +/- vBU?

d_vbu{1}=[cell(details.n_subjs,4); repmat({zeros(6,6)},1,4)];  d_vbu{2}=d_vbu{1}; % Rej_vBu+,  Rej_vBu-, NonRej_vBU+, NonRej_vBU-
d_vbu{3}=nan(details.n_subjs,4);  d_vbu{4}=d_vbu{3};
for t=1:2
    for s=1:details.n_subjs
        ws.tt=d_trialstats{s,t+1};
        
        for e=1:6
            ee=7-e;
            for n=1:6
                ws.t= ws.tt(ws.tt(:, scol.EnvThreat)*6==e & ws.tt(:, scol.NTokens)/2==n, :);
                ws.tpos=ws.t(ws.t(:, scol.vBestUnchosen)>=0, :);
                ws.tneg=ws.t(ws.t(:, scol.vBestUnchosen)<0, :);
                ws.tpos_rej=ws.tpos(ws.tpos(:, scol.Choice)==2, :);
                ws.tpos_nonrej=ws.tpos(ws.tpos(:, scol.Choice)~=2, :);
                ws.tneg_rej=ws.tneg(ws.tneg(:, scol.Choice)==2, :);
                ws.tneg_nonrej=ws.tneg(ws.tneg(:, scol.Choice)~=2, :);
                d_vbu{t}{s,1}(ee, n)=size(ws.tpos_rej,1);
                d_vbu{t}{s,3}(ee, n)=size(ws.tpos_nonrej,1);
                d_vbu{t}{s,2}(ee, n)=size(ws.tneg_rej,1);
                d_vbu{t}{s,4}(ee, n)=size(ws.tneg_nonrej,1);
            end
        end
        
        d_vbu{2+t}(s,1)= sum(ws.tt(:, scol.vBestUnchosen)>=0 & ws.tt(:, scol.Choice)==2) ;
        d_vbu{2+t}(s,2)= sum(ws.tt(:, scol.vBestUnchosen)<0 & ws.tt(:, scol.Choice)==2);
        d_vbu{2+t}(s,3)= sum(ws.tt(:, scol.vBestUnchosen)>=0 & ws.tt(:, scol.Choice)~=2) ;
        d_vbu{2+t}(s,4)= sum(ws.tt(:, scol.vBestUnchosen)<0 & ws.tt(:, scol.Choice)~=2);
        
        
        for c=1:4
            d_vbu{t}{end,c}=d_vbu{t}{end,c}+d_vbu{t}{s,c};
        end
    end
    for c=1:4
            d_vbu{t}{end,c}= d_vbu{t}{end,c}./details.n_subjs;
    end
end
close all hidden, figure('color','w'); k=1; titles={'[vBU+] Reject'; '[vBU-] Reject'; '[vBU+] NonReject';  '[vBU-] NonReject'};  titles=[titles; titles]; f.FontSize=20;
for t=1:2
    for c=1:4
        subplot(3,4,k);
        imagesc(d_vbu{t}{end,c}); colorbar, axis square; 
        title(titles{k}, 'FontSize', f.FontSize); 
        caxis([0 15])
        xlabel('NTokens'); ylabel('EnvThreat'); axis off; k=k+1;
    end
end



figure('color','w')


% subplot(3,4,k:k+1);  
subplot(1,2,1)
ws.vars= reshape(std(d_vbu{3})./sqrt(details.n_subjs), 2,2)';  ws.means= reshape(mean(d_vbu{3}) ,2,2)';
barwitherr(ws.vars, ws.means);
set(gca, 'xticklabel', {'Reject' 'Accept/Explore'}, 'FontSize', f.FontSize, 'FontName', f.FontName);  legend({'+ Counterfactual value';'- Counterfactual value'}, 'FontSize', f.FontSize, 'FontName', f.FontName); k=k+2; 
ylabel('No. of trials', 'FontSize', f.FontSize, 'FontName', f.FontName)
% title('[Exp] No. trials split by vBU valence and Rejecting', 'FontSize', f.FontSize); 
title('Experimental condition', 'FontSize', f.FontSize, 'FontName', f.FontName); 
ws.vars= reshape(std(d_vbu{4})./sqrt(details.n_subjs), 2,2)';  ws.means= reshape(mean(d_vbu{4}) ,2,2)';  % 
% subplot(3,4,k:k+1);  
subplot(1,2,2)
barwitherr(ws.vars, ws.means);
set(gca, 'xticklabel', {'Reject' 'Accept/Explore'}, 'FontSize', f.FontSize, 'FontName', f.FontName);  
legend({'+ Counterfactual value';'- Counterfactual value'}, 'FontSize', f.FontSize, 'FontName', f.FontName)
ylabel('No. of trials', 'FontSize', f.FontSize, 'FontName', f.FontName)
title('Control condition', 'FontSize', f.FontSize, 'FontName', f.FontName); 
% title('[Ctrl] No. trials split by vBU valence and Rejecting', 'FontSize', f.FontSize); 


%% No trials w inoptim choice 

d_ci=nan(details.n_subjs,8); ch={'Acc';'Rej';'Expl';'All'};    % cF then ct; acc, rej, expl, all
for s=1:details.n_subjs, k=1; 
    for t=1:2
        for c=1:4
            d_ci(s,k)= mean(d_choiceinoptim{s,t+1}{c}(:)); k=k+1; 
        end
    end
end

% Compare cF and ct
for c=1:4 
    [wc.h wc.p wc.ci wc.stats]= ttest(d_ci(:,c)-d_ci(:, c+4));
    disp([ ch{c} ' : t(' num2str(wc.stats.df) ')=' num2str(wc.stats.tstat,3) '  ,p=' num2str(wc.p,3)]);
end

figure('color','w'), f.FontSize=25; f.FontName='PT Sans';    
barwitherr(reshape(std(d_ci)./sqrt(details.n_subjs),4,2), reshape(mean(d_ci),4,2))
xlabel('Choice','FontSize', f.FontSize, 'FontName', f.FontName), ylabel('Mean no. Trials','FontSize', f.FontSize, 'FontName', f.FontName)
title('Non-optimal choice, V(Chosen)<V(Best Unchosen)','FontSize', f.FontSize, 'FontName', f.FontName)
legend('Exp', 'Ctrl')
set(gca,'FontSize', f.FontSize, 'FontName', f.FontName)
set(gca,'xticklabel',ch)

% 
% ,'FontSize', f.FontSize, 'FontName', f.FontName)
% 
% ,'FontSize', f.FontSize, 'FontName', f.FontName)

