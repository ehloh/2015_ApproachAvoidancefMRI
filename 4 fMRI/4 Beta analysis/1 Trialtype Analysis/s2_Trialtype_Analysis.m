% Analyse Trialtype betas
clear all;close all hidden; clc
% clear all;clc

% Request
request.do_glm=0;
request.do_glm_stepwise=0;
request.do_pca=0;
request.var_conweights=0;
request.plot_tstat=1; 
request.beh_correlations=0;  % Manualy specify in whichever section
%
request.TaskspaceCanonical=1;  % Canonical (vs subject-specified-warped) task space?
request.meancentre=1;
%
request.specificsubjects= {};

request.FLmodel='t1_1_Trialtype';  % ##### WHICH FL MODEL #########
% request.FLmodel='t2_1_TrialtypeNc';

% Where are the ROI images/beta file  ###### WHICH ROIS ##########
% request.roibatch='c13 battery'; 
% request.roibatch='c13 amyg'; 
% request.roibatch='Subfields 95'; 
request.roibatch='Anat HPC Amyg'; 
% request.roibatch='c14 battery'; 
% request.roibatch='v28gChoice battery';
% request.roibatch='c14 HPC'; 
% request.roibatch='v3g HPC'; 


% Behavioural modelling results ###### BEH MODELLING FILES ##########
where.modres='3 Hierarchical';
request.behres_cf='res_hierarfitmodels_cF (16-Apr-2015) bpji8_L10968';
request.behres_ct='res_hierarfitmodels_ct (16-Apr-2015) bpji11_L11293';

for o1=1:1 % General settings and specifications 
    
    % Where
    log.AntsType='_Basic';  path(pathdef);  
    w.w=pwd; if strcmp(w.w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/1 Explore fMRI'; where.experiment_folder='/Users/EleanorL/Desktop/2 EXPLORE fMRI data';  where.data_brain=[where.experiment_folder filesep '1 Brain data'];  where.modellingdata=     '/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis/4 Fit computational models/2 Analysis inputs';   where.modfol='/Users/EleanorL/Dropbox/SCRIPPS/2 Explore experiment/3 Analysis/4 Fit computational models';
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI'; where.experiment_folder='G:\2 [Explore]';   where.modfol='C:\Users\e.loh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis\4 Fit computational models'; 
        where.data_brain='G:\2 [Explore]\1 Brain data'; addpath(where.where);  where.modellingdata='C:\Users\eloh\Dropbox\SCRIPPS\2 Explore experiment\3 Analysis\4 Fit computational models\2 Analysis inputs';   end;  addpath(where.where); addpath(where.modfol)
    
    
    
    
    if isdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep request.FLmodel log.AntsType])==0; mkdir([where.experiment_folder filesep '2 Second level results s4Ants' filesep request.FLmodel log.AntsType]); end
    request.SLfol=[where.experiment_folder filesep '2 Second level results s4Ants' filesep request.FLmodel log.AntsType ]; if isdir(request.SLfol)==0; mkdir(request.SLfol); end
    request.FLfol=['m_' request.FLmodel log.AntsType  ' Contrasted'];
    request.where_betas= [request.SLfol filesep request.roibatch];
    addpath([where.modfol filesep '1 Value functions'])
    where.modvals4fmri  = [where.experiment_folder filesep '2 Second level results' filesep '2 Behavioural model details' filesep];
    
    
    
%     C:\Users\e.loh\Dropbox\SCRIPPS\1 Explore fMRI\7 Beta analysis\1 Trialtype Analysis
    
    
    
    
    
    
    
    % Behavioural modelling results 
    request.behfilename='All data (09-May-2014)'; % From modelling folder
    request.task={'cF' 'ct'}; 
    request.allvars={'Subject' 's'; 'EnvThreat' 'ET';'NTokens' 'N'; 'pLoss' 'pL'; 'Entropy' 'Un'; 'EntropyNTok' 'UnN';  'EV' 'EV'; 'PavConflict' 'PCon'; 'EVConflict' 'EVCon';'EnvThreatOrig' 'ETo';  'vChosen' 'vCh'; 'vBestUnchosen' 'vBU';  'vBestUnchosen_neg' 'vBUn'; 'vBestUnchosen_pos' 'vBUp';   };
    addpath([where.modfol filesep '1 Value functions' ]);
    addpath(genpath([where.modfol filesep '1 Value functions' ]));
    addpath([where.modfol filesep '1 Value functions'  filesep 'b']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bp']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bm']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bjm']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bi']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bpi']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bmi']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bpm']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bpmi']);
    addpath([where.modfol filesep '1 Value functions'  filesep 'bpjm']);
    
    
    % Subjects
    w=load([where.data_brain filesep 'datalogexplore_allsubs.mat']); log.datalog=w.datalog;
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, request.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    request.n_subjs=log.n_subjs ;
    
    %
    if isempty(strfind(request.FLmodel, 'Chunk'))==1; log.n_levels=6;  request.trialchar={'t';[]};
    else log.n_levels=3;  request.trialchar={'e';'n'};
    end
    if strcmp(request.FLmodel(1), 'f')==1; log.n_choices=3; else log.n_choices=1; end, request.n_choices=log.n_choices;
    request.subchecks_forflex=0; 
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(request.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end; disp(' ');
    disp(['             First level model:              ' request.FLmodel])
    disp(['             No. of EnvThreat/NTokens levels:       ' num2str(log.n_levels)]); 
    disp(['             Cells divided by how many choices:     ' num2str(log.n_choices)]); disp(' ')
    disp(['             FL folder:      ' request.FLfol])
    disp(['             SL folder:      ' request.SLfol]); disp(' ');
%     input('Hit Enter to start      ')
    disp('=======================================================')
end

%% (0) Settings + Load all materials 
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per imagesc plot)
% 
request.ROIs={
%     'HPC_aL_sc'; 'HPC_aR_sc' % c13g battery ################
%     'HPC_aR_stc'; 'HPC_aL_stc'
    
%     'Striatum_R';'DLPFC_L';'DLPFC_R';
%     'HPC_aL_stc'; 'HPC_aR_stc'; 'HPC_aL_sc'; 'HPC_aR_sc'   % c13 battery ################
%     'BA10'; 'Striatum_R'
%     'DLPFC_BA10_L'; 'DLPFC_BA10_R'    
'rHPC_aL';'rHPC_aR'; 'rAmygdala_L';'rAmygdala_R';       % Anat ########5
%     'HPC_aL_atc';'HPC_aR_atc';  % c14 battery ################
%     'HPC_aL_ac';'HPC_aR_ac';
%     'DLPFC_BA10_L';'DLPFC_BA10_R';
%     'HPC_aL_atc';'HPC_aR_atc';  % v28g battery ################
%     'HPC_aL_apc';'HPC_aR_apc';
%     'DLPFC_L_c';'DLPFC_R_c';
%     'rlPFC_L_c';'Parietal_c';
%     'HPC_aL_tv001'; 'HPC_aR_tv001' % v3g battery ################
%     'HPC_aL_mev05fwe'; 'HPC_aR_mev05fwe'   
% 'HPC_L_s001atc'; 'HPC_R_s001atc'    %     c14 HPC (Masking) ################
% 'HPC_L_s001ac'; 'HPC_R_s001ac'
% 'HPC_L_s01ac'; 'HPC_R_s01ac'
% 'HPC_L_s01atc'; 'HPC_R_s01atc'
% 'HPC_L_s05atc'; 'HPC_R_s05atc'
};  

% (1) Load beta matrix
[log d_roimatrix_cF d_roimatrix_ct] = f_LoadTrialtypeBetamatrix(request, log);

% (2) Load behavioural files
for o1=1:1 % Columns 
    col.EnvThreat=1;
    col.NTokens=2;
    col.pLoss=3;
    col.Entropy=4;
    col.EntropyNTok=5;
    col.EV=6;
    col.EVGain=7;
    col.EVLoss=8;
    col.EVConflict=9;
    col.PavConflict=10;
    %
    col.Subject=11;
    col.pLossAversive=12;
    col.EnvThreatOrig=13;
    col.NTokensOrig=14;
    col.b=15;
    col.simb=16;  % simulated b (from glm)
    col.Task=17;  % For beh modelling scripts
    col.Trialnum=18;
    col.VExplore=19;
    col.Choice=20;
    col.Misc=21;
    
    % Value columns
    col.vChosen=25;
    col.vBestUnchosen= col.vChosen+1;
    col.vWorstUnchosen= col.vChosen+2;
    col.vBestUnchosen_pos= col.vChosen+3;
    col.vBestUnchosen_neg= col.vChosen+4;
    col.BestUnchosen_Is= col.vChosen+5;
    col.vBestUnchosen_UnchoChoice(1)= col.vChosen+6;
    col.vBestUnchosen_UnchoChoice(2)= col.vChosen+7;
    col.vBestUnchosen_UnchoChoice(3)= col.vChosen+8;
    col.vGamble= col.vChosen+11;
    col.vGamblepos= col.vChosen+12;
    col.vGambleneg= col.vChosen+13;
    col.Modalchoice=col.vChosen+14;  % If draw
    col.vModalchoice=col.vChosen+15;
    col.vSee=col.vChosen+16;
    col.vNoSee=col.vChosen+17;
    col.vAccept=col.vChosen+18;
    col.vReject=col.vChosen+19;
    col.vExplore=col.vChosen+20;
    col.vBestNonmodal =col.vChosen+21;
end
w.cf=load([where.modfol filesep '2 Analysis inputs' filesep where.modres filesep request.behres_cf '.mat']);
w.ct=load([where.modfol filesep '2 Analysis inputs' filesep where.modres filesep request.behres_ct '.mat']);
disp(['Models used: '  w.cf.r_res{1,1} ' for cF, ' w.ct.r_res{1,1} ' for ct'])
log.modres.cf=[w.cf.r_res{1,1}, w.cf.details.models(strcmp(w.cf.details.models(:,1),  w.cf.r_res{1,1}),3) {w.cf.r_res{1,2}(:,4:end)}];
log.modres.ct=[w.ct.r_res{1,1}, w.ct.details.models(strcmp(w.ct.details.models(:,1),  w.ct.r_res{1,1}),3) {w.ct.r_res{1,2}(:,4:end)}];
if isempty(request.specificsubjects)==0
    log.modres.cf{3} = log.modres.cf{3}(cell2mat(cellfun(@(x)find(strcmp(w.cf.details.subjects, x)), request.specificsubjects, 'UniformOutput',0)), :);
    log.modres.ct{3}=  log.modres.ct{3}(cell2mat(cellfun(@(x)find(strcmp(w.ct.details.subjects, x)), request.specificsubjects, 'UniformOutput',0)), :);
end
if request.TaskspaceCanonical, [ d_design] = fpar_taskspacevar( request, col,  []); else [ d_design] = fpar_taskspacevar( request, col, {log.modres.cf log.modres.ct}); end
d_design = [[log.subjects;{'GroupSpace';'GroupPar'}] d_design];
d_beh=load([where.modvals4fmri log.modres.cf{1}(1:strfind(log.modres.cf{1}, '_')-1) log.modres.ct{1}(1:strfind(log.modres.ct{1}, '_')-1) filesep spm_select('List', [where.modvals4fmri log.modres.cf{1}(1:strfind(log.modres.cf{1}, '_')-1) log.modres.ct{1}(1:strfind(log.modres.ct{1}, '_')-1)] , 'Model values.*.mat')]);
for e=1:6  % Load pre-calculated values 
    for n=1:6
        for t=1:2
            wc.n=find(d_design{1,t+1}(:, col.EnvThreatOrig).*6==e & d_design{1,t+1}(:, col.NTokens)./2==n);
            for s=1:log.n_subjs
                ss= strcmp( d_beh.details.subjects, log.subjects{s});
                
                d_design{s,t+1}(wc.n,  col.vChosen) = d_beh.d_vchosenand{ss,t+1}{1}(e,n);
                d_design{s,t+1}(wc.n,  col.vBestUnchosen) = d_beh.d_vchosenand{ss,t+1}{2}(e,n);
                d_design{s,t+1}(wc.n,  col.vBestUnchosen) = d_beh.d_vchosenand{ss,t+1}{2}(e,n); 
                
                % Add more quantities here 
            end
        end
    end
end
for t=1:2 % Manually specify variables 
    for s=1:log.n_subjs
        d_design{s,t+1}(:, [col.vBestUnchosen_pos col.vBestUnchosen_neg])=0;
        d_design{s,t+1}(  find(d_design{s,t+1}(:, col.vBestUnchosen)>0), col.vBestUnchosen_pos) = d_design{s,t+1}(  find(d_design{s,t+1}(:, col.vBestUnchosen)>0), col.vBestUnchosen);
        d_design{s,t+1}(  find(d_design{s,t+1}(:, col.vBestUnchosen)<0), col.vBestUnchosen_neg) = d_design{s,t+1}(  find(d_design{s,t+1}(:, col.vBestUnchosen)<0), col.vBestUnchosen);
    end
end
if request.beh_correlations==1;
    [wb.n wb.t wb.r]=xlsread([where.where filesep '1 Behavioural data' filesep 'Group behaviour files' filesep 'Behavioural profile for correlation.xlsx'], 'BehInhibExplore');
    d_behscores=[cellstr(wb.t(1:size(wb.n,1)+1,1))  [wb.t(1,2:end) ;  num2cell(wb.n)]];
    if isempty(request.specificsubjects)==0, d_behscores =[d_behscores(1,:);  d_behscores(cellfun(@(x)find(strcmp(d_behscores(:,1), x)), request.specificsubjects), :)]; end
end

% (3) Compile betas and task-space variables
request.rescaleIVs=1;
d_betas=[log.rois repmat({[]}, log.n_rois,2)];
for s=1:log.n_subjs
    % IVs
    for v=1:size(request.allvars,1)
        eval(['ws.dcf(:,col.' request.allvars{v,1} ')=d_design{s,1+1}(:, col.' request.allvars{v,1} ');'])
        eval(['ws.dct(:,col.' request.allvars{v,1} ')=d_design{s,1+2}(:, col.' request.allvars{v,1} ');'])
    end
    ws.dcf(:,col.Subject)=s; ws.dct(:,col.Subject)=s;
    ws.dcf(:,col.NTokensOrig)=ws.dcf(:,col.NTokens); ws.dct(:,col.NTokensOrig)=ws.dct(:,col.NTokens);
     
    % Rescale IVs
    if request.rescaleIVs
        for v=1:size(request.allvars,1)
            if sum(strcmp(request.allvars{v}, {'Subject'; 'EnvThreatOrig'; 'NTokensOrig';})) ==0
                eval(['ws.v=ws.dcf(:, col.' request.allvars{v,1}  ');'])
                ws.v=ws.v - min(ws.v); % Adjust floor
                ws.v= ws.v./ max(ws.v);  % Scale by ceiling
                eval(['ws.dcf(:, col.' request.allvars{v,1}  ')=ws.v;'])
                %             disp([request.allvars{v} '  - cF range: ' num2str(min(ws.v)) '    -    '  num2str(max(ws.v))])
                
                if sum(strcmp(request.allvars{v}, {'EVConflict', 'PavConflict', 'vBestUnchosen_neg'}))==0
                    eval(['ws.v=ws.dct(:, col.' request.allvars{v,1}  ');'])
                    ws.v=ws.v - min(ws.v); % Adjust floor
                    ws.v= ws.v./ max(ws.v);  % Scale by ceiling
                    eval(['ws.dct(:, col.' request.allvars{v,1}  ')=ws.v;'])
                    %             disp([request.allvars{v} '  - ct range: ' num2str(min(ws.v)) '    -    '  num2str(max(ws.v))])
                end
            end
        end
    end
    
    for r=1:length(log.rois)
        wr.cf=d_roimatrix_cF(:,:,1, r,s);
        wr.cf=flipud(wr.cf);  % Revert EnvThreat
        wr.ct=d_roimatrix_ct(:,:,1, r,s);
        wr.ct=flipud(wr.ct);  % Revert EnvThreat
        wr.dcf=ws.dcf;  wr.dct=ws.dct; 
        wr.dcf(:,col.b)=wr.cf(:); wr.dct(:,col.b)=wr.ct(:); 
        
        %
        d_betas{r,2}=[d_betas{r,2}; wr.dcf];
        d_betas{r,3}=[d_betas{r,3}; wr.dct];
        wr=[];
    end
    ws=[];
end

% (4) Check that lineup between variables and betas is correct
request.do_checklineup=0;   % ONLY works if variable re-scaling if off
if request.do_checklineup
    request.checksub=randi(log.n_subjs);     
    request.checkvaronce={'EnvThreat';'NTokens';}; request.checkvarbothtasks={'pLoss';'PavConflict';  'vChosen'; 'vBestUnchosen'; 'EVConflict'};
    request.checkvar =  [request.checkvaronce; request.checkvarbothtasks];  d_check =[log.rois cell(log.n_rois,2)];  % d_check{roi, task}{Variable}
    for t=1:2
        wt.subd= d_betas{r, t+1}(d_betas{r, t+1}(:, col.Subject)==request.checksub, :);
        for r=1:log.n_rois
            for e=1:6
                for n=1:6
                    wt.d= wt.subd(wt.subd(:, col.EnvThreatOrig).*6==e &  wt.subd(:, col.NTokensOrig)./2==n, :); 
                    for v=1:length(request.checkvar )
                        eval(['d_check{r,t+1}{col.' request.checkvar{v} '}(7-e,n) = unique(wt.d(:, col.' request.checkvar{v} '));'])
                        
                    end
                    d_check{r,t+1}{col.b}(7-e,n) = unique(wt.d(:, col.b)); 
                end
            end
        end
    end
    
    f.subplotcols=  length(request.checkvaronce)  + length(request.checkvarbothtasks)*2 +1 +2;  % Plot 
    f.subplot_VerHorz=[0.01 0.005];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.05 0.01];  f.figwidth= 1400; f.figheight=log.n_rois*150;    f.nancol=[0 0 0]; k=1;  fm=figure('Position',[100,70,f.figwidth,f.figheight], 'Color',[1 1 1]);
    for r=1:log.n_rois
        subtightplot(log.n_rois, f.subplotcols,  k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1; 
        text(-0.5,0.5, log.rois{r},'FontSize', 15); axis off
        
        for v=1:length(request.checkvaronce)        
            subtightplot(log.n_rois, f.subplotcols,  k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
            eval([' imagescnan(d_check{r,1+1}{col.' request.checkvaronce{v} '}, ''nancolor'', f.nancol)']), axis off, axis square, colorbar
            if r==1; title(request.checkvaronce{v}), end
        end
        k=k+1; 
        
        for t=1:2
            for v=1:length(request.checkvarbothtasks)
                subtightplot(log.n_rois, f.subplotcols,  k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
                eval([' imagescnan(d_check{r,t+1}{col.' request.checkvarbothtasks{v} '}, ''nancolor'', f.nancol)']), axis off, axis square, colorbar
                if r==1; title(['[' request.task{t} ']  '  request.checkvarbothtasks{v}]), end
            end
            if t==1, k=k+1; end
        end
        
    end
end

%%


% STEPWISE GLM (Fixed only)
% r_swglm{t}: (1) name, (2) ivs + estimates, (3) mdl, (4) simulated b, (5)  t stat (6) sig t stat (7) Fixed-fx plot
if request.do_glm_stepwise
    
    close all hidden
    request.do_stepwiseglm_plot=1;
    r_swglm = repmat({[log.rois cell(log.n_rois,6)]},1,2); 
    request.do_stepwiseglm_cols={'roi';'IvsEsts';'Mdl';'sim b (mixed)';'sim b tstat';'sim b tstat-sig'; 'sim b (fixed)'};
    allivs={'EnvThreat'; 'NTokens'; 'pLoss';  'Entropy'; 'EV'; 
%         'PavConflict';  
%         'EVConflict';
    }; request.psig=0.05;
    col.allivcols=[];  for i=1:length(allivs), eval(['col.allivcols=[col.allivcols col.' allivs{i} '];']), end
    w.DummySub= sortrows(repmat([eye(log.n_subjs) (1:log.n_subjs)'],36,1), log.n_subjs+1); w.DummySub     = w.DummySub(:, 1:end-1); 
    for r=1:log.n_rois % Run 
        for t=1:2
            wr.x=[d_betas{r, t+1}(:, col.allivcols)  w.DummySub];
            
             mdl=stepwiseglm(wr.x, d_betas{r, t+1}(:, col.b), 'linear', 'Upper','linear', 'Criterion', 'bic',  'CategoricalVars',  length(allivs)+1:  length(allivs)+log.n_subjs, 'PredictorVars', 1:length(allivs)) ;
            
%             mdl = GeneralizedLinearModel.stepwise(wr.x, d_betas{r, t+1}(:, col.b), 'constant', 'Criterion', 'bic',  'CategoricalVars',  length(allivs)+1:  length(allivs)+log.n_subjs, 'PredictorVars', 1:length(allivs)) ;
            
            
            
            r_swglm{t}{r,2}(:,  1)= mdl.CoefficientNames';
%             r_swglm{t}{r,2}(:,2)=  num2cell(double( mdl.Coefficients(:,1) ));   % Mac
            r_swglm{t}{r,2}(:,2)=  num2cell(table2array( mdl.Coefficients(:,1) ));   % PC 
            r_swglm{t}{r,3}=mdl; 
            
            % Name predictors
            for p=2:size(r_swglm{t}{r,2},1)
                if isempty(strfind(r_swglm{t}{r,2}{p,1}, ':'))
                    r_swglm{t}{r,2}{p,1} = allivs{str2num(r_swglm{t}{r,2}{p,1}(2:end))};
                else % Interaction terms 
                    wr.i1=  r_swglm{t}{r,2}{p,1}(1: strfind(r_swglm{t}{r,2}{p,1}, ':')-1); 
                    wr.i2= r_swglm{t}{r,2}{p,1}( strfind(r_swglm{t}{r,2}{p,1}, ':')+1:end); 
                    r_swglm{t}{r,2}{p,1} = [ allivs{str2num(wr.i1(2:end))} '*' allivs{str2num(wr.i2(2:end))}];
                end
            end
            wr=[]; 
        end
    end
    for r=1:log.n_rois % Compile for plot
        for t=1:2
            
            % Simulations + sort by subject + plot like it's mixed 
            wr.b=predict(r_swglm{t}{r,3}, [d_betas{r, t+1}(:, col.allivcols)  w.DummySub]);  % Simulated betas
            wr.b=[wr.b sortrows(repmat(1:log.n_subjs, 1,6*6)')]; 
            wr.bs=nan(log.n_subjs,6*6);
            for s=1:log.n_subjs
                wr.bs(s,:)= wr.b( wr.b(:, end)==s, 1); 
            end
            r_swglm{t}{r,4}  = flipud(reshape(mean(wr.bs),6,6)); 
            [wr.h wr.p wr.ci wr.stats]= ttest(wr.bs);
            r_swglm{t}{r,5}  = flipud(reshape(wr.stats.tstat,6,6)); 
            wr.stats.tstat(wr.p>request.psig)=nan; 
            r_swglm{t}{r,6}  = flipud(reshape(wr.stats.tstat,6,6)); 
            
            % Fixed plot?? (#7)
            wr.bm =reshape(wr.b(:,1), 6,6, log.n_subjs);
            r_swglm{t}{r,7} = flipud(mean(wr.bm,3));
            wr=[]; 
        end
    end
    
    % Table output
    rt_swglm{1}=  [[ {'Exp' }; log.rois] [[{'(Intercept)'} allivs']; repmat({' '}, log.n_rois, length(allivs)+1)]]; rt_swglm{2}=rt_swglm{1}; rt_swglm{2}{1,1}='Ctrl'; 
    for r=1:log.n_rois
        for t=1:2
            for p=1:size(r_swglm{t}{r,2},1)
                if sum(strcmp(rt_swglm{t}(1,:), r_swglm{t}{r,2}{p,1}))==1
                    rt_swglm{t}{r+1,  find(strcmp(rt_swglm{t}(1,:), r_swglm{t}{r,2}{p,1}))} = num2str(r_swglm{t}{r,2}{p,2}, 3);
                else % non-standard IVs (interaction terms)
                    rt_swglm{t} = [rt_swglm{t} [r_swglm{t}{r,2}{p,1}; repmat({' '}, log.n_rois,1)]];
                    rt_swglm{t}{r+1, end} = num2str(r_swglm{t}{r,2}{p,2},3);
                end
            end
        end
    end
    openvar a, a=repmat({' '}, log.n_rois*2+3, 20);  a(1:log.n_rois+1, 1:size(rt_swglm{1},2))= rt_swglm{1};  a(log.n_rois+3:log.n_rois*2+3, 1:size(rt_swglm{2},2))= rt_swglm{2}; 
    
       
    % PLOT
    if request.do_stepwiseglm_plot
        request.do_stepwiseglm_plotwhich=[4 5]; % see details on r_swglm for details of identity
%         f.subplotcols=4;  f.figwidth= 600; f.figheight=500; f.subplot_VerHorz=[0.1 0.05];  f.fig_BotTop=[0.05 0.05];  f.fig_LeftRight=[0.02 0.05];  f.nancol=[0 0 0];
        f.subplotcols=4;  f.figwidth= 600; f.figheight=600; f.subplot_VerHorz=[0.03 0.03];  f.fig_BotTop=[0.05 0.05];  f.fig_LeftRight=[0.02 0.05];  f.nancol=[0 0 0];
        for w=1:length(request.do_stepwiseglm_plotwhich)
            figure('Name', ['GLM-sim betas ('  request.do_stepwiseglm_cols{request.do_stepwiseglm_plotwhich(w)}  ')'], 'color', 'w', 'Position',[50+w*50,100,f.figwidth,f.figheight]); f.FontSize=15; k=1;
            for r=1:log.n_rois
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                text(0.2,0.5, log.rois{r}, 'FontSize', f.FontSize);
                if r==1; title(request.do_stepwiseglm_cols{request.do_stepwiseglm_plotwhich(w)} , 'FontSize', f.FontSize);end
                axis off, k=k+1;
                
                % Plot cF
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( r_swglm{1}{r,  request.do_stepwiseglm_plotwhich(w)}, 'NaNColor', f.nancol);  
                set(gca, 'XTick', [])
                try colorbar; end
                if r==1; title('sim cF', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;
                
                % Plot ct
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( r_swglm{2}{r,  request.do_stepwiseglm_plotwhich(w)}, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
                if r==1; title('sim ct', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;
                
                % Plot cF - ct
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( r_swglm{1}{r,  request.do_stepwiseglm_plotwhich(w)} - r_swglm{2}{r,  request.do_stepwiseglm_plotwhich(w)} , 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
                if r==1; title('sim cF - ct', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;
            end
        end
    end
    
end


%% (1) GLM on betas
%  d_roimatrix= EnvThreat x NTokens x choice x roi x subject
%                       EnvThreat in reverse order (i.e. as per imagesc plot)


if request.do_glm
    
    
close all hidden
    request.do_glmplot=1;
    request.do_simfromglm=1;
    % WHICH IVS?
    request.psig=0.05;
    ivs={
        'EnvThreat'; 'NTokens';
        'pLoss'; 
        'Entropy'; 'EV';
        
%         'vChosen';'vBestUnchosen';
%         'vBestUnchosen_neg';
%         'PavConflict'
%         'EVConflict'
        };
%     ivs={'Entropy';'PavConflict';};
%     ivs={'Entropy';'pLoss';};
%     ivs={'Entropy';'EVConflict';};
%     ivs={'Entropy';'vBestUnchosen_neg'};
%     ivs={'EntropyNTok';'vBestUnchosen_neg'};

    ivs={'vChosen';'vBestUnchosen'};
    
    % MANUAL alterantions
    disp('Manual alterations to task space !!'); 
    for r=1:log.n_rois
        
        d_betas{r,3}(:, col.PavConflict)=0;
        d_betas{r,3}(:, col.EVConflict)=0;
        
        d_betas{r,2}(:, col.pLossAversive)=d_betas{r,2}(:, col.pLoss); 
        d_betas{r,3}(:, col.pLossAversive)=0; 
    end
    
    % Fixed GLM
    col.ivcols=[];  for i=1:length(ivs), eval(['col.ivcols=[col.ivcols col.' ivs{i} '];']), end
    col.ivcols_fixed = [col.ivcols col.Subject];
    r_fglm{1}= [[{' '};  log.rois]  [[{'Cst'} ivs']; cell(log.n_rois, length(ivs)+1)]]; r_fglm{2}= r_fglm{1};
    for r=1:length(log.rois)
        disp([log.rois{r} '     ----------------------------------'])
        
        %     % Residuals, both tasks
        %     disp('Residuals, both tasks:');
        %     wt.dboth= [d_betas{r,2} ;  d_betas{r,3}];
        %     wt.dboth(:, col.pLossAversive)=wt.dboth(:, col.pLoss);
        %     wt.dboth(size(d_betas{r,2},1)+1:end , col.pLossAversive) =0;
        %     [wt.b1 wt.bint wt.rr1 wt.rint wt.stats1]=regress(wt.dboth(:, col.b),  [ones(size(wt.dboth,1),1)   wt.dboth(:, col.Entropy)]); wt.stats1
        %     [wt.b2 wt.bint wt.rr2 wt.rint wt.stats2]=regress(wt.rr1,  [ones(size(wt.dboth,1),1)   wt.dboth(:, col.ivcols_fixed)]); wt.stats2
        
        for t=1:2
            disp(['Cond ' request.task{t} ':'])
            
            %         % Residuals, single task
            %         [wt.b1 wt.bint wt.rr wt.rint wt.stats1]=regress(d_betas{r,t+1}(:,col.b),  [ones(size(d_betas{r,t+1},1),1)   d_betas{r,t+1}(:, col.Entropy)]); wt.stats1
            %         [wt.b2 wt.bint wt.rr2 wt.rint wt.stats2]=regress(wt.rr,  [ones(size(d_betas{r,t+1},1),1)   d_betas{r,t+1}(:, col.ivcols_fixed)]); wt.stats2
            
            
            % Predict what?
            wt.y= d_betas{r,t+1}(:,col.b);  % glm on Betas
            
            % Execute GLM #######
            [wt.beta, wt.d, wt.s]  =glmfit(d_betas{r,t+1}(:,col.ivcols_fixed), wt.y);
            wt.beta=num2cell(wt.beta); wt.beta(wt.s.p>request.psig)= repmat({' '}, sum(wt.s.p>request.psig),1);
            wt.res=wt.beta;
            %         wt.res=num2cell(wt.s.p);
            
            % Record
            r_fglm{t}(r+1,2:end) = wt.res(1:end-1)';
            wt=[];
        end
        
        
%         BOTH 
%         
%         wr.dboth= [d_betas{r,2} ;  d_betas{r,3}];
%         wr.y=wr.dboth(:, col.b);
%         [wr.beta, wr.d, wr.s]  =glmfit(wr.dboth(:,col.ivcols_fixed), wr.y);
%         disp(wr.beta(2:end-1)')
%         disp(wr.s.p(2:end-1)')
    end
%     openvar r_fglm{1}, openvar r_fglm{2}
    
    % First + Second level GLM 
    d_betab= [log.rois  repmat({nan(log.n_subjs,length(ivs))}, log.n_rois,2) ];
    d_betabcon= [log.rois  repmat({nan(log.n_subjs,1)}, log.n_rois,2) ]; % constant
    % r_glm{1}= [[{' '};  log.rois]  [[{'Cst'} ivs']; cell(log.n_rois, length(ivs)+1)]]; r_glm{2}= r_glm{1}; r_glmp=r_glm;
    r_glm{1}= [[{' '};  log.rois]  [  ivs'; cell(log.n_rois, length(ivs))]  ]; r_glm{2}= r_glm{1};  r_glm{3}= r_glm{1}; r_glmp=r_glm;
    for r=1:length(log.rois)
        for t=1:2
            for s=1:log.n_subjs
                ws.d=d_betas{r,t+1}(d_betas{r,t+1}(:, col.Subject)==s, :);
                [ws.b ws.dd ws.stats]= glmfit(ws.d(:, col.ivcols), ws.d(:,col.b));
                
%                 
%                 % Stopped here - get R2
%                 ws.predb= glmval(ws.b, ws.d(:, col.ivcols),'identity'); 
                d_betab{r,t+1}(s,:)=  ws.b(2:end)' ;
%                 
%                 
%                 ws.stats
%                 error
                d_betabcon{r,t+1}(s)=ws.b(1);
                
                
                wt=[];
            end
            
            % Second level: one-sample ttests per iv
            [wt.h wt.p wt.ci wt.st]=ttest(d_betab{r,t+1});
            wt.tstat = num2cell(wt.st.tstat); wt.tstat( wt.p>request.psig) = repmat({' '}, 1, sum(wt.p>request.psig));
            r_glm{t}(r+1, 2:end) =  wt.tstat;
            r_glmp{t}(r+1, 2:end) =  num2cell(wt.p);
            
            wt=[];
        end
    end
    
    % Plot Mixed GLM
    request.plotrange=[-10 25];
%     request.plotrange=[];
    f.subplotcols=4;  f.figwidth= 1000;f.figheight=800;f.subplot_VerHorz=[0.05 0.05];  f.fig_BotTop=[0.05 0.05];  f.fig_LeftRight=[0.01 0.2];  f.nancol=[0 0 0];
    if request.do_glmplot, figure('color', 'w', 'Position',[200,70,f.figwidth,f.figheight]); end,  f.FontSize=15; k=1;
    request.ivnames = cellfun(@(x)  request.allvars{ strcmp(request.allvars(:,1), x),2}, ivs, 'UniformOutput',0);
    for r=1:length(log.rois)
        % Second level cF > ct
        [wt.h wt.p wt.ci wt.st]=ttest( d_betab{r,1+1} -  d_betab{r,1+2 });
        wt.tstat = num2cell(wt.st.tstat); wt.tstat( wt.p>request.psig) = repmat({' '}, 1, sum(wt.p>request.psig));
        r_glm{3}(r+1, 2:end) =  wt.tstat;
        r_glmp{3}(r+1, 2:end) =  num2cell(wt.p);
        
        % PLOT
        if request.do_glmplot
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;  % Name of ROI
            text(0.5, 0.8,  log.rois{r},'FontSize', f.FontSize), axis off;
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;    % cF
            barwitherr(std(d_betab{r,1+1})./sqrt(log.n_subjs), mean(d_betab{r,1+1}), 'y');
            ylim(request.plotrange)
            %     barwitherr(std(d_betab{r,1+1}), mean(d_betab{r,1+1}), 'y')
            
            if r==1; title(['[cF] ' ] ,'FontSize', f.FontSize); end, xlim([0 length(ivs)+1]);
            set(gca, 'xticklabel',  request.ivnames) %  set(gca, 'xticklabel', [{'Cst'}; ivs])
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;    % ct
            barwitherr( std(d_betab{r,2+1})./sqrt(log.n_subjs), mean(d_betab{r,2+1}), 'y');
            ylim(request.plotrange)
            %     barwitherr( std(d_betab{r,2+1}), mean(d_betab{r,2+1}), 'y')
            
            if r==1; title(['[ct] ' ] ,'FontSize', f.FontSize), end; xlim([0 length(ivs)+1]);
            set(gca, 'xticklabel',  request.ivnames) %  set(gca, 'xticklabel', [{'Cst'}; ivs])
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;   % cF > ct
            barwitherr( std( d_betab{r,1+1} -  d_betab{r,2+1} )./sqrt(log.n_subjs), mean(d_betab{r,1+1} -  d_betab{r,2+1}), 'y');
            ylim(request.plotrange)
            %         barwitherr( std( d_betab{r,1+1} -  d_betab{r,2+1} ), mean(d_betab{r,1+1} -  d_betab{r,2+1}), 'y')
            if r==1; title(['[cF - ct] ' ] ,'FontSize', f.FontSize), end, xlim([0 length(ivs)+1]);
            set(gca, 'xticklabel',  request.ivnames); %  set(gca, 'xticklabel', [{'Cst'}; ivs])
        end
    end
    openvar r_glm{1}, openvar r_glm{2}, openvar r_glm{3}
    
    % Simulate betas from GLM 
    if request.do_simfromglm 
        r_simt = [log.rois repmat({nan(6,6)},log.n_rois,3)];   r_simm = [log.rois repmat({zeros(6,6)},log.n_rois,3)]; r_simts= r_simt ;
        d_simglm= [log.rois  repmat({nan(log.n_subjs, 6*6)}, log.n_rois,2)];
    
        for r=1:log.n_rois
            for s=1:log.n_subjs   % [1st + 2nd level GLM] Compile simulated data per subject
                for t=1:2
                    % Nans are not a problem - they are zero'd out
                    d_simglm{r,t+1}(s,:) = (repmat(d_betabcon{r,t+1}(s) , 6*6,1) +   sum(d_betas{r, t+1}(d_betas{r, t+1}(:, col.Subject) ==s, col.ivcols) .*repmat(d_betab{r,t+1}(s,:) , 6*6,1)  ,2))';
                    % w.des(:, [col.EnvThreatOrig col.NTokens col.ivcols]) = d_betas{r,t+1}(d_betas{r, t+1}(:, col.Subject) ==s, [col.EnvThreatOrig col.NTokens col.ivcols]);  % Capture design in the same dimensions
                end
            end
            
            % mean
            r_simm{r,1+1} = flipud(reshape(mean(d_simglm{r,1+1}),6,6));
            r_simm{r,2+1} = flipud(reshape(mean(d_simglm{r,2+1}),6,6));
            r_simm{r,3+1} = flipud(reshape(mean(  d_simglm{r,1+1} - d_simglm{r,1+2} ),6,6));
            % subplot(1,2,1),  imagesc( flipud(reshape( w.des(:, col.EnvThreatOrig),6,6))), axis square
            % subplot(1,2,2),  imagesc( flipud(reshape( w.des(:, col.NTokens),6,6))), axis square
            
            % t stat + t stat sig
            [wr.h wr.p wr.ci wr.stat]= ttest(d_simglm{r,1+1});  % cF
            r_simt{r,1+1} = flipud(reshape(wr.stat.tstat,6,6));
            wr.stat.tstat(wr.p>request.psig)= nan; 
            r_simts{r,1+1} = flipud(reshape(wr.stat.tstat,6,6)); 
            [wr.h wr.p wr.ci wr.stat]= ttest(d_simglm{r,1+2});  % ct 
            r_simt{r,1+2} = flipud(reshape(wr.stat.tstat,6,6)); 
            wr.stat.tstat(wr.p>request.psig)= nan; 
            r_simts{r,1+2} = flipud(reshape(wr.stat.tstat,6,6)); 
            [wr.h wr.p wr.ci wr.stat]= ttest( d_simglm{r,1+1} -  d_simglm{r,1+2});  % cF > ct
            r_simt{r,1+3} = flipud(reshape(wr.stat.tstat,6,6));
            wr.stat.tstat(wr.p>request.psig)= nan; 
            r_simts{r,1+3} = flipud(reshape(wr.stat.tstat,6,6));
            
        end
        
        % PLOT
        request.do_glmplotwhat={'r_simt'; 'r_simts'}; % r_simm= mean, r_simt= tstat, r_simts= tstat sig
        f.subplotcols=4; f.figwidth= 600;f.figheight=500; f.subplot_VerHorz=[0.03 0.05];  
        f.fig_BotTop=[0.02 0.05];  f.fig_LeftRight=[0.02 0.02];  f.nancol=[0 0 0];
        for w=1:length(request.do_glmplotwhat)
            figure('Name', ['GLM-sim betas (' request.do_glmplotwhat{w} ')'], 'color', 'w', 'Position',[50+w*50,100,f.figwidth,f.figheight]); f.FontSize=15; k=1;
            eval(['d_plot=' request.do_glmplotwhat{w} ';'])
            for r=1:log.n_rois
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                text(0.2,0.5, log.rois{r}, 'FontSize', f.FontSize);
                axis off, k=k+1;

                % Plot cF
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( d_plot{r,1+1}, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
                if r==1; title('sim cF', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;

                % Plot ct
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( d_plot{r,1+2}, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
                if r==1; title('sim ct', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;

                % Plot cF - ct
                subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagescnan( d_plot{r,1+3}, 'NaNColor', f.nancol);  colorbar; set(gca, 'XTick', [])
                if r==1; title('sim cF - ct', 'FontSize', f.FontSize);end
                axis square, axis off,k=k+1;
            end
        end
    end
    disp('##########################################'), disp('GLM plot x names:'), disp( [{'Cst'}; ivs])

    % Behaviour correlations
    if request.beh_correlations
        request.beh4corr={'st.State';'st.Trait'; 'per.cF_Reject';'m.cF_j'; 'm.j_cFmct'; 'infoexnull.cF'}; 
        wb.d= cellfun(@(x) cell2mat(d_behscores(2:end, find(strcmp(d_behscores(1,:), x) ) )), request.beh4corr, 'UniformOutput',0);
        d_behcor= [request.beh4corr wb.d];  % Col 1=behaviour name, Col 2= beh scores, Col 3=results
        request.betab_corr=ivs;
        
        % MANUALLY SPECIFY
        r_behcorr = [cell(2, 1) [request.beh4corr';  cell(1, length(request.beh4corr))]]; k=2;
        for r=1:log.n_rois
            r_behcorr{k,1}= log.rois{r};  k=k+1;
            for t=1:2
                for i=1:length(request.betab_corr)
                    r_behcorr {k,1}= [request.task{t} ' ' request.betab_corr{i}];
                    for b=1:size(d_behcor,1)
                        wb.beh= d_behcor{b,2};
                        [wi.r wi.p]= corr( d_betab{r,t+1}(:, strcmp( ivs,  request.betab_corr{i})), wb.beh);
                        if wi.p<0.001
                            r_behcorr {k, b+1}=  ['*** r=' num2str(wi.r) ', p=' num2str(wi.p)] ;
                        elseif wi.p<0.01
                            r_behcorr {k, b+1}=  ['** r=' num2str(wi.r) ', p=' num2str(wi.p)] ;
                        elseif wi.p<request.psig
                            r_behcorr {k, b+1}=  ['* r=' num2str(wi.r) ', p=' num2str(wi.p)] ;
                        else  r_behcorr {k, b+1}=' ';
                        end
                        wi=[];  
                    end, k=k+1;
                end
            end
        end
        openvar r_behcorr
    end
    
    
    
end


%% (3) On-sample stats on cell betas + Plot

if request.plot_tstat
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
%     d_roinames=request.ROIs;
    d_roinames=log.rois;
    %
    k=1; f.FontSize=25;  f.subplotcols=4;  f.subplot_VerHorz=[0.005 0.03]; f.fig_BotTop=[0.001 0.035]; f.fig_LeftRight=[0.1 0.1];  f.nancol=[0 0 0];
    figure('color','w', 'Position', [0 50 600 1200])
    for r=1:log.n_rois
        wr.d=  [d_plot{2}(:,:,r) d_plot{1}(:,:,r)];
        wr.range=[min(wr.d(:)) max(wr.d(:))];
        
        subtightplot(log.n_rois , f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        %     text(0.2,0.5, log.rois{r}, 'FontSize', f.FontSize);
        text(0.5,0.5, d_roinames{r}, 'FontSize', f.FontSize);
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

% error('DONE! :)')

for o=1:1 % CURRENTLY UNUSED :) 
%% (3) PCA

if request.do_pca
    request.psig=0.05;
    ivs={
        'EnvThreat'; 'NTokens'; 'pLoss'; 'Entropy'; 
%         'EV';
%         'PavConflict'; 'EVConflict'
        };
    col.ivcols=[];  for i=1:length(ivs), eval(['col.ivcols=[col.ivcols col.' ivs{i} '];']), end,
    request.ivnames = cellfun(@(x)  request.allvars{ strcmp(request.allvars(:,1), x),2}, ivs, 'UniformOutput',0);
    
    % PCA + plot
    request.pca.n_components=[4 4]; % MANUALLY specify from plots. # cF, # ct
    figure('color','w', 'position',[100 400 1200 600])
    for t=1:2 % Identify/set up components  
        % Identify components (in IVs)
        %     Looking at COEFS tells you how much each of your original raw IVs loads onto each of the components. from the relative loading etc, you then need
        %     to decide what the psychological character of the component is.
        %     Components are order according to increasing amount of variance accounted for individually (1st component = most variance)
        wt.x = d_betas{r,t+1}(:, col.ivcols);
        [wt.coeff, wt.scores, wt.latent] = princomp(wt.x ); % find what the principal components are - a mixture of your original raw IVs
        subplot(2,2, t)
        %     subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
        plot(cumsum(wt.latent)/sum(wt.latent)) % how much variance explained by n components (no particular order) - allows you to see how many components you are looking for
        xlim([1  length(ivs)]), set(gca,'xtick', 1:length(ivs))
        xlabel('No. components'), ylabel('Cumulative variance explained')
        title([ '[' request.task{t} '] Variance by # components'])
        
        % IV loadings
        subplot(2,2, 2 +t)
        %     subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
        bar(wt.coeff(:,1:request.pca.n_components(t)))
        set(gca,'xticklabel', request.ivnames), xlim([0  length(ivs)+1])
        ylabel('Loading coefficient'); title([ '[' request.task{t} '] IV loading per component (assumed ' num2str(request.pca.n_components(t)) ' comps)']);
        %     if t==1, k=k+1; end
        
        % Store
        d_pcacomp{1, t}= wt.coeff;
        d_pcacomp{2, t}= wt.scores;
        d_pcacomp{3, t}= wt.latent;
        
        wt=[];
    end
    f.subplotcols=3;  f.subplot_VerHorz=[0.1 0.15];  f.fig_BotTop=[0.05 0.05];  f.fig_LeftRight=[0.05 0.05];  f.figwidth= 800;   f.figheight=1000; f.nancol=[0 0 0];
    figure('color', 'w', 'Position',[100,100,f.figwidth,f.figheight]); f.FontSize=15; k=1;
    for r=1:log.n_rois
        subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;  % Name of ROI
        text(0.5, 0.8,  log.rois{r},'FontSize', f.FontSize), axis off;
        for t=1:2
            wr.y = d_betas{r,t+1}(:, col.b);
            wr.comp_ivs= d_pcacomp{2, t}(:,1:request.pca.n_components(t)); % SCORES(1:N_COMPS) =  IVs, of each of the different principal components. See cumsum figure above to determine how many components to include
            
            % Redo GLM
            [wr.b, wr.dev, wr.stats] = glmfit(wr.comp_ivs, wr.y);
            %     wr.b = num2cell(wr.b); %     repmat({0}, sum(wr.stats.p > request.psig),1);
            wr.b ( wr.stats.p > request.psig)=0;
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;
            bar(wr.b), xlim([0 request.pca.n_components(t)+1]); set(gca,'xtick', 1:request.pca.n_components(t))
            xlabel('Component no.'), ylabel('Effect of comp on fMRI betas (sig only)')
            if r==1; title([ '[' request.task{t} '] Component betas']), end
        end
    end
end


%%



%% (4) Apply contrast weights
%	roimatrix: EnvThreat x NTokens x Choice x ROI

if request.var_conweights
    request.var_conweights_plot=1;
    request.psig=0.05;
    ivs={
        'EnvThreat'; 'NTokens';
        'pLoss';'Entropy'; 'EntropyNTok'; 
        'EV';
        'vChosen';'vBestUnchosen';
        'vBestUnchosen_neg';
        'PavConflict'
        'EVConflict'
        };
%     ivs={'Entropy';'pLoss';};
%     ivs={'Entropy'};
%     ivs={'EntropyNTok';'PavConflict'};
%     ivs={'EntropyNTok';'vBestUnchosen_neg'};
        
    
    
    % Apply contrast weights + ttest
    d_weightb= [log.rois repmat({nan(log.n_subjs,length(ivs))}, log.n_rois,2)];
    r_weightb{1}=  [[{' ' } ivs']; [log.rois cell(log.n_rois,length(ivs))]];  r_weightb{2}= r_weightb{1}; r_weightb{3}=r_weightb{1};
    for r=1:log.n_rois  
        
        % Apply weights to get subject-level 
        for s=1:log.n_subjs
            for t=1:2
                for v=1:length(ivs)
                    eval(['d_weightb{r,t+1}(s,v)= sum(d_betas{r,t+1}(d_betas{r,t+1}(:,col.Subject)==s, col.' ivs{v} ').*d_betas{r,t+1}(d_betas{r,t+1}(:,col.Subject)==s, col.b));'])
                end
            end
        end
        
        % t-test
        for t=1:2
            [wr.h wr.p wr.ci wr.stats]=ttest(d_weightb{r,t+1});
            wr.stats.tstat=num2cell(wr.stats.tstat);
            wr.stats.tstat(wr.p>request.psig)=  repmat({' '}, 1, sum(wr.p>request.psig));
            wr.res=wr.stats.tstat;
%             wr.res=num2cell(wr.p);

            % Output
            r_weightb{t}(r+1, 2:end) = wr.res;
        end
        
        % Comparing cF and ct 
        [wr.h wr.p wr.ci wr.stats]=ttest( d_weightb{r,1+1} - d_weightb{r,2+1} );
        wr.stats.tstat=num2cell(wr.stats.tstat);
        wr.stats.tstat(wr.p>request.psig)=  repmat({' '}, 1, sum(wr.p>request.psig));
        wr.res=wr.stats.tstat;
        %             wr.res=num2cell(wr.p);
        
        % Output
        r_weightb{3}(r+1, 2:end) = wr.res;
        
        wr=[];
    end
    openvar r_weightb{1}, openvar r_weightb{2}, openvar r_weightb{3} 
    
    % Plot
    f.subplotcols=4;  f.subplot_VerHorz=[0.1 0.05];  f.fig_BotTop=[0.05 0.05];  f.fig_LeftRight=[0.02 0.05];  f.figwidth= 800;   f.figheight=1000; f.nancol=[0 0 0];
    figure('Name', 'Variable tracking (via beta weighting)', 'color', 'w', 'Position',[200,100,f.figwidth,f.figheight]); f.FontSize=15; k=1;
    request.ivnames = cellfun(@(x)  request.allvars{ strcmp(request.allvars(:,1), x),2}, ivs, 'UniformOutput',0);
    if request.var_conweights_plot
        for r=1:log.n_rois
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;  % Name of ROI
            text(0.5, 0.8,  log.rois{r},'FontSize', f.FontSize), axis off;
            
            % cF 
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;    % cF
            barwitherr(std(d_weightb{r,1+1})./sqrt(log.n_subjs),  mean(d_weightb{r,1+1}), 'y')            
%             barwitherr(std(d_weightb{r,1+1}),  mean(d_weightb{r,1+1}), 'y')
            if r==1; title(['[cF] ' ] ,'FontSize', f.FontSize); end, xlim([0.5 length(ivs)+0.5])
            set(gca, 'xticklabel',  request.ivnames) %  set(gca, 'xticklabel', [{'Cst'}; ivs])
            
            % ct
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;    % ct
            barwitherr(std(d_weightb{r,2+1})./sqrt(log.n_subjs),  mean(d_weightb{r,2+1}), 'y')            
%             barwitherr(std(d_weightb{r,1+1}),  mean(d_weightb{r,1+1}), 'y')
            if r==1; title(['[ct] ' ] ,'FontSize', f.FontSize); end, xlim([0.5 length(ivs)+0.5])
            set(gca, 'xticklabel',  request.ivnames) %  set(gca, 'xticklabel', [{'Cst'}; ivs])
            
            
            % cF > ct
            subtightplot(log.n_rois, f.subplotcols, k,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight); k=k+1;   % cF > ct
            barwitherr(std(d_weightb{r,1+1} - d_weightb{r,2+1} )./sqrt(log.n_subjs),  mean(d_weightb{r,1+1} - d_weightb{r,2+1} ), 'y')            
            if r==1; title(['[cF > ct] ' ] ,'FontSize', f.FontSize); end, xlim([0.5 length(ivs)+0.5])
            set(gca, 'xticklabel',  request.ivnames) %  set(gca, 'xticklabel', [{'Cst'}; ivs])
        end
    end
    
end




%% (2) Plot mean betas  
%	roimatrix: EnvThreat x NTokens x Choice x ROI
%   n_roimatrix: e x n x choice x task

request.plotacrossdesign=0;

if request.plotacrossdesign==1
    m_roimatrix_cF=nan(log.n_levels,log.n_levels,log.n_choices, log.n_rois); m_roimatrix_ct=m_roimatrix_cF;  m_roimatrix_comb=m_roimatrix_cF;  n_roimatrix=nan(log.n_levels,log.n_levels,log.n_choices,2); % e x n x choice x task
    for e=1:log.n_levels  % Compile
        for n=1:log.n_levels
            for c=1:log.n_choices
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
    
    % Plot settings
    log.betarange=[]; %[-1 1]; %[]; %[0 10];
    if strcmp(request.FLmodel(1), 'f')==1 &&  isempty(log.betarange)==1   % Beta range? (standardize within ROIs; overridden by requests)
        log.meanbrange=cell(log.n_rois, 2);   for r=1:log.n_rois % tweak here to calculate ranges for cF and ct separately
            %             log.meanbrange{r,1}=  [min( min(min(m_roimatrix_cF(:,:, :,r))))    max( max(max(m_roimatrix_cF(:,:, :,r)))) ];
            %             log.meanbrange{r,2}=  [min( min(min(m_roimatrix_ct(:,:, :,r))))     max( max(max(m_roimatrix_ct(:,:, :,r)))) ];
            log.meanbrange{r,1}=   [min( min(min([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))    max( max(max([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))  ];
            log.meanbrange{r,2}=log.meanbrange{r,1};
            disp('Natural range for betas (for means plot): '); disp(['   '  log.rois{r} ':    '   num2str(log.meanbrange{r,1}(1)) ' to ' num2str(log.meanbrange{r,1}(2))])
        end
        log.betarange =log.meanbrange;
    elseif strcmp(request.FLmodel(1), 't')==1 &&  isempty(log.betarange)==1
        log.meanbrange=cell(log.n_rois, 2); for r=1:log.n_rois
            log.meanbrange{r,1}=   [min( min(min([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))    max( max(max([m_roimatrix_cF(:,:, :,r) m_roimatrix_ct(:,:, :,r)])))  ];
            log.meanbrange{r,2}=log.meanbrange{r,1};
            disp('Natural range for betas (for means plot): '); disp(['   '  log.rois{r} ':    '   num2str(log.meanbrange{r,1}(1)) ' to ' num2str(log.meanbrange{r,1}(2))])
        end ; log.betarange =log.meanbrange;
    elseif strcmp(request.FLmodel(1), 'f')==1;  log.betarange = repmat({log.betarange}, log.n_rois, 2);
    end
    
    % Plot
    f.figwidth= 1400; f.figheight=log.n_rois*250+150;   f.figheight=1000; f.nancol=[0 0 0]; if log.n_choices==3; f.subplotcols=8; f.plotspace=3;  f.plotrowsadd=2; f.subplot_VerHorz=[0.01 0.005];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.15 0.01];  else f.subplotcols=4; f.plotspace=2; f.plotrowsadd=0;   f.subplot_VerHorz=[0.01 0.005]; f.fig_BotTop=[0.001 0.025]; f.fig_LeftRight=[0.3 0.3];   end
    fm=figure('Name', [request.FLmodel ': Mean betas for each ROI [Accept, Reject, Explore, NoBomb, Bomb, Explore]   (min. no. subjects per cell = ' num2str(request.minsubs_cell) ')'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    for r=1:log.n_rois % Betas for Roi x task
        % cF
        subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        text(-0.5,0.5, log.rois{r},'FontSize', 15); axis off
        for c=1:log.n_choices
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_cF(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.n_choices==1; title('cF');  elseif r==1; title(['cF ' log.choicecF{c}(1:end-1)]); end
            if isempty(log.betarange)==0; caxis(log.betarange{r,1}); end
            %             colormap('summer')
        end
        
        % ct
        %         if log.n_choices==3
        %             subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+f.plotspace,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        %             text(0.3,0.5, log.rois{r},'FontSize', 8); axis off
        %         end
        for c=1:log.n_choices
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (r-1)*f.subplotcols+1+f.plotspace+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_ct(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.n_choices==1; title('ct');  elseif r==1; title(['ct ' log.choicect{c}(1:end-1)]); end
            if isempty(log.betarange)==0; caxis(log.betarange{r,2}); end
            %             colormap('summer')
        end
    end
    for o1=1:1 % No. of subjects (add to figure)
        log.range_nsubs=[8 log.n_subjs];
        log.range_ntrials=[8 14];
        if log.n_choices==3;
            t=1; c=1;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'No. of subjects included' ,'FontSize', 8); axis off
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagesc(n_roimatrix(:,:,c,t), log.range_nsubs);  axis off; axis square; colorbar
            if log.n_choices==3
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
            if log.n_choices==3
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
            for c=1:log.n_choices
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+(t-1)*3+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(squeeze(s_ntrialmatrix{log.n_subjs+1}(t,c,:,:)), log.range_ntrials);  axis off; axis square; colorbar
            end
            t=2;
            subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+5,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            text(0.3,0.5, 'Mean no. of trials per cell' ,'FontSize', 8); axis off
            for c=1:log.n_choices
                subtightplot(log.n_rois+f.plotrowsadd, f.subplotcols,  (log.n_rois+1)*f.subplotcols+f.plotspace+c+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
                imagesc(squeeze(s_ntrialmatrix{log.n_subjs+1}(t,c,:,:)), log.range_ntrials);  axis off; axis square; colorbar
            end
        end
    end
    f.figheight=log.n_rois*250+150;   f.figheight=1000; f.nancol=[0 0 0];  % Betas x roi (combine tasks)
    if log.n_choices==3; f.figwidth= 1400; f.subplotcols=4; f.plotspace=3; f.subplot_VerHorz=[0.01 0.005];  f.fig_BotTop=[0.001 0.025];  f.fig_LeftRight=[0.15 0.01]; f.nancol=[0 0 0];
    else f.figwidth= 700; f.subplotcols=2;  f.subplot_VerHorz=[0.01 0.005]; f.fig_BotTop=[0.001 0.025]; f.fig_LeftRight=[0.3 0.3];
    end
    fm=figure('Name', [request.FLmodel ': Mean betas for each ROI, combined across tasks    (min. no. subjects per cell = ' num2str(request.minsubs_cell) ')'], 'NumberTitle', 'off', 'Position',[200,70,f.figwidth,f.figheight]); set(gcf,'Color',[1 1 1]);
    for r=1:log.n_rois
        subtightplot(log.n_rois, f.subplotcols,  (r-1)*f.subplotcols+1,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
        text(0.3,0.5, log.rois{r},'FontSize', 15); axis off
        for c=1:log.n_choices
            subtightplot(log.n_rois, f.subplotcols,  (r-1)*f.subplotcols+1+c,  f.subplot_VerHorz,  f.fig_BotTop, f.fig_LeftRight);
            imagescnan(m_roimatrix_comb(:,:,c,r), 'NaNColor', f.nancol); axis off; axis square; colorbar
            if r==1 && log.n_choices==1; title('Combined cF ct');  elseif r==1; title(['cFct ' log.choicecF{c}(1:end-1)]); end
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
            wr.cf=[]; wr.ct=[]; r_anova{r+1,1}= log.rois{r};
            
            % Compile sample
            for ee=1:6
                e=log.n_levels+1-ee;
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

end


