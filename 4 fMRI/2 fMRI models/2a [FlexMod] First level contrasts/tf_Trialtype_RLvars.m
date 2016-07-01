function [Con log] = tf_Trialtype_RLvars(where, log)
% Set up 1st-level contrasts for all RL variables. For use with Trialtype
% models. Assumes regressors are labelled e.g. cF_t1-1
%
%   RL variables:  
%         - EntropyNTok (O)
%         - VExplore (V) 
%         - Entropy (U)
%         - EV (L)
%         - pLoss (P)
%         - EnvThreat (E)
%         - NTokens (N)
%         - EntropyNTok (O)
%
%
%   For each variable, set up: cF only, ct only, both tasks, cF-ct, ct-cF
%
% Weights applied are according to their real specification in the behavioural models
% See fpar_conflict and fpar_control for details. No choice regressors 
%
% Second level models included:
%   (1) Task Anova (1-way FLEXIBLE Factorial)
%   (2) Add sensible contrasts to the Factorial Anova
%   - These 2nd-level models are set up for each RL variable, separately
%
% ---------------------------------------------------------------------------------------

% Which RL variables to include?
RLvariables_same={'Entropy'; 'pLoss'; 'EnvThreat'; 'NTokens'; 'EntropyNTok'}; % Specified the same across both tasks
RLvariables_diff={'VExplore'; 'EV'};
RLvariables=vertcat(RLvariables_same,RLvariables_diff);

%% (1) General changeable settings

% Folders
where.subFLest=cell(log.n_subjs,1);
for s=1:log.n_subjs
    where.subFLest=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated' filesep];
end

% Where are first level contrasts?
log.firstlevelfolder=['m_' log.firstlevel_contraststype ' Contrasted' filesep];

%% (2) Compile instructions for first-level contrasts (apply parameter values)

if  log.execute_contrasts==1
    
    % Get parameter values (read from function used to run computational models)
    col.TrialType_EnvThreat=1;  % Columns marking EnvThreat and NTokens levels (for reading Trial-type name). Other columns are for building instructions (weights)
    col.TrialType_NTokens=2;
    col.EnvThreat=3;
    col.NTokens=4;
    col.pLoss=5;
    col.Entropy=6;
    col.VExplore=7;
    col.EV=8;
    col.EntropyNTok=9;
    tasks={'cF'; 'ct'};
    pars=nan*zeros(6*6,9);
    pars(:,col.EnvThreat)=ceil((1:36)/6)';
    pars(:,col.NTokens)=2*repmat((1:6)',6,1);
    pars(:,col.TrialType_EnvThreat)=pars(:,col.EnvThreat);
    pars(:,col.TrialType_NTokens)=pars(:,col.NTokens)/2;
    
    % Set up parameters + mean centre params (except first 2 cols, assumed trialtype markers)
    [ parvals_cF] = fpar_conflict( pars, col);
    parvals_cF=[parvals_cF(:,[1 2])  parvals_cF(:, 3:end)- repmat(mean(parvals_cF(:, 3:end)), 36, 1)];
    [ parvals_ct] = fpar_control( pars, col); % Control-task parvals_cF - only used for EV
    parvals_ct=[parvals_ct(:,[1 2])  parvals_ct(:, 3:end)- repmat(mean(parvals_ct(:, 3:end)), 36, 1)];
   
    % (A) Build instructions for parameters that are specified the same across tasks ####################
    instruc.singlepar=RLvariables_same; k=1;
    for p=1:length(instruc.singlepar) % Combined  
        SearchCriteria{k,1}=[instruc.singlepar{p} '_both'];
        %
        wp.par=[num2cell(parvals_cF) cell(6*6,2)]; j=size(wp.par,2);
        wp.par(:,j)=cellfun(@(x,y)[{['cF_t' num2str(x) '-' num2str(y)] ; ['ct_t' num2str(x) '-' num2str(y)]}], wp.par(:,col.TrialType_EnvThreat), wp.par(:,col.TrialType_NTokens),'UniformOutput', 0);
        eval(['wp.par(:,j-1)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
        %
        SearchCriteria{k,2}=wp.par(:,j-1:j); SearchCriteria{k,3}=p;
        wp=[]; k=k+1;
    end
    
    for t=1:2 % Single task
        for p=1:length(instruc.singlepar)
            SearchCriteria{k,1}=[instruc.singlepar{p} '_' tasks{t}];
            %
            wp.par=[num2cell(parvals_cF) cell(6*6,2)]; j=size(wp.par,2);
            wp.par(:,j)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,col.TrialType_EnvThreat), wp.par(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par(:,j-1)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
            %
            SearchCriteria{k,2}=wp.par(:,j-1:j); SearchCriteria{k,3}=p;
            %
            wp=[]; k=k+1;
        end
    end
    for t=1:2 % Task comparisons  
        for p=1:length(instruc.singlepar)
            SearchCriteria{k,1}=[instruc.singlepar{p} '_' tasks{t} '-' tasks{3-t}];
            wp.par=[num2cell(parvals_cF) cell(6*6,4)]; j=size(wp.par,2);
            
            % Positive task
            wp.par(:,j-2)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,col.TrialType_EnvThreat), wp.par(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par(:,j-3)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
            
            % Negative task
            wp.par(:,j)=cellfun(@(x,y)[{[tasks{3-t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,col.TrialType_EnvThreat), wp.par(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par(:,j-1)= num2cell(cell2mat(wp.par(:, col.'  instruc.singlepar{p,1} '))*-1);'])
            
            %
            SearchCriteria{k,2}=vertcat(wp.par(:,j-3:j-2),wp.par(:,j-1:j)); SearchCriteria{k,3}=p;
            wp=[]; k=k+1;
        end
    end
    SearchCriteria=sortrows(SearchCriteria,3); SearchCriteria(:,3)=num2cell(0);
    
    % (B) Build instructions for parameters that are specified differently across tasks ####################
    instruc.splitpar=RLvariables_diff;
    for p=1:length(instruc.splitpar)  % Combined  
        SearchCriteria{k,1}=[instruc.splitpar{p} '_both'];
        
        % Conflict task parameters
        wp.par_cF=[num2cell(parvals_cF) cell(6*6,2)]; j=size(wp.par_cF,2);
        wp.par_cF(:,j)=cellfun(@(x,y)[{['cF_t' num2str(x) '-' num2str(y)]}], wp.par_cF(:,col.TrialType_EnvThreat), wp.par_cF(:,col.TrialType_NTokens),'UniformOutput', 0);
        eval(['wp.par_cF(:,j-1)=wp.par_cF(:, col.'  instruc.splitpar{p,1} ');'])
        
        % Control task parameters
        wp.par_ct=[num2cell(parvals_ct) cell(6*6,2)]; jj=size(wp.par_ct,2);
        wp.par_ct(:,jj)=cellfun(@(x,y){['ct_t' num2str(x) '-' num2str(y)]}, wp.par_ct(:,col.TrialType_EnvThreat), wp.par_ct(:,col.TrialType_NTokens),'UniformOutput', 0);
        eval(['wp.par_ct(:,jj-1)=wp.par_ct(:, col.'  instruc.splitpar{p,1} ');'])
        
        %
        SearchCriteria{k,2}=vertcat(wp.par_cF(:,j-1:j),wp.par_ct(:,jj-1:jj)); SearchCriteria{k,3}=p;
        wp=[]; k=k+1;
    end
    for t=1:2 % Single task  
        for p=1:length(instruc.splitpar)
            SearchCriteria{k,1}=[instruc.splitpar{p} '_' tasks{t}];
            %
            switch t
                case 1;  wp.par=[num2cell(parvals_cF) cell(6*6,2)]; j=size(wp.par,2);
                case 2;  wp.par=[num2cell(parvals_ct) cell(6*6,2)]; j=size(wp.par,2);
            end
            %
            wp.par(:,j)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,col.TrialType_EnvThreat), wp.par(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par(:,j-1)=wp.par(:, col.'  instruc.splitpar{p,1} ');'])
            %
            SearchCriteria{k,2}=wp.par(:,j-1:j); SearchCriteria{k,3}=p;
            %
            wp=[]; k=k+1;
        end
    end
    for t=1:2 %  Task comparisons  
        for p=1:length(instruc.splitpar)
            SearchCriteria{k,1}=[instruc.splitpar{p} '_' tasks{t} '-' tasks{3-t}];
            switch t
                case 1;  
                    wp.par1=[num2cell(parvals_cF) cell(6*6,2)]; j=size(wp.par1,2);
                    wp.par2=[num2cell(parvals_ct) cell(6*6,2)]; jj=size(wp.par2,2);
                case 2;  
                    wp.par1=[num2cell(parvals_ct) cell(6*6,2)]; j=size(wp.par1,2);
                    wp.par2=[num2cell(parvals_cF) cell(6*6,2)]; jj=size(wp.par2,2);
            end
            
            % Positive task
            wp.par1(:,j)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par1(:,col.TrialType_EnvThreat), wp.par1(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par1(:,j-1)=wp.par1(:, col.'  instruc.splitpar{p,1} ');'])
            
            % Negative task
            wp.par2(:,jj)=cellfun(@(x,y)[{[tasks{3-t} '_t' num2str(x) '-' num2str(y)]}], wp.par2(:,col.TrialType_EnvThreat), wp.par2(:,col.TrialType_NTokens), 'UniformOutput', 0);
            eval(['wp.par2(:,jj-1)=wp.par2(:, col.'  instruc.splitpar{p,1} ');'])
            
            %
            SearchCriteria{k,2}=vertcat(wp.par1(:,j-1:j), wp.par2(:,jj-1:jj)); SearchCriteria{k,3}=p;
            %
            wp=[]; k=k+1;
        end
    end
    SearchCriteria(:,3)=[];
    
    % Forward to building contrasts
    Con.SearchCriteria=SearchCriteria;
end


%% (3) Get regressor weights for First-level contrast
% Search criteria is fed into function xfx_2GetRegweightsFromSearchCriteria,
% to derive regressor weights for each contrast. Format of search criteria
% must therefore follow format specified in this search function.

if  log.execute_contrasts==1
    disp('####### Applying search criteria to derive regressor weights for contrasts  ##########'); disp(' ');
    RegLists=cell(log.n_subjs,1); ConInstruc=cell(log.n_subjs,1);
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
        [RegLists{s}] = xfx_1LoadSimpleReglist(where, log,s); % Load (edited) list of regressors
        [ConInstruc{s}] = xfx_2GetRegweightsFromSearchCriteria(SearchCriteria, RegLists{s});
    end
end

%% (4) Generate batch for First-level contrasts

if  log.execute_contrasts==1
    
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.firstlevelfolder];
        [matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,log, ws, s,ConInstruc{s});
    end
    
    % Record/output details
    Con.RegLists=RegLists;
    Con.SearchCriteria=SearchCriteria;
    Con.ConInstruc=ConInstruc;
    Con.matlabbatch=matlabbatch;
end

%% (5) Set up the 2nd level models (executed next) that go with this model
% Establish instructions for 2nd level analyses to run ################
% Col 1=Name of 2nd-level model, Col 2=Which analysis function, Col 3=Inputs to function
% Instructions are executed in main script to run all possible 2nd level models


% keyboard

% Fetching details from (already-contrasted) first level model - sampling 1st subject
disp(' ################ Prepping for 2nd level ##################')
c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.firstlevelfolder 'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run, for each RL variable ----------------
%   1: Task x Variable flexible factorial
%   2: Add sensible contrasts to the 2x2 (Task x Variable) Factorial
Nmodels=length(RLvariables);
Analysis=cell(Nmodels,3); m=0;

% ################################################################

for p=1:length(RLvariables)
    
    % [Model #P-1a] Task x Explore ################
    for o1=1:1
        m=m+1; Analysis{m,1}=[RLvariables{p} ' TaskAnova'];
        Analysis{m,2}='NFlexFactorial';
        Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
        
        % Specify factorial design
        Analysis{m,3}.NFactorialCells=4;
        Analysis{m,3}.Design{1}.name='Task';      % Factor 1: NCells (2x2=4)
        Analysis{m,3}.Design{1}.dept=0;
        Analysis{m,3}.Design{1}.variance=1;
        Analysis{m,3}.Design{1}.gmsca=0;
        Analysis{m,3}.Design{1}.ancova=0;
        Analysis{m,3}.Design{2}.name='Subject';      % Factor 2: Subjects
        Analysis{m,3}.Design{2}.dept=0;
        Analysis{m,3}.Design{2}.variance=0;
        Analysis{m,3}.Design{2}.gmsca=0;
        Analysis{m,3}.Design{2}.ancova=0;
        Analysis{m,3}.Covar=[]; % No covariates for now
        
        % Factors to include in model
        Analysis{m,3}.ModelFac=[1 2]; % Factor no.s (corresponding to above, i.e. factor 2=Subject)
        
        % Assign contrast images to correct design cells (via search instruction)
        Analysis{m,3}.Cells={[RLvariables{p} '_cF']             1 ; % Cell assignment can also be multi-factorial, e..g [1 2]
                                            [RLvariables{p} '_ct']             2 ; };
        Analysis{m,3}.Cells_FacAssignment=cell2mat(Analysis{m,3}.Cells(:,2));
        for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
            if length(  find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1}))   )==1
                Analysis{m,3}.Cells{i,3}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
            else
                error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
            end
        end
        
    end
    
    % [Model #P-1b] Add sensible contrasts to the Task x Explore factorial
    for o1=1:1
        m=m+1; Analysis{m,1}=['Adding Factorial contrasts to ' RLvariables{p} ' TaskAnova'];
        Analysis{m,2}='AddContrasts';
        %
        Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep RLvariables{p} ' TaskAnova'];
        Analysis{m,3}.DeleteOldCon=1; % There aren't any existing contrasts to begin with, using FlexFactorial
        
        % Weights for subject covariates
        %       covariate weights= ones(1, log.n_subjs)*(no. of regressors weighted)/(log.n_subjs)
        %       e.g. [1 1 0 0 ones(1,13)*(2)/13]
        SubCovarWeights=ones(1,log.n_subjs)/(log.n_subjs);
        
        % Contrast instructions
        Analysis{m,3}.ConInstruc={['ME ' RLvariables{p} ' all']        [1 1 SubCovarWeights*2]                                   2;
                                                    'ME task'                                  [1 -1]                                2;
                                                    'cF Only'                                   [1 0 SubCovarWeights]       1;
                                                    'ct Only'                                   [0 1 SubCovarWeights]       1;
                                                    'cF-ct'                                      [1 -1]                                 1;
                                                    'ct-cF'                                      [-1 1]                                 1;
                                                    'Identity'                                  [1 0 SubCovarWeights;0 1 SubCovarWeights]                                   2;
                                                    };
    end
    
    % [Model #P-2a] One-sample ttest for the variable itself ################
    for o1=1:1
        m=m+1; Analysis{m,1}=[RLvariables{p} '_both_onesampleT'];
        Analysis{m,2}='onesampleT';
        %
        Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
        Analysis{m,3}.Covar=[];
        if length(     find(strcmp(Con.ConNames(:,1), [RLvariables{p} '_both']))    )==1 % Fetch correct con image
            Analysis{m,3}.ConImage=Con.ConNames{find(strcmp(Con.ConNames(:,1), [RLvariables{p} '_both'])),2};
        else  error('When fetching contrast images for 1 sampleT, # of con-image matches for a particular cell ~=1')
        end
    end
    
    % [Model #P-2b] Add 2nd-level contrasts to the one-sample ttest
    for o1=1:1
        m=m+1; Analysis{m,1}=['Adding 2nd-level contrasts to ' RLvariables{p} ' onesampleT'];
        Analysis{m,2}='AddContrasts';
        %
        Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep RLvariables{p} '_both_onesampleT'];
        Analysis{m,3}.DeleteOldCon=1; % There aren't any existing contrasts to begin with, using FlexFactorial
        
        % Contrast instructions
        Analysis{m,3}.ConInstruc={'Pos'     1       1;
                                                   'Neg'     -1       1;};
    end
    
end

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

% Save details to master 2nd-level folder
details.RLvariables=RLvariables; 
details.log=log;
save([where.firstlevel_resultsfolder filesep 'details_2ndlevelmodel.mat'], 'details')


end

