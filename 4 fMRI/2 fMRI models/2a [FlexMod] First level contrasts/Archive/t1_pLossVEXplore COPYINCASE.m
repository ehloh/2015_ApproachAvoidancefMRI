function [Con log] = t1_pLossVEXplore(where, log)
% [Con log] = t1_pLossVEXplore(where, log)
% Set up 1st-level contrasts for: 
%         (a) pLoss (both, cF, ct, cF-ct, ct-cF)
%         (b) VExplore (both, cF, ct, cF-ct, ct-cF)
%         (c) pLoss-VEXplore (both, cF, ct)
%         (d) VExplore-pLoss (both, cF, ct)
%
% Weights applied are according to their real specification 
% in the CONFLICT behavioural models. Control task trials still
% follow the Conflict parameter spec. Choice is ignored in contrasts,
% but spec+est model already controls for choice. Script can be 
% adapted for other RL variables in the future, but only 2 pars can be
% compared
%
% Second level models included:
%   (1) Task x Explore (2x2 FLEXIBLE Factorial)
%   (2) Add sensible contrasts to the 2x2 (Task x Explore) Factorial
%
% Output 'Con.ConInstruc' follows specifications of RegWeights output from the
% function xfx_2GetRegweightsFromSearchCriteria, details the (subject-specific) 
% regerssor weights for each contrast for each subject
%
% ---------------------------------------------------------------------------------------

% Check: Is this the correct Contrasts table for the chosen onsets model?
onsetsmodels4thiscontrasttype={'m_t1_Trialtype'};
if sum(strcmp(log.onsetsmodel,onsetsmodels4thiscontrasttype))~=1
    error('Wrong contrasts type chosen for this onsets model!')
end

%% (1) General changeable settings

% Folders
where.subFLest=cell(log.n_subjs,1);
where.subFLfol=cell(log.n_subjs,1);
for s=1:log.n_subjs
    where.subFLest=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Estimated' filesep];
    where.subFLfol=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel(1:3) '     ' log.firstlevel_contraststype filesep];
end

%% (2) Compile instructions for first-level contrasts (apply parameter values)

if  log.execute_contrasts==1
    % Get parameter values (read from function used to run computational models)
    col.EnvThreat=1;
    col.NTokens=2;
    col.pLoss=3;
    col.Entropy=4;
    col.VExplore=5;
    col.EV=6;
    col.OutcomeMagnitude=7;
    col.OutcomeMean=8;
    col.OutcomeVariance=9;
    tasks={'cF'; 'ct'};
    parvals=nan*zeros(6*6,9);
    parvals(:,1)=ceil((1:36)/6)';
    parvals(:,2)=2*repmat((1:6)',6,1);
    [ parvals] = fpar_conflict( parvals, col);
    parvals(:,col.NTokens)=parvals(:,col.NTokens)/2;
    
    % (1) Build instructions for each parameter
    instruc.singlepar={'pLoss'; 'VExplore'}; k=1;
    for p=1:length(instruc.singlepar) % Combined
        SearchCriteria{k,1}=[instruc.singlepar{p} '_both'];
        %
        wp.par=[num2cell(parvals(:,1:6)) cell(6*6,2)];
        wp.par(:,8)=cellfun(@(x,y)[{['cF_t' num2str(x) '-' num2str(y)] ; ['ct_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
        eval(['wp.par(:,7)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
        %
        SearchCriteria{k,2}=wp.par(:,7:8); SearchCriteria{k,3}=p;
        wp=[]; k=k+1;
    end
    for t=1:2 % Single task
        for p=1:length(instruc.singlepar)
            SearchCriteria{k,1}=[instruc.singlepar{p} '_' tasks{t}];
            %
            wp.par=[num2cell(parvals(:,1:6)) cell(6*6,2)];
            wp.par(:,8)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
            eval(['wp.par(:,7)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
            %
            SearchCriteria{k,2}=wp.par(:,7:8); SearchCriteria{k,3}=p;
            wp=[]; k=k+1;
        end
    end
    for t=1:2 % Task comparisons
        for p=1:length(instruc.singlepar)
            SearchCriteria{k,1}=[instruc.singlepar{p} '_' tasks{t} '-' tasks{3-t}];
            wp.par=[num2cell(parvals(:,1:6)) cell(6*6,4)];
            
            % Positive task
            wp.par(:,8)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
            eval(['wp.par(:,7)=wp.par(:, col.'  instruc.singlepar{p,1} ');'])
            
            % Negative task
            wp.par(:,10)=cellfun(@(x,y)[{[tasks{3-t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
            eval(['wp.par(:,9)= num2cell(cell2mat(wp.par(:, col.'  instruc.singlepar{p,1} '))*-1);'])
            
            %
            SearchCriteria{k,2}=vertcat(wp.par(:,7:8),wp.par(:,9:10)); SearchCriteria{k,3}=p;
            wp=[]; k=k+1;
        end
    end
    SearchCriteria=sortrows(SearchCriteria,3); SearchCriteria(:,3)=num2cell(0);
    
    % (2) Build instructions comparing parameters
    for p=1:length(instruc.singlepar)  % Combined across task
        SearchCriteria{k,1}=[instruc.singlepar{p} '-' instruc.singlepar{3-p} '_both'];
        wp.par=[num2cell(parvals(:,1:6)) cell(6*6,4)];
        
        % Positive parameter
        wp.par(:,8)=cellfun(@(x,y)[{['cF_t' num2str(x) '-' num2str(y)] ; ['ct_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
        eval(['wp.par(:,7)=wp.par(:, col.'  instruc.singlepar{p} ');'])
        
        % Negative parameter
        wp.par(:,10)=cellfun(@(x,y)[{['cF_t' num2str(x) '-' num2str(y)] ; ['ct_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
        eval(['wp.par(:,9)= num2cell(cell2mat(wp.par(:, col.'  instruc.singlepar{3-p} '))*-1);'])
        
        %
        SearchCriteria{k,2}=vertcat(wp.par(:,7:8),wp.par(:,9:10)); SearchCriteria{k,3}=p;
        wp=[]; k=k+1;
    end
    for t=1:2 % Single task
        for p=1:length(instruc.singlepar)
            SearchCriteria{k,1}=[instruc.singlepar{p} '-' instruc.singlepar{3-p} '_' tasks{t}];
            wp.par=[num2cell(parvals(:,1:6)) cell(6*6,4)];
            
            % Positive parameter
            wp.par(:,8)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
            eval(['wp.par(:,7)=wp.par(:, col.'  instruc.singlepar{p} ');'])
            
            % Negative parameter
            wp.par(:,10)=cellfun(@(x,y)[{[tasks{t} '_t' num2str(x) '-' num2str(y)]}], wp.par(:,1), wp.par(:,2),'UniformOutput', 0);
            eval(['wp.par(:,9)= num2cell(cell2mat(wp.par(:, col.'  instruc.singlepar{3-p} '))*-1);'])
            
            %
            SearchCriteria{k,2}=vertcat(wp.par(:,7:8),wp.par(:,9:10)); SearchCriteria{k,3}=p;
            wp=[]; k=k+1;
        end
    end
    SearchCriteria=sortrows(SearchCriteria,3);SearchCriteria(:,3)=[];
    Con.SearchCriterias=SearchCriteria;
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
        [ConInstruc{s}] = xfx_2GetRegweightsFromSearchCriteria(SearchCriterias, RegLists{s});
    end
end

%% (4) Generate batch for First-level contrasts

if  log.execute_contrasts==1
    
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel(1:4) ' Contrasted   ' log.firstlevel_contraststype filesep];
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


error('stopped here - haven''t set up 2nd level models yet')

% 
% % Fetch details from the contrasted first level model
% disp(' ################ Prepping for 2nd level ##################')
% disp('Fetching details from (already-contrasted) first level model - sampling 1st subject')
% c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.onsetsmodel(1:4) ' Contrasted   ' log.firstlevel_contraststype filesep 'SPM.mat']);
% for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
% disp(Con.ConNames);
% 
% % Models to run? 
% %   1: Task x Explore factorial
% %   2: Add sensible contrasts to the 2x2 (Task x Explore) Factorial
% Nmodels=1;
% Analysis=cell(Nmodels,3);
% 
% % ################################################################
% 
% % [Model #1] Task x Explore ################
% for o1=1:1 
%     m=1; Analysis{m,1}='TaskxExplore';
%     Analysis{m,2}='NFlexFactorial';
%     Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
%     
%     % Specify factorial design
%     Analysis{m,3}.NFactorialCells=4;
%     Analysis{m,3}.Design{1}.name='FacCell';      % Factor 1: NCells (2x2=4)
%     Analysis{m,3}.Design{1}.dept=1;
%     Analysis{m,3}.Design{1}.variance=1;
%     Analysis{m,3}.Design{1}.gmsca=0;
%     Analysis{m,3}.Design{1}.ancova=0;
%     Analysis{m,3}.Design{2}.name='Subject';      % Factor 2: Subjects
%     Analysis{m,3}.Design{2}.dept=0;
%     Analysis{m,3}.Design{2}.variance=0;
%     Analysis{m,3}.Design{2}.gmsca=0;
%     Analysis{m,3}.Design{2}.ancova=0;
%     
%     % Factors to include in model 
%     Analysis{m,3}.ModelFac=[1 2]; % Factor no.s (corresponding to above, i.e. factor 2=Subject)  
%     
%     % Assign contrast images to correct design cells (via search instruction)
%     Analysis{m,3}.Cells=cell(2*2, 2); 
%     Analysis{m,3}.Cells={'cF_Explore'             1 ; % Cell assignment can also be multi-factorial, e..g [1 2]
%                                       'cF_NonExplore'       2;
%                                       'ct_Explore'              3;
%                                       'ct_NonExplore'        4;};
%     Analysis{m,3}.Cells_FacAssignment=cell2mat(Analysis{m,3}.Cells(:,2));
%     for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
%         if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
%             Analysis{m,3}.Cells{i,3}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
%         else
%             error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
%         end
%     end
%     
%     % Subject covariates (Exploration parameters from model fits)
%     modelname='val2_b_fecx';
%     explore_pars=nan*zeros(log.n_subjs,2);
%     for s=1: log.n_subjs
%         ws=SubjectsOK{s,2}.b.cF_modelpars.pars;
%         %
%         explore_pars(s,1)=ws{find(strcmp(ws(:,1), modelname)==1),7};
%         explore_pars(s,2)=ws{find(strcmp(ws(:,1), modelname)==1),8};
%         %
%         ws=[];    
%     end
%     Analysis{m,3}.Covar(1).Type='single'; % 1st exploration parameter
%     Analysis{m,3}.Covar(1).Ncells=4; 
%     Analysis{m,3}.Covar(1).CovSingle_Vector=explore_pars(:,1);
%     Analysis{m,3}.Covar(1).CovSingle_Name='Explore1';
%     Analysis{m,3}.Covar(2).Type='single'; % 2nd exploration parameter
%     Analysis{m,3}.Covar(2).Ncells=4; 
%     Analysis{m,3}.Covar(2).CovSingle_Vector=explore_pars(:,2);
%     Analysis{m,3}.Covar(2).CovSingle_Name='Explore2';
%    
%     % #############################################
%     for o2=1:1 % [DISUSED] Using the [Non-flexible] Full factorial
% %     Analysis{m,3}.Design{1}.name='Task';      % Factor 1: Task
% %     Analysis{m,3}.Design{1}.levels=2;
% %     Analysis{m,3}.Design{1}.dept=0;
% %     Analysis{m,3}.Design{1}.variance=1;
% %     Analysis{m,3}.Design{1}.gmsca=0;
% %     Analysis{m,3}.Design{1}.ancova=0;
% %     Analysis{m,3}.Design{2}.name='Explore';      % Factor 2: Explore
% %     Analysis{m,3}.Design{2}.levels=2;
% %     Analysis{m,3}.Design{2}.dept=0;
% %     Analysis{m,3}.Design{2}.variance=1;
% %     Analysis{m,3}.Design{2}.gmsca=0;
% %     Analysis{m,3}.Design{2}.ancova=0;
% %     
% %     % Assign contrast images to correct design cells (via search instruction)
% %     Analysis{m,3}.Cells=cell(2*2, 3); 
% %     Analysis{m,3}.Cells={'cF_Explore'             1 1;
% %                                       'cF_NonExplore'       1 2;
% %                                       'ct_Explore'              2 1;
% %                                       'ct_NonExplore'        2 2;};
% %     for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
% %         if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
% %             Analysis{m,3}.Cells{i,1}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
% %         else
% %             error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
% %         end
% %     end
% %     
% %     % Subject covariates (Exploration parameters from model fits)
% %     modelname='val2_b_fecx';
% %     explore_pars=nan*zeros(log.n_subjs,2);
% %     for s=1: log.n_subjs
% %         ws=SubjectsOK{s,2}.b.cF_modelpars.pars;
% %         %
% %         explore_pars(s,1)=ws{find(strcmp(ws(:,1), modelname)==1),7};
% %         explore_pars(s,2)=ws{find(strcmp(ws(:,1), modelname)==1),8};
% %         %
% %         ws=[];    
% %     end
% %     Analysis{m,3}.Covar(1).Type='single';
% %     Analysis{m,3}.Covar(1).Ncells=4; % 1st exploration parameter
% %     Analysis{m,3}.Covar(1).CoVector=explore_pars(:,1);
% %     Analysis{m,3}.Covar(1).CovNames=cellstr([char(repmat({'explore1_'}, log.n_subjs,1)) char(log.subjects)]);    
% %     Analysis{m,3}.Covar(2).Ncells=4; % 2nd exploration parameter
% %     Analysis{m,3}.Covar(2).CoVector=explore_pars(:,2);
% %     Analysis{m,3}.Covar(2).CovNames=cellstr([char(repmat({'explore2_'}, log.n_subjs,1)) char(log.subjects)]);    
%     
%     end
% end
% 
% % [Model #2] Add sensible contrasts to the Task x Explore factorial
% for o1=1:1 
%     m=m+1; Analysis{m,1}='Adding Factorial contrasts to TaskxExplore';
%     Analysis{m,2}='AddContrasts';
%     %
%     Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep 'TaskxExplore'];
%     Analysis{m,3}.DeleteOldCon=1; % There aren't any existing contrasts to begin with, using FlexFactorial
%       
%     % Weights for subject covariates
%     %       covariate weights= ones(1, log.n_subjs)*(no. of regressors weighted)/(log.n_subjs)
%     %       e.g. [1 1 0 0 ones(1,13)*(2)/13]
%     SubCovarWeights=ones(1,log.n_subjs)/(log.n_subjs);
%     SubCovarsForExplorePar=[0 0];
%     
%     
%     % Contrast instructions
%     Analysis{m,3}.ConInstruc={'ME_Task'                                 [1 1 -1 -1]         2;                                  % Factorial F contrasts
%                                                 'ME_ExploreOr'                          [1 -1 1 -1]         2;
%                                                 'Interaction_TaskxExpploreOr'    [1 -1 -1 1]         2;
%                                                 
%     
%                                                 % T-contrasts
%                                                 'ConflictTask'                [1 1 0 0   SubCovarsForExplorePar      SubCovarWeights*2]           1;     % Main factor means
%                                                 'ControlTask'                [0 0 1 1   SubCovarsForExplorePar      SubCovarWeights*2]           1;
%                                                 'Explore'                       [1 0 1 0   SubCovarsForExplorePar      SubCovarWeights*2]           1;
%                                                 'NonExplore'                [0 1 0 1   SubCovarsForExplorePar      SubCovarWeights*2]           1;
%                                                 'cF_Explore'                [1 0 0 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;        % Cell means
%                                                 'ct_Explore'                [0 0 1 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;
%                                                 'cF_NonExplore'           [0 1 0 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;
%                                                 'ct_NonExplore'             [0 0 0 1   SubCovarsForExplorePar      SubCovarWeights*1]           1;
%                                                 
%                                                 % Contrast comparisons (t-test) do not incldue subject-covariate weightings
%                                                 'Task_cF-ct'             [1 1 -1 -1]           1;        % Factor effects
%                                                 'Task_ct-cF'             [-1 -1 1 1]           1;
%                                                 'Explore-Non'             [1 -1 1 -1]           1;
%                                                 'Non-Explore'             [-1 1 -1 1]           1;
%                                                 'cF_Explore-Non'        [1 -1 0 0]           1;   % Cell effects
%                                                 'cF_Non-Explore'        [-1 1 0 0]           1;
%                                                 'ct_Explore-Non'        [0 0 1 -1]           1;
%                                                 'ct_Non-Explore'        [0 0 -1 1]           1;
%                                                 'cF_Explore_vs'         [3 -1 -1 -1]           1;   % Cell versus
%                                                 'cF_NonExplore_vs'    [-1 3 -1 -1]           1;
%                                                 'ct_Explore_vs'          [-1 -1 3 -1]           1;
%                                                 'ct_NonExplore_vs'     [-1 -1 -1 3]           1;
%                                                 'TaskxExp_Inter1'       [1 -1 -1 1]           1; % Interaction terms
%                                                 'TaskxExp_Inter2'       [-1 1 1 -1]           1;};
%                                             
%                                             
% end
% 
% % ################################################################
% 
% % Output to main script
% Con.SecondLevelAnalysis=Analysis;

end

