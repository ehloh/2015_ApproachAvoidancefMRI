function [Con log] = f4_ExploreOr_MovWind(where, log)
% Compare Explore vs NonExplore, including all cells where % Explore is
% between 35-75%. Include all subjects for which no. qualifying cells <1
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

%% (1) General changeable settings

settings.range_explore=[0.35 0.75];
settings.min_numcells=1;

%% (2) Generate subject-specific search criteria (First level contrasts)
% Search criteria is fed into function xfx_2GetRegweightsFromSearchCriteria,
% to derive regressor weights for each contrast. Format of search criteria
% must therefore follow format specified in this search function.

disp('########### Identifying included subjects + specifying subject-specific search criteria (to find regressors for weighting) #####################')
SubjectsOK=cell(log.n_subjs,5); % Col 2=Behavioural stats,  Col 3=N cells included, Col 4=Include subject, Col 5= SearchCriteria
RegLists=cell(log.n_subjs,1); 

for s=1:log.n_subjs % Check subject stats + compile search criteria
    disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
    [RegLists{s}] = xfx_1LoadSimpleReglist(where, log,s); % Load (edited) list of regressors
    ws.b=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.subjects{s} '_behstats.mat']);
    if s==1; input('Accessing variable ''_behstats.mat'' in FL folder. Is this correct for all subjects? '); end
        
    % Log Percentage
    ws.Conflict_p=cell(6,6); ws.Conflict_include=nan*zeros(6,6);
    ws.Control_p=cell(6,6); ws.Control_include=nan*zeros(6,6);
    for e=1:6 
        for n=1:6
            % Conflict task
            ws.Conflict_p{e,n}(1)=ws.b.cF_stats{e,n}.p_Explore;
            ws.Conflict_p{e,n}(2)=ws.b.cF_stats{e,n}.p_Accept+ws.b.cF_stats{e,n}.p_Reject;
            if ws.Conflict_p{e,n}(1)>=settings.range_explore(1) && ws.Conflict_p{e,n}(1)<=settings.range_explore(2)
                ws.Conflict_include(e,n)=1;
            else
                ws.Conflict_include(e,n)=0;
            end
            
            % Control task
            ws.Control_p{e,n}(1)=ws.b.ct_stats{e,n}.p_Explore;
            ws.Control_p{e,n}(2)=ws.b.ct_stats{e,n}.p_NoBomb+ws.b.ct_stats{e,n}.p_Bomb;
            if ws.Control_p{e,n}(1)>=settings.range_explore(1) && ws.Control_p{e,n}(1)<=settings.range_explore(2)
                ws.Control_include(e,n)=1;
            else
                ws.Control_include(e,n)=0;
            end
        end
    end
    % NOTE: EnvThreat is NOT in the same format as Trialtype
    % (e.g. cF_Accept_t1-2 is the LOWEST EnvThreat. But in beh_stats, it
    % is in visualization space (i.e. cell 1,1=Highest thread, 2
    % Tokens). This need to be corrected before the Trial-type search
    % criteria is constructed
    
    % Is this subject ok?
    SubjectsOK{s,1}=log.subjects{s};
    SubjectsOK{s,2}=ws;
    SubjectsOK{s,3}=sum(sum(ws.Conflict_include));
    if SubjectsOK{s,3}>=settings.min_numcells
        SubjectsOK{s,4}=1;
    else SubjectsOK{s,4}=0;
    end
    
    % ####### Format for the search function, if this subject will be included ###########
    if  log.execute_contrasts==1 && SubjectsOK{s,4}==1
        
        % Which cells?
        [ws.Cells2Include(:,1) ws.Cells2Include(:,2)]=find(ws.Conflict_include);
        ws.Cells2Include(:,1)=7-ws.Cells2Include(:,1); % Convert EnvThreat from Visualization space to Trial-type space
        
        % Construct search criteria!
        ws.Crit=cell(4,2); c=1;

        % [Contrast 1] Conflict task, Explore ################
        ws.Crit{c,1}='cF_Explore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(ws.Cells2Include,1)
            tt=['_t' num2str(ws.Cells2Include(i,1)) '-' num2str(ws.Cells2Include(i,2))];
            %
            ws.Crit{c,2}{2}{v}=['cF_Explore' tt];
            v=v+1;
        end
        c=c+1;
        
        % [Contrast 2] Conflict task, NonExplore ################
        ws.Crit{c,1}='cF_NonExplore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(ws.Cells2Include,1)
            tt=['_t' num2str(ws.Cells2Include(i,1)) '-' num2str(ws.Cells2Include(i,2))];
            %
            ws.Crit{c,2}{2}{v}=['cF_Accept' tt];
            v=v+1;
            %
            ws.Crit{c,2}{2}{v}=['cF_Reject' tt];
            v=v+1;
        end
        c=c+1;

        % [Contrast 3] Control task, Explore ################
        ws.Crit{c,1}='ct_Explore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(ws.Cells2Include,1)
            tt=['_t' num2str(ws.Cells2Include(i,1)) '-' num2str(ws.Cells2Include(i,2))];
            %
            ws.Crit{c,2}{2}{v}=['ct_Explore' tt];
            v=v+1;
        end
        c=c+1;
        
        % [Contrast 4] Control task, NonExplore ################
        ws.Crit{c,1}='ct_NonExplore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(ws.Cells2Include,1)
            tt=['_t' num2str(ws.Cells2Include(i,1)) '-' num2str(ws.Cells2Include(i,2))];
            %
            ws.Crit{c,2}{2}{v}=['ct_NoBomb' tt];
            v=v+1;
            %
            ws.Crit{c,2}{2}{v}=['ct_Bomb' tt];
            v=v+1;
        end
        c=c+1;
    else
        ws.Crit=[];
    end
    
    % Save subject details
    SubjectsOK{s,5}=ws.Crit;
    ws=[];
end

% Alter list to include subjects who pass the criteria (log.subjects, log.n_subjs, SearchCriterias, RegLists)
log.old_subjects=log.subjects;
log.subjects=SubjectsOK(find(cell2mat(SubjectsOK(:,4))==1),1); log.n_subjs=length(log.subjects);
SearchCriterias=SubjectsOK(find(cell2mat(SubjectsOK(:,4))==1),5);
RegLists=RegLists(find(cell2mat(SubjectsOK(:,4))==1));

% Display for confirnation
disp('Criteria for subject inclusion:' ); disp(settings)
disp(' '); disp(['Removing bad  subjects  --------------------']); disp(SubjectsOK(find(cell2mat(SubjectsOK(:,4))==0),1))
disp([num2str(log.n_subjs) ' subjects included:']); disp(log.subjects); input('OK?  ');

%% (3) Get regressor weights for First-level contrast

if  log.execute_contrasts==1
    disp('####### Applying search criteria to derive regressor weights for contrasts  ##########'); disp(' ');
    ConInstruc=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##################'])
        [ConInstruc{s}] = xfx_2GetRegweightsFromSearchCriteria(SearchCriterias{s}, RegLists{s});
    end
end

%% (4) Generate batch for First-level contrasts

if  log.execute_contrasts==1
    
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep  'm_' log.firstlevel_contraststype ' Contrasted' filesep];
        [matlabbatch] = xfx_3Setup1stlevelContrastbatch(where,log, ws, s,ConInstruc{s});
    end
    
    % Record/output details
    Con.RegLists=RegLists;
    Con.SearchCriterias=SearchCriterias;
    Con.ConInstruc=ConInstruc;
    Con.SubjectsOK=SubjectsOK;
    Con.matlabbatch=matlabbatch;
end

%% (5) Set up the 2nd level models (executed next) that go with this model
% Establish instructions for 2nd level analyses to run ################
% Col 1=Name of 2nd-level model, Col 2=Which analysis function, Col 3=Inputs to function
% Instructions are executed in main script to run all possible 2nd level models

% Fetch details from the contrasted first level model
disp(' ################ Prepping for 2nd level ##################')
disp('Fetching details from (already-contrasted) first level model - sampling 1st subject')
c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep 'm_' log.firstlevel_contraststype ' Contrasted' filesep 'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run? 
%   1: Task x Explore factorial
%   2: Add sensible contrasts to the 2x2 (Task x Explore) Factorial
Nmodels=1;
Analysis=cell(Nmodels,3);

% ################################################################

% [Model #1] Task x Explore ################
for o1=1:1 
    m=1; Analysis{m,1}='TaskxExplore';
    Analysis{m,2}='NFlexFactorial';
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
    
    % Specify factorial design
    Analysis{m,3}.NFactorialCells=4;
    Analysis{m,3}.Design{1}.name='FacCell';      % Factor 1: NCells (2x2=4)
    Analysis{m,3}.Design{1}.dept=1;
    Analysis{m,3}.Design{1}.variance=1;
    Analysis{m,3}.Design{1}.gmsca=0;
    Analysis{m,3}.Design{1}.ancova=0;
    Analysis{m,3}.Design{2}.name='Subject';      % Factor 2: Subjects
    Analysis{m,3}.Design{2}.dept=0;
    Analysis{m,3}.Design{2}.variance=0;
    Analysis{m,3}.Design{2}.gmsca=0;
    Analysis{m,3}.Design{2}.ancova=0;
    
    % Factors to include in model 
    Analysis{m,3}.ModelFac=[1 2]; % Factor no.s (corresponding to above, i.e. factor 2=Subject)  
    
    % Assign contrast images to correct design cells (via search instruction)
    Analysis{m,3}.Cells=cell(2*2, 2); 
    Analysis{m,3}.Cells={'cF_Explore'             1 ; % Cell assignment can also be multi-factorial, e..g [1 2]
                                      'cF_NonExplore'       2;
                                      'ct_Explore'              3;
                                      'ct_NonExplore'        4;};
    Analysis{m,3}.Cells_FacAssignment=cell2mat(Analysis{m,3}.Cells(:,2));
    for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
        if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
            Analysis{m,3}.Cells{i,3}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
        else
            error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
        end
    end
    
    % Subject covariates (Exploration parameters from model fits)
    modelname='val2_b_fecx';
    explore_pars=nan*zeros(log.n_subjs,2);
    for s=1: log.n_subjs
        ws=SubjectsOK{s,2}.b.cF_modelpars.pars;
        %
        explore_pars(s,1)=ws{find(strcmp(ws(:,1), modelname)==1),7};
        explore_pars(s,2)=ws{find(strcmp(ws(:,1), modelname)==1),8};
        %
        ws=[];    
    end
    Analysis{m,3}.Covar(1).Type='single'; % 1st exploration parameter
    Analysis{m,3}.Covar(1).Ncells=4; 
    Analysis{m,3}.Covar(1).CovSingle_Vector=explore_pars(:,1);
    Analysis{m,3}.Covar(1).CovSingle_Name='Explore1';
    Analysis{m,3}.Covar(2).Type='single'; % 2nd exploration parameter
    Analysis{m,3}.Covar(2).Ncells=4; 
    Analysis{m,3}.Covar(2).CovSingle_Vector=explore_pars(:,2);
    Analysis{m,3}.Covar(2).CovSingle_Name='Explore2';
   
    % #############################################
    for o2=1:1 % [DISUSED] Using the [Non-flexible] Full factorial
%     Analysis{m,3}.Design{1}.name='Task';      % Factor 1: Task
%     Analysis{m,3}.Design{1}.levels=2;
%     Analysis{m,3}.Design{1}.dept=0;
%     Analysis{m,3}.Design{1}.variance=1;
%     Analysis{m,3}.Design{1}.gmsca=0;
%     Analysis{m,3}.Design{1}.ancova=0;
%     Analysis{m,3}.Design{2}.name='Explore';      % Factor 2: Explore
%     Analysis{m,3}.Design{2}.levels=2;
%     Analysis{m,3}.Design{2}.dept=0;
%     Analysis{m,3}.Design{2}.variance=1;
%     Analysis{m,3}.Design{2}.gmsca=0;
%     Analysis{m,3}.Design{2}.ancova=0;
%     
%     % Assign contrast images to correct design cells (via search instruction)
%     Analysis{m,3}.Cells=cell(2*2, 3); 
%     Analysis{m,3}.Cells={'cF_Explore'             1 1;
%                                       'cF_NonExplore'       1 2;
%                                       'ct_Explore'              2 1;
%                                       'ct_NonExplore'        2 2;};
%     for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
%         if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
%             Analysis{m,3}.Cells{i,1}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
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
%     Analysis{m,3}.Covar(1).Type='single';
%     Analysis{m,3}.Covar(1).Ncells=4; % 1st exploration parameter
%     Analysis{m,3}.Covar(1).CoVector=explore_pars(:,1);
%     Analysis{m,3}.Covar(1).CovNames=cellstr([char(repmat({'explore1_'}, log.n_subjs,1)) char(log.subjects)]);    
%     Analysis{m,3}.Covar(2).Ncells=4; % 2nd exploration parameter
%     Analysis{m,3}.Covar(2).CoVector=explore_pars(:,2);
%     Analysis{m,3}.Covar(2).CovNames=cellstr([char(repmat({'explore2_'}, log.n_subjs,1)) char(log.subjects)]);    
    
    end
end

% [Model #2] Add sensible contrasts to the Task x Explore factorial
for o1=1:1 
    m=m+1; Analysis{m,1}='Adding Factorial contrasts to TaskxExplore';
    Analysis{m,2}='AddContrasts';
    %
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep 'TaskxExplore'];
    Analysis{m,3}.DeleteOldCon=1; % There aren't any existing contrasts to begin with, using FlexFactorial
      
    % Weights for subject covariates
    %       covariate weights= ones(1, log.n_subjs)*(no. of regressors weighted)/(log.n_subjs)
    %       e.g. [1 1 0 0 ones(1,13)*(2)/13]
    SubCovarWeights=ones(1,log.n_subjs)/(log.n_subjs);
    SubCovarsForExplorePar=[0 0];
    
    
    % Contrast instructions
    Analysis{m,3}.ConInstruc={'ME_Task'                                 [1 1 -1 -1]         2;                                  % Factorial F contrasts
                                                'ME_ExploreOr'                          [1 -1 1 -1]         2;
                                                'Interaction_TaskxExpploreOr'    [1 -1 -1 1]         2;
                                                
    
                                                % T-contrasts
                                                'ConflictTask'                [1 1 0 0   SubCovarsForExplorePar      SubCovarWeights*2]           1;     % Main factor means
                                                'ControlTask'                [0 0 1 1   SubCovarsForExplorePar      SubCovarWeights*2]           1;
                                                'Explore'                       [1 0 1 0   SubCovarsForExplorePar      SubCovarWeights*2]           1;
                                                'NonExplore'                [0 1 0 1   SubCovarsForExplorePar      SubCovarWeights*2]           1;
                                                'cF_Explore'                [1 0 0 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;        % Cell means
                                                'ct_Explore'                [0 0 1 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;
                                                'cF_NonExplore'           [0 1 0 0   SubCovarsForExplorePar      SubCovarWeights*1]           1;
                                                'ct_NonExplore'             [0 0 0 1   SubCovarsForExplorePar      SubCovarWeights*1]           1;
                                                
                                                % Contrast comparisons (t-test) do not incldue subject-covariate weightings
                                                'Task_cF-ct'             [1 1 -1 -1]           1;        % Factor effects
                                                'Task_ct-cF'             [-1 -1 1 1]           1;
                                                'Explore-Non'             [1 -1 1 -1]           1;
                                                'Non-Explore'             [-1 1 -1 1]           1;
                                                'cF_Explore-Non'        [1 -1 0 0]           1;   % Cell effects
                                                'cF_Non-Explore'        [-1 1 0 0]           1;
                                                'ct_Explore-Non'        [0 0 1 -1]           1;
                                                'ct_Non-Explore'        [0 0 -1 1]           1;
                                                'cF_Explore_vs'         [3 -1 -1 -1]           1;   % Cell versus
                                                'cF_NonExplore_vs'    [-1 3 -1 -1]           1;
                                                'ct_Explore_vs'          [-1 -1 3 -1]           1;
                                                'ct_NonExplore_vs'     [-1 -1 -1 3]           1;
                                                'TaskxExp_Inter1'       [1 -1 -1 1]           1; % Interaction terms
                                                'TaskxExp_Inter2'       [-1 1 1 -1]           1;};
                                            
                                            
end

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

end

