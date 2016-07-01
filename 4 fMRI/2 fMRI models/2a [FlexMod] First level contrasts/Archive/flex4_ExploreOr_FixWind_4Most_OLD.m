function [Con log] = flex4_ExploreOr_FixWind_4Most(where, log)
% [Con log] = flex4_ExploreOr_FixWind_4Most(where, log)
% Compare Explore vs NonExplore, limited only to 4 most expected cells 
% EnvThreat 3 & 4, NTokens 10 & 12).
% Exclude subjects who have <30% Explore overall in these cells (not per
% cell)
%
% Second level models included:
%   (1) Task x Explore (2x2 Factorial)
%   (2) Add sensible contrasts to the 2x2 (Task x Explore) Factorial
%
% Output 'Con.ConInstruc' follows specifications of RegWeights output from the
% function f_2GetRegweightsFromSearchCriteria, details the (subject-specific) 
% regerssor weights for each contrast for each subject
%
% ---------------------------------------------------------------------------------------

% Check: Is this the correct Contrasts table for the chosen onsets model?
onsetsmodels4thiscontrasttype={'m_f1_ChoicexTrialtype'};
if sum(strcmp(log.onsetsmodel,onsetsmodels4thiscontrasttype))~=1
    error('Wrong contrasts type chosen for this onsets model!')
end

%% (1) General changeable settings

settings.min_percent_events=0.2;

%% (2) Generate subject-specific search criteria (First level contrasts)
% Search criteria is fed into function f_2GetRegweightsFromSearchCriteria,
% to derive regressor weights for each contrast. Format of search criteria
% must therefore follow format specified in this search function.

disp('########### Identifying included subjects + specifying subject-specific search criteria (to find regressors for weighting) #####################')
SubjectsOK=cell(log.n_subjs,3); % Col 2=% Explore Conflict,  Col 3=% Explore Control, Col 4=Subject ok?, Col 5= Search criteria
RegLists=cell(log.n_subjs,1); 
Cell2Include={3 5; 3 6; 4 5; 4 6}; % Col 1=Env (in REVERSE trial-type meaning. Altered when constructing search criteria), Col 2=NTokens

for s=1:log.n_subjs % Check subject stats + compile search criteria
    disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------'])
    [RegLists{s}] = f_1LoadSimpleReglist(where, log,s); % Load (edited) list of regressors
    ws.b=load([where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.subjects{s} '_behstats.mat']);
    
    % Counters
    ws.Conf_n_Exp=0;
    ws.Conf_n_Or=0;
    ws.Ctrl_n_Exp=0;
    ws.Ctrl_n_Or=0;
    for i=1:size(Cell2Include,1)
        % NOTE: EnvThreat is NOT in the same format as Trialtype 
        % (e.g. cF_Accept_t1-2 is the LOWEST EnvThreat. But in beh_stats, it
        % is in visualization space (i.e. cell 1,1=Highest thread, 2 Tokens)
        
        % Conflict task
        wc.conf=ws.b.cF_stats{Cell2Include{i,1},Cell2Include{i,2}};
        ws.Conf_n_Exp=ws.Conf_n_Exp+wc.conf.n_Explore;
        ws.Conf_n_Or=ws.Conf_n_Or+wc.conf.n_Accept+wc.conf.n_Reject;
        % Control task
        wc.ctrl=ws.b.ct_stats{Cell2Include{i,1},Cell2Include{i,2}};
        ws.Ctrl_n_Exp=ws.Ctrl_n_Exp+wc.ctrl.n_Explore;
        ws.Ctrl_n_Or=ws.Ctrl_n_Or+wc.ctrl.n_NoBomb+wc.ctrl.n_Bomb;
        
        %
        wc=[];
    end
    
    % Is this subject ok?
    SubjectsOK{s,1}=log.subjects{s};
    SubjectsOK{s,2}=ws.Conf_n_Exp/(ws.Conf_n_Exp+ws.Conf_n_Or);
    SubjectsOK{s,3}=ws.Ctrl_n_Exp/(ws.Ctrl_n_Exp+ws.Ctrl_n_Or);
    ws.ConflictOK=SubjectsOK{s,2}>=settings.min_percent_events && SubjectsOK{s,2}<=(1-settings.min_percent_events);
    ws.ControlOK=SubjectsOK{s,3}>=settings.min_percent_events && SubjectsOK{s,3}<=(1-settings.min_percent_events);
    if ws.ConflictOK==1 % && ws.ControlOK==1 % Include Conflict & Control as criteria
        SubjectsOK{s,4}=1;
    else
        SubjectsOK{s,4}=0;
    end
    
    % ####### Format for the search function, if this subject will be included ###########
    if  log.execute_contrasts==1 && SubjectsOK{s,4}==1
        ws.Crit=cell(4,2); c=1;
        % NOTE: So far, EnvThreat has NOT been in the same format as Trialtype
        % e.g. cF_Accept_t1-2 is the LOWEST EnvThreat. But in beh_stats, it
        % is in visualization space (i.e. cell 1,1=Highest thread, 2
        % Tokens). This need to be corrected before the Trial-type search
        % criteria is constructed, HERE. 

        % [Contrast 1] Conflict task, Explore ################
        ws.Crit{c,1}='cF_Explore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(Cell2Include,1)
            tt=['_t' num2str(7-Cell2Include{i,1}) '-' num2str(Cell2Include{i,2})];
            %
            ws.Crit{c,2}{2}{v}=['cF_Explore' tt];
            v=v+1;
        end
        c=c+1;

        % [Contrast 2] Conflict task, NonExplore ################
        ws.Crit{c,1}='cF_NonExplore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value    
        for i=1:size(Cell2Include,1)
            tt=['_t' num2str(7-Cell2Include{i,1}) '-' num2str(Cell2Include{i,2})];
            %
            ws.Crit{c,2}{2}{v}=['cF_Accept' tt];
            v=v+1;
            ws.Crit{c,2}{2}{v}=['cF_Reject' tt];
            v=v+1;
        end
        c=c+1;

        % [Contrast 3] Control task, Explore ################
        ws.Crit{c,1}='ct_Explore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value
        for i=1:size(Cell2Include,1)
            tt=['_t' num2str(7-Cell2Include{i,1}) '-' num2str(Cell2Include{i,2})];
            %
            ws.Crit{c,2}{2}{v}=['ct_Explore' tt];
            v=v+1;
        end
        c=c+1;

        % [Contrast 4] Control task, NonExplore ################
        ws.Crit{c,1}='ct_NonExplore';
        ws.Crit{c,2}{1}=1; v=1; % Weight value    
        for i=1:size(Cell2Include,1)
            tt=['_t' num2str(7-Cell2Include{i,1}) '-' num2str(Cell2Include{i,2})];
            %
            ws.Crit{c,2}{2}{v}=['ct_NoBomb' tt];
            v=v+1;
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

% Alter list to include stochastic subjects only (log.subjects, log.n_subjs, SearchCriterias, RegLists)
log.old_subjects=log.subjects;
log.subjects=SubjectsOK(find(cell2mat(SubjectsOK(:,4))==1),1); log.n_subjs=length(log.subjects);
SearchCriterias=SubjectsOK(find(cell2mat(SubjectsOK(:,4))==1),5);
RegLists=RegLists(find(cell2mat(SubjectsOK(:,4))==1));

% Display for confirnation
disp(' '); disp('Subject-inclusion is based on CONFLICT task performance only, not Control task ##############'); 
disp(' '); disp('Removing deterministic subjects  --------------------'); disp(SubjectsOK(find(cell2mat(SubjectsOK(:,4))==0),1))
disp([num2str(log.n_subjs) ' subjects included (criteria for inclusion: >' num2str(settings.min_percent_events) '): ']);
disp(log.subjects); input('OK?  ');

%% (3) Get regressor weights for First-level contrast

if log.execute_contrasts==1
    disp('####### Applying search criteria to derive regressor weights for contrasts  ##########'); disp(' ');
    ConInstruc=cell(log.n_subjs,1);
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##################'])
        [ConInstruc{s}] = f_2GetRegweightsFromSearchCriteria(SearchCriterias{s}, RegLists{s});
    end
end

%% (4) Generate batch for First-level contrasts

if  log.execute_contrasts==1
    
    input('Continue to execute contrasts?    ')
    disp('############### Executing contrasts #####################')
    
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '   -  '  log.subjects{s} ' ##########'])
        ws.where_contrasts=[where.data_brain filesep log.subjects{s} filesep '2 First level' filesep log.onsetsmodel ' Contrasted ' log.firstlevel_contraststype filesep];
        [matlabbatch] = f_3Setup1stlevelContrastbatch(where,log, ws, s,ConInstruc{s});
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
c=load([where.data_brain filesep log.subjects{1} filesep '2 First level' filesep log.onsetsmodel ' Contrasted ' log.firstlevel_contraststype filesep 'SPM.mat']);
for i=1:size(c.SPM.xCon,2); Con.ConNames{i,1}=c.SPM.xCon(i).name; Con.ConNames{i,2}=c.SPM.xCon(i).Vcon.fname; end
disp(Con.ConNames);

% Models to run? 
%   1: Task x Explore factorial
%   2: Add sensible contrasts to the 2x2 (Task x Explore) Factorial
Nmodels=2;
Analysis=cell(Nmodels,3);

% ################################################################

% [Model #1] Task x Explore ################
for o1=1:1 
    m=1; Analysis{m,1}='TaskxExplore';
    Analysis{m,2}='NFactorial';
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep Analysis{m,1}];
    %
    Analysis{m,3}.Design{1}.name='Task';      % Factor 1: Task
    Analysis{m,3}.Design{1}.levels=2;
    Analysis{m,3}.Design{1}.dept=0;
    Analysis{m,3}.Design{1}.variance=1;
    Analysis{m,3}.Design{1}.gmsca=0;
    Analysis{m,3}.Design{1}.ancova=0;
    Analysis{m,3}.Design{2}.name='Explore';      % Factor 2: Explore
    Analysis{m,3}.Design{2}.levels=2;
    Analysis{m,3}.Design{2}.dept=0;
    Analysis{m,3}.Design{2}.variance=1;
    Analysis{m,3}.Design{2}.gmsca=0;
    Analysis{m,3}.Design{2}.ancova=0;
    
    % Assign contrast images to correct design cells (via search instruction)
    Analysis{m,3}.Cells=cell(2*2, 3); 
    Analysis{m,3}.Cells={'cF_Explore'             1 1;
                                      'cF_NonExplore'       1 2;
                                      'ct_Explore'              2 1;
                                      'ct_NonExplore'        2 2;};
    for i=1:size(Analysis{m,3}.Cells,1) % Replace cell names with contrast images
        if length(find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})))==1
            Analysis{m,3}.Cells{i,1}=Con.ConNames{find(strcmp(Con.ConNames(:,1),Analysis{m,3}.Cells{i,1})),2};
        else
            error('When fetching contrast images for each cell, # of con-image matches for a particular cell ~=1')
        end
    end
    
    % Subject covariates
    Analysis{m,3}.Covar(1).Ncells=4;
    Analysis{m,3}.Covar(1).CoVector=ones(log.n_subjs,1);
    Analysis{m,3}.Covar(1).CovNames=cellstr([char(repmat({'sub_'}, log.n_subjs,1)) char(log.subjects)]);    
end

% [Model #2] Add sensible contrasts to the Task x Explore factorial
for o1=1:1 
    m=m+1; Analysis{m,1}='Adding Factorial contrasts to TaskxExplore';
    Analysis{m,2}='AddContrasts';
    %
    Analysis{m,3}.Directory=[where.firstlevel_resultsfolder filesep 'TaskxExplore'];
    Analysis{m,3}.DeleteOldCon=0;
    SubCovarWeights=ones(1,log.n_subjs)/(log.n_subjs/2);
    Analysis{m,3}.ConInstruc={'ConflictTask'                [1 1 0 0   SubCovarWeights];     % Main factor means
                                                'ControlTask'                [0 0 1 1   SubCovarWeights];
                                                'Explore'                       [1 0 1 0   SubCovarWeights];
                                                'NonExplore'                [0 1 0 1   SubCovarWeights];
                                                'cF_Explore'                [1 0 0 0   SubCovarWeights];        % Cell means
                                                'ct_Explore'                [0 0 1 0   SubCovarWeights];
                                                'cF_NonExplore'           [0 1 0 0   SubCovarWeights];
                                                'ct_NonExplore'             [0 0 0 1   SubCovarWeights];
                                                'Task_cF-ct'             [1 1 -1 -1   SubCovarWeights];        % Factor effects
                                                'Task_ct-cF'             [-1 -1 1 1   SubCovarWeights];
                                                'Explore-Non'             [1 -1 1 -1   SubCovarWeights];
                                                'Non-Explore'             [-1 1 -1 1   SubCovarWeights];
                                                'cF_Explore-Non'        [1 -1 0 0   SubCovarWeights];   % Cell effects
                                                'cF_Non-Explore'        [-1 1 0 0   SubCovarWeights];
                                                'ct_Explore-Non'        [0 0 1 -1   SubCovarWeights];
                                                'ct_Non-Explore'        [0 0 -1 1   SubCovarWeights];
                                                'cF_Explore_vs'         [3 -1 -1 -1   SubCovarWeights];   % Cell versus
                                                'cF_NonExplore_vs'    [-1 3 -1 -1   SubCovarWeights];
                                                'ct_Explore_vs'          [-1 -1 3 -1   SubCovarWeights];
                                                'ct_NonExplore_vs'     [-1 -1 -1 3   SubCovarWeights];
                                                'TaskxExp_Inter1'       [1 -1 -1 1   SubCovarWeights]; % Interaction terms
                                                'TaskxExp_Inter2'       [-1 1 1 -1   SubCovarWeights];};
                                            
end

% ################################################################

% Output to main script
Con.SecondLevelAnalysis=Analysis;

end

